"""Independent full-energy, analytic remainder and coupled comparison diagnostics."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_hybrid_core_regulator_comparison import (
    audit,
    domain,
    dynamics,
    operators,
    quadratic,
    source,
)
from p8_vacuum_affine_selfconsistent_finite_feedback import symmetry
from p8_vacuum_affine_weyl_operator_comparison import (
    derivatives as original_derivatives,
)

pytestmark = pytest.mark.filterwarnings("error::RuntimeWarning")
PACKETS = (source.data, quadratic.data, domain.data, operators.data, dynamics.data)


@pytest.mark.parametrize("packet", PACKETS)
def test_all_entire_exact_packet_residuals_and_gates(packet):
    data = packet()
    for name, value in data["checks"].items():
        entries = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.factor(entry) == 0 for entry in entries), name
    assert all(bool(value) for value in data["gates"].values())


@cache
def independent_radial_partitions():
    # Partial Bell recurrence from singleton weight4 and pair weight2.
    # This is not the explicit factorial partition sum in the implementation.
    bell = [[s.Integer(1)]]
    sums = [s.Integer(1)]
    for n in range(1, 207):
        row = [s.Integer(0)] * (n + 1)
        for k in range(1, n + 1):
            first = bell[n - 1][k - 1]
            second = bell[n - 2][k - 1] if n >= 2 and k - 1 < len(bell[n - 2]) else 0
            row[k] = 4 * first + 2 * (n - 1) * second
        bell.append(row)
        sums.append(sum(row[k] * 256**k * s.factorial(k) ** 2 for k in range(n + 1)))
    return sums


@pytest.mark.parametrize("order", range(207))
def test_every_radial_partition_by_independent_Bell_recurrence(order):
    literal = independent_radial_partitions()[order]
    assert literal == operators.radial_partition(order)
    assert literal <= 2048**order * s.factorial(order) ** 3


@pytest.mark.parametrize("order", range(206))
def test_every_phase_ratio_and_all_parameter_derivatives(order):
    for amplitude in (source.H_AMPLITUDE, source.F_AMPLITUDE):
        for j in range(3):
            ratio = operators.jet(order + 1, amplitude, j) / operators.jet(
                order, amplitude, j
            )
            assert ratio == 2056 * (order + 1) ** 3 / source.RADIUS
            assert ratio < s.Rational(1, 10**9)


@pytest.mark.parametrize("order", range(197))
def test_all_original_phase_product_coefficients_retained(order):
    assert operators.coefficient(order) == original_derivatives.coefficient(order)


@pytest.mark.parametrize(
    "normal",
    ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, -1, 2), (2, 3, 5), (-3, 2, 4)),
)
@pytest.mark.parametrize("seed_id", range(3))
def test_full_constraint_and_DeWitt_cancellation_for_independent_modes(normal, seed_id):
    k = s.Matrix(normal)
    radial = k * k.T / k.dot(k)
    transverse = s.eye(3) - radial
    A = s.Matrix([[2 + seed_id, 3, -1], [3, -2, 4 - seed_id], [-1, 4 - seed_id, 1]])
    TT = transverse * A * transverse - transverse * s.trace(transverse * A) / 2
    r, S = s.symbols("independent_trace independent_source")
    longitudinal = (3 * S - r) * (radial - s.eye(3) / 3) / 4
    full = r * s.eye(3) / 6 + TT + longitudinal
    assert (full * k - S * k / 2).applyfunc(s.expand) == s.zeros(3, 1)
    assert s.simplify(s.trace(TT)) == 0 and TT * k == s.zeros(3, 1)
    assert s.expand(s.trace(TT * longitudinal)) == 0
    density = 2 * s.trace(full * full) - s.trace(full) ** 2
    assert s.expand(density - 2 * s.trace(TT * TT) + r * S / 2 - 3 * S * S / 4) == 0
    assert s.diff(s.expand(density), r, 2) == 0
    # Deleting either of the cancelling scalar momentum pieces is wrong.
    assert s.diff(-r * r / 12, r, 2) != 0
    assert s.diff(2 * s.trace(longitudinal * longitudinal), r, 2) != 0


def test_literal_H2_independently_by_multivariate_chain_rule():
    packet = quadratic.literal_hessian()
    rules = packet["whole_density_substitution"]
    eps = next(
        atom
        for value in rules.values()
        for atom in value.free_symbols
        if atom.name == "independent_variation_parameter"
    )
    v = next(
        atom
        for value in rules.values()
        for atom in value.free_symbols
        if atom.name == "volume_field"
    )
    zero = {z: value.subs(eps, 0) for z, value in rules.items()}
    first = {z: s.diff(value, eps).subs(eps, 0) for z, value in rules.items()}
    second = {z: s.diff(value, eps, 2).subs(eps, 0) for z, value in rules.items()}
    H = source.q.HAMILTONIAN
    independent = s.Rational(9, 2) * v * v * H.subs(zero, simultaneous=True)
    for z in source.q.COORDS:
        Hz = s.diff(H, z).subs(zero, simultaneous=True)
        independent += Hz * (second[z] / 2 + 3 * v * first[z])
        for w in source.q.COORDS:
            independent += (
                s.diff(H, z, w).subs(zero, simultaneous=True) * first[z] * first[w] / 2
            )
    assert s.expand(independent - packet["whole_original_H2"]) == 0


def test_entire_original_constraint_is_full_lapse_stationarity_equation():
    q = source.q
    assert (
        s.expand(q.eliminate_N_primitive(s.diff(q.HAMILTONIAN, q.N) - q.CONSTRAINT))
        == 0
    )
    C1 = quadratic.literal_hessian()["whole_literal_C1"]
    assert C1.has(q.Hclock)
    assert not s.expand(quadratic.literal_hessian()["whole_original_H2"]).has(q.Hclock)


@pytest.mark.parametrize(
    "CN,C1,H2",
    ((3, 7, 11), (-4, 5, -2), (s.Rational(5, 2), s.Rational(2, 3), s.Rational(7, 5))),
)
def test_complete_implicit_quadratic_schur_correction_sign(CN, C1, H2):
    n, z = s.symbols("independent_lapse_shift independent_phase")
    H = s.Rational(CN, 2) * n * n + C1 * n * z + H2 * z * z
    root = s.solve(s.diff(H, n), n)[0]
    on_shell = s.expand(H.subs(n, root))
    assert on_shell == (H2 - s.Rational(C1 * C1, 2 * CN)) * z * z
    assert abs(s.Rational(C1 * C1, 2 * CN)) <= s.Rational(C1 * C1, 4)
    assert s.diff(on_shell - H.subs(n, 0), z, 2) != 0


@cache
def independently_integrated_diagonal_EH():
    # An exact metric, not a linearized Ricci formula imported from the module.
    x, eps, v, tau = s.symbols("x eps v tau", real=True)
    profiles = (v, v - tau / 2, v + tau / 2)
    g = s.diag(*(s.exp(2 * eps * a * s.cos(x)) for a in profiles))
    inv = g.inv()

    def d(value, index):
        return s.diff(value, x) if index == 0 else s.Integer(0)

    Gamma = [
        [
            [
                s.simplify(
                    sum(
                        inv[a, b] * (d(g[b, j], i) + d(g[b, i], j) - d(g[i, j], b))
                        for b in range(3)
                    )
                    / 2
                )
                for j in range(3)
            ]
            for i in range(3)
        ]
        for a in range(3)
    ]
    Ricci = s.Matrix(
        3,
        3,
        lambda i, j: sum(
            d(Gamma[a][i][j], a)
            - d(Gamma[a][i][a], j)
            + sum(
                Gamma[a][a][b] * Gamma[b][i][j] - Gamma[a][j][b] * Gamma[b][i][a]
                for b in range(3)
            )
            for a in range(3)
        ),
    )
    curvature = s.simplify(s.trace(inv * Ricci))
    density = s.exp(3 * eps * v * s.cos(x))
    EH2 = s.simplify(s.diff(density * curvature, eps, 2).subs(eps, 0) / 2)
    R2 = s.simplify(s.diff(curvature, eps, 2).subs(eps, 0) / 2)
    R1 = s.simplify(s.diff(curvature, eps).subs(eps, 0))
    average = lambda value: (
        s.integrate(s.expand_trig(value), (x, 0, 2 * s.pi)) / (2 * s.pi)
    )
    return v, tau, average(EH2), average(R2), average(3 * v * s.cos(x) * R1)


def test_complete_exact_metric_EH_keeps_tensor_and_density_contacts():
    v, tau, total, curvature, density_contact = independently_integrated_diagonal_EH()
    assert s.expand(total - v * v + tau * tau / 4) == 0
    assert s.expand(curvature + 5 * v * v + tau * tau / 4) == 0
    assert s.expand(density_contact - 6 * v * v) == 0
    assert s.expand(total - curvature - density_contact) == 0
    assert s.diff(total, tau, 2) != 0


@pytest.mark.parametrize("normal", ((1, 0, 0), (1, 2, 2), (2, -1, 2)))
def test_all_Proca_polarizations_and_directional_divergence_curl(normal):
    k = s.Matrix(normal)
    L = k * k.T / k.dot(k)
    PT = s.eye(3) - L
    mass, omega, zeta, kap = s.symbols("mass omega zeta kap", positive=True)
    B = PT + (omega / mass) * L
    B_inverse = PT + (mass / omega) * L
    Q = B / s.sqrt(omega * kap * zeta)
    P = B_inverse * s.sqrt(omega * zeta / kap)
    assert (Q * P.T - s.eye(3) / kap).applyfunc(s.simplify) == s.zeros(3)
    assert k.T * PT == s.zeros(1, 3)
    assert (k.T * P - mass * s.sqrt(zeta / (omega * kap)) * k.T).applyfunc(
        s.simplify
    ) == s.zeros(1, 3)
    cross = s.Matrix([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    assert cross * L == s.zeros(3)
    assert (cross * Q - cross / s.sqrt(omega * kap * zeta)).applyfunc(
        s.simplify
    ) == s.zeros(3)


@pytest.mark.parametrize("radius", (source.RADIUS, source.ANALYSIS_RADIUS))
def test_all_full_phase_images_fit_common_unmodified_source_domain(radius):
    b = domain.bounds(radius)
    assert len(b["images"]) == len(b["averages"]) == 12
    assert b["D0"] == sum(
        b[key] for key in ("full_Dtrace", "full_Dshape", "full_Dvector", "full_Dmatter")
    )
    assert b["whole_center_residual"] > 12 * 10**4 * source.previous.HOMOGENEOUS_RADIUS
    assert (
        source.previous.HOMOGENEOUS_RADIUS + max(b["images"]) + 8 * source.TIME
        < source.previous.COORDINATE
    )
    assert (
        b["whole_center_residual"] / 3 + b["contraction"] * source.previous.LAPSE_RADIUS
        < source.previous.LAPSE_RADIUS
    )


@pytest.mark.parametrize(
    "point", (s.Rational(1, 10), s.Rational(1, 4), s.Rational(2, 5))
)
def test_whole_analytic_remainder_is_bounded_not_deleted(point):
    z = s.Symbol("independent_radial_variable")
    f = z * z / (1 - z / 2)
    quadratic = s.diff(f, z, 2).subs(z, 0) * z * z / 2
    actual = s.factor(f - quadratic).subs(z, point)
    # On complex |z|<=1, |f|<=2. This bounds the WHOLE tail on point<1.
    cauchy = 2 * point**3 / (1 - point)
    assert 0 < actual < cauchy
    assert s.diff(f, z, 3).subs(z, 0) != 0
    assert (
        s.factor(f - quadratic - actual) != 0
    )  # A numerical tail is not a new symbol.


@pytest.mark.parametrize("order", range(7))
def test_exact_heat_identity_by_coefficient_telescoping(order):
    t, D = s.symbols("independent_heat_time independent_generator")
    P = sum((-t * D) ** j / s.factorial(j) for j in range(order + 1))
    # Remove the common exponential and compare every polynomial coefficient.
    derivative = s.Poly(s.diff(P, t) + D * P, t)
    assert derivative.degree() == order
    for power in range(order + 1):
        target = (
            (-1) ** order * D ** (order + 1) / s.factorial(order)
            if power == order
            else 0
        )
        assert s.expand(derivative.nth(power) - target) == 0


@pytest.mark.parametrize("degree", (2, 4, 8, 12, 14, 16, 18))
def test_full_positive_heat_integral_remainder_for_higher_polynomials(degree):
    x, t = s.symbols("independent_x independent_theta")
    polynomial = x**degree

    def D(value):
        return s.diff(value, x, 2) / 4

    def heat(value, time):
        current = value
        result = value
        for k in range(1, degree // 2 + 1):
            current = D(current)
            result += time**k * current / s.factorial(k)
        return s.expand(result)

    powers = [polynomial]
    for _ in range(7):
        powers.append(D(powers[-1]))
    P6 = sum((-1) ** j * powers[j] / s.factorial(j) for j in range(7))
    remainder = -s.integrate(t**6 * heat(powers[7], t), (t, 0, 1)) / s.factorial(6)
    assert s.expand(polynomial - heat(P6, 1) - remainder) == 0
    if degree >= 14:
        assert remainder != 0


def test_compact_core_vanishing_does_not_pass_to_full_heat_remainder():
    y, delta = s.symbols("heat_integration_coordinate positive_small_displacement")
    weight = s.exp(-y * y)
    kernel = s.simplify(
        sum(
            (-1) ** j * s.diff(weight, y, 2 * j) / (4**j * s.factorial(j))
            for j in range(7)
        )
        / weight
    )
    p = s.Poly(s.expand(kernel.subs(y, 1 + delta)), delta)
    center = p.nth(0)
    variation = sum(
        abs(p.nth(j)) * s.Rational(1, 100) ** j for j in range(1, p.degree() + 1)
    )
    assert abs(center) > variation and center != 0
    # Any nonzero positive smooth bump supported in (1,1.01) vanishes on
    # the core (-1/2,1/2), as do all its derivatives. Its convolution with
    # this sign-definite heat*P6 kernel at0 is nonzero. Since b(0)=0,
    # r6(0)=-heat(P6 b)(0) is nonzero despite the entire compact symbol jet.
    assert "generally does NOT" in operators.data()["whole_cutoff_difference_warning"]


@pytest.mark.parametrize(
    "name,amplitude", (("H", source.H_AMPLITUDE), ("F", source.F_AMPLITUDE))
)
@pytest.mark.parametrize("order", range(3))
def test_both_operator_bounds_and_ordering_difference_retain_all_terms(
    name, amplitude, order
):
    b = operators.bounds(amplitude, order)
    terms = b["whole_all_eight_heat_terms"]
    assert len(terms) == 8 and all(term > 0 for term in terms)
    assert b["whole_finite_coherent_polynomial_bound"] == sum(terms[:7])
    assert b["whole_retained_Weyl_heat_remainder"] == 10**112 * terms[7] > 0
    assert (
        b["whole_both_complete_operator_bound"] == sum(terms[:7]) + 10**112 * terms[7]
    )
    assert (
        b["whole_Weyl_minus_calibrated_operator_difference"]
        == sum(terms[2:7]) + 10**112 * terms[7]
    )


@pytest.mark.parametrize(
    "a,b,c,e,f",
    (
        (
            s.Rational(1, 10),
            s.Rational(1, 5),
            s.Rational(1, 4),
            s.Rational(1, 7),
            s.Rational(1, 11),
        ),
        (
            s.Rational(1, 4),
            s.Rational(1, 10),
            s.Rational(1, 3),
            s.Rational(1, 13),
            s.Rational(1, 17),
        ),
    ),
)
def test_full_coupled_two_row_inverse_not_identical_background_histories(a, b, c, e, f):
    A = s.Matrix([[1 - a, -b], [-c, 1]])
    actual = A.inv() * s.Matrix([f, e])
    y = (b * e + f) / (1 - a - b * c)
    z = c * y + e
    assert actual == s.Matrix([y, z])
    assert y > f / (1 - a) and z > e
    assert s.det(A) == 1 - a - b * c


@pytest.mark.parametrize("angle", ("0.01", "0.2", "0.7"))
def test_positive_coherent_effect_probability_and_covariant_transport(angle):
    with mp.workdps(60):
        tail = mp.mpf("0.0001")
        seed = mp.matrix([mp.sqrt(1 - tail), mp.sqrt(tail)])
        theta = mp.mpf(angle)
        U = mp.matrix([[mp.cos(theta), -mp.sin(theta)], [mp.sin(theta), mp.cos(theta)]])
        effect = mp.matrix([[0, 0], [0, 1]])
        state = U * seed
        probability = mp.re((state.H * effect * state)[0])
        assert 0 < probability <= tail + 2 * mp.norm(state - seed)
        V = mp.matrix([[1, mp.j], [mp.j, 1]]) / mp.sqrt(2)
        physical = V * state
        transported = V * effect * V.H
        assert abs((physical.H * transported * physical)[0] - probability) < mp.mpf(
            "1e-55"
        )
        assert abs((physical.H * effect * physical)[0] - probability) > mp.mpf("1e-4")


@pytest.mark.parametrize("phase", (s.pi / 2, s.pi, 3 * s.pi / 2))
def test_distinct_global_scalar_phases_do_not_define_small_physical_vector_distance(
    phase,
):
    seed = s.Matrix([1, 0])
    left = seed
    right = s.exp(s.I * phase) * seed
    assert s.simplify(((left - right).H * (left - right))[0]) >= 2
    readout = s.Matrix([[3, 1], [1, 4]])
    assert (
        s.simplify((left.H * readout * left)[0] - (right.H * readout * right)[0]) == 0
    )
    assert "OWN" in dynamics.data()["whole_scalar_phase_and_covariance_boundary"]


@pytest.mark.parametrize("index", range(24))
def test_new_bounds_do_not_replace_actual_cubic_covariance_or_symmetry(index):
    U = symmetry.representations()[index]
    assert U * U.T == s.eye(48)
    assert source.RADIUS == symmetry.seed.CORE
    assert audit.MODEL == audit.previous.MODEL
    assert source.H_AMPLITUDE < audit.previous.quantum.H_AMPLITUDE


@pytest.mark.parametrize(
    "api,value",
    (
        (operators.radial_partition, True),
        (operators.radial_partition, 1.0),
        (operators.radial_partition, s.Float(1)),
        (domain.bounds, float(source.RADIUS)),
    ),
)
def test_guards_run_before_cache_lookup_after_valid_warmup(api, value):
    if api is domain.bounds:
        api(source.RADIUS)
    else:
        api(1)
    with pytest.raises((TypeError, ValueError)):
        api(value)


@pytest.mark.parametrize("bad", (True, False, 0.0, 1.0, s.Float(0), s.Float(1)))
def test_cached_operator_parameter_order_cannot_bypass_exact_guard(bad):
    operators.bounds(source.H_AMPLITUDE, 0)
    operators.bounds(source.H_AMPLITUDE, 1)
    with pytest.raises((TypeError, ValueError)):
        operators.bounds(source.H_AMPLITUDE, bad)


@pytest.mark.parametrize(
    "name,call,args",
    audit.bad_cases(),
    ids=lambda value: value if isinstance(value, str) else None,
)
def test_all_unsupported_inputs_and_scope_upgrades_are_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_complete_original_frontier_and_new_comparison_scope():
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 130
    assert len({row["id"] for row in audit.matching()}) == 130
    assert audit.ITEM["id"] == audit.matching()[-1]["id"]
    assert "OPEN" in audit.observable()["not_established"]
    assert "No exact support" in audit.observable()["not_established"]
