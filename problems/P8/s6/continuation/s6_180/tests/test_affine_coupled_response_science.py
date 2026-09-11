"""Independent literal tree, covariance, primitive and reconstruction checks."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_coupled_response import (
    audit,
    classical,
    inverse,
    matching,
    quantum,
)

PACKETS = list(audit.packets())


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_all_exact_residuals(name):
    value = audit.residuals()[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    ), name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_unsupported_inverse_scope_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("name", PACKETS)
def test_each_quantitative_packet(name):
    assert all(audit.packets()[name].get("gates", {}).values())


def test_actual_boundary_omission_changes_bounce_lapse_constraint():
    d = classical.clock_coefficients()
    J = classical.quadratic()["J"]
    missing = s.factor(classical.weighted(d["primitive_NN"]) / 2)
    assert missing.subs(classical.u, 0) == -9
    assert J.subs(classical.u, 0) == s.Rational(243, 160)
    assert (J + missing).subs(classical.u, 0) < 0
    assert J + missing != J


@pytest.mark.parametrize(
    "time",
    (-s.Rational(1, 2), -s.Rational(1, 4), 0, s.Rational(1, 4), s.Rational(1, 2)),
)
def test_literal_two_velocity_legendre_lapse_pivot(time):
    q = classical.quadratic()
    n, vd, sd, pv, ps = s.symbols("n vd sd pv ps", real=True)
    theta, J, w = [q[key].subs(classical.u, time) for key in ("theta", "J", "w")]
    L = (
        -3 * vd * vd
        + 6 * theta * n * vd
        + (J + w * w / 2 - 3 * theta * theta) * n * n
        + sd * sd / 2
        + w * n * sd
    )
    solution = s.solve([s.diff(L, vd) - pv, s.diff(L, sd) - ps], (vd, sd))
    ham = s.factor((pv * vd + ps * sd - L).subs(solution, simultaneous=True))
    assert s.factor(-s.diff(ham, n, 2) / 2 - J) == 0
    assert J > 0
    assert abs(classical.adapted()["A"].subs(classical.u, time)) >= s.Rational(
        6144, 15625
    )


@pytest.mark.parametrize("sector", ("T", "L"))
@pytest.mark.parametrize("D", (2, 3, 4, 5))
def test_literal_full_dimensional_constraint_and_current(sector, D):
    N, a, k, m, A0, A, V = s.symbols("N a k m A0 A V", positive=True)
    if sector == "T":
        L = (
            a ** (D - 2) * V * V / (2 * N)
            - N * (a ** (D - 4) * k * k + a ** (D - 2) * m * m) * A * A / 2
        )
    else:
        L = (
            a ** (D - 2) * (V - k * A0) ** 2 / (2 * N)
            + a**D * m * m * A0 * A0 / (2 * N)
            - N * a ** (D - 2) * m * m * A * A / 2
        )
        temporal = s.solve(s.diff(L, A0), A0)[0]
        L = s.factor(L.subs(A0, temporal))
    g2 = s.factor(s.diff(L, V, 2))
    frequency2 = s.factor(-s.diff(L, A, 2) / g2)
    current_pair = s.factor(
        a * s.diff(g2, a) / g2 + a * s.diff(g2 * frequency2, a) / (g2 * frequency2)
    )
    high = s.limit(current_pair, k, s.oo)
    assert high == 2 * (D - 3 if sector == "T" else D - 1)
    # Under t -> psi(t), N=psi' and A=A0(psi), the full
    # coordinate action is multiplied by psi', including the constraint.
    lapse = s.symbols("positive_time_jacobian", positive=True)
    changed = s.factor(L.subs({N: lapse, V: lapse * V}, simultaneous=True))
    assert s.factor(changed - lapse * L.subs(N, 1)) == 0


@pytest.mark.parametrize("c,b", ((1, 0), (0, 1), (2, -3)))
def test_covariant_cosmological_and_Einstein_density_uncancelled_Ward(c, b):
    u = quantum.u
    H = quantum.H
    eta, w = quantum.eta, quantum.w
    n = s.diff(eta, u)
    v = w + H * eta
    rho = -c - 6 * b * H * H
    P = c + 4 * b * s.diff(H, u) + 6 * b * H * H
    drho = -12 * b * H * (s.diff(v, u) - H * n)
    dP = 4 * b * (
        s.diff(v, u, 2) - H * s.diff(n, u) - 2 * s.diff(H, u) * n
    ) + 12 * b * H * (s.diff(v, u) - H * n)
    R = -12 * b * H * s.diff(w, u)
    Q = 4 * b * s.diff(w, u, 2) + 12 * b * H * s.diff(w, u)
    physical = s.Matrix([-drho - 3 * rho * v, 3 * dP + 3 * P * (n + 3 * v)])
    transformed = s.Matrix([-quantum.D(physical[0]) + H * physical[1], physical[1]])
    target = (
        quantum.data()["complete_adapted_vector_response"]
        .subs(
            {quantum.rho: rho, quantum.P: P, quantum.R: R, quantum.Q: Q},
            simultaneous=True,
        )
        .doit()
    )
    assert (transformed - target).applyfunc(s.simplify) == s.zeros(2, 1)
    if c == 1 and b == 0:
        assert quantum.data()["local_time_channel"].subs(quantum.P, 1) != 0


@pytest.mark.parametrize("degree", (4, 6, 9))
def test_variable_fourth_coefficient_full_primitive_commutator(degree):
    t, r = s.symbols("t r", positive=True)
    A = lambda z: 1 + z + 2 * z * z + z**5
    eta = lambda z: z**degree
    direct = s.integrate((t - r) ** 3 * A(r) * s.diff(eta(r), r, 4) / 6, (r, 0, t))
    remainder = s.diff((t - r) ** 3 * A(r) / 6, r, 4)
    exact = A(t) * eta(t) + s.integrate(remainder * eta(r), (r, 0, t))
    assert s.expand(direct - exact) == 0
    assert s.expand(direct - A(t) * eta(t)) != 0


@pytest.mark.parametrize("order", range(4))
def test_literal_lower_local_primitive_without_endpoint_loss(order):
    t, r = s.symbols("t r", positive=True)
    coeff = 1 + r + r**3
    z = r**7
    direct = s.integrate((t - r) ** 3 * coeff * s.diff(z, r, order) / 6, (r, 0, t))
    kernel = (-1) ** order * s.diff((t - r) ** 3 * coeff / 6, r, order)
    expected = s.integrate(kernel * z, (r, 0, t))
    assert s.expand(direct - expected) == 0
    assert s.diff(kernel, t).is_polynomial(t, r)


@pytest.mark.parametrize("col,order", [(i, j) for i in (0, 1) for j in range(4)])
def test_actual_first_row_continuous_derivative_envelopes(col, order):
    coef = classical.adapted()["complete_derivative_coefficients"]["0" + str(col)][
        order
    ]
    bounds = inverse.first_row()[
        "actual_classical_first_row_coefficient_derivative_envelopes"
    ][(col, order)]
    for j, upper in bounds.items():
        actual = s.diff(coef, classical.u, j)
        for point in (
            -s.Rational(1, 2),
            -s.Rational(1, 4),
            0,
            s.Rational(1, 4),
            s.Rational(1, 2),
        ):
            assert abs(actual.subs(classical.u, point)) <= upper


def test_explicit_local_lapse_recovery_not_total_inverse_norm():
    p = inverse.first_row()
    C1 = p["physical_lapse_reconstruction_C1_upper"]
    assert s.Rational(14976) < C1 < s.Rational(15000)
    assert (
        C1
        == p["actual_A_derivative_envelopes"][1]
        + p["kernel_diagonal_row_upper"]
        + p["kernel_first_output_derivative_integral_upper"]
    )
    assert p["quantum_first_row_added_bound"] > 0
    assert "not for the full inverse" in p["C1_bound_uses"]
    assert "not numerically evaluated" in inverse.data()["uncomputed"]


@pytest.mark.parametrize(
    "B", (s.Rational(1, 1000), s.S.One, s.Integer(1000), s.Integer(10) ** 800)
)
def test_constructive_weight_for_every_positive_majorant(B):
    weight = (4 * B + 1) ** 2
    contraction = 2 * B / (4 * B + 1)
    assert weight > 1
    assert 0 < contraction < s.Rational(1, 2)
    assert 1 / (1 - contraction) < 2
    assert s.sqrt(weight) == 4 * B + 1


@pytest.mark.parametrize("weight", (1, 2, 10, 1000))
def test_weighted_weak_log_integral_bound_independent_quadrature(weight):
    with mp.workdps(45):
        l = mp.mpf(weight)
        value = mp.quad(
            lambda r: mp.exp(-l * r) * (1 - mp.log(r)), [0, mp.mpf("0.01"), 1]
        )
        assert 0 < value <= (2 + mp.log(l)) / l
        assert (2 + mp.log(l)) / l <= 2 / mp.sqrt(l) + mp.mpf("1e-40")


@pytest.mark.parametrize("p", (0, 1, -1, -3, 1 + 2j, -5 + 1j))
def test_same_massive_range_gap_not_full_coupled_spectrum(p):
    with mp.workdps(45):
        H = 4 + mp.quad(
            lambda y: y * y * (3 - 2 * y * y + 3 * y**4) * p / (4 + p * (1 - y * y)),
            [0, 1],
        )
        assert mp.re(H) >= mp.mpf(16) / 15
        assert abs(1 / H) <= mp.mpf(15) / 16
    assert "only the active range inverse" in matching.data()["inverse_kernel"]


def test_nonstationary_reference_residual_not_prepared_inverse_input():
    d = inverse.data()
    assert "not a prepared forcing" in d["nonstationary_boundary"]
    assert (
        "nonzero quantum one-point stress" in audit.observable()["background_boundary"]
    )
    assert len(audit.frontier()) == 9 and len(audit.matching_frontier()) == 36
    assert len(audit.residuals()) == 109 and audit.scalar_entry_count() == 140
    assert len(audit.gates()) == 24 and audit.rejected_inputs() == 101
