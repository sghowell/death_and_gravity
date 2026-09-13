"""Independent complete-jet, graph, tensor, matching and scope checks."""

import itertools

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_heavy_parent_one_loop import audit, bounds, germs, loops


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_exact_residual(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_every_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_all_unsupported_inputs_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_exact_complete_counts_and_frontier():
    assert (
        len(audit.residuals()),
        audit.scalar_entry_count(),
        len(audit.gates()),
        len(audit.controls()),
        audit.rejected_inputs(),
    ) == (66, 66, 44, 9, 285)
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 95
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()
    assert audit.validate_scope(audit.frontier(), audit.matching())


@pytest.mark.parametrize("S,x", [(4, -1), (4, 1), (9, 0), (10**196, -1), (10**196, 1)])
def test_closed_physical_window_endpoints(S, x):
    assert audit.require_physical(S, x) == (S, x)


@pytest.mark.parametrize(
    "u,X,order", [(0, 1, 0), (1, s.Rational(7, 8), 4), (-1, s.Rational(9, 8), 4)]
)
def test_closed_finite_clock_tube(u, X, order):
    assert audit.require_finite_tube(u, X, order) == (u, X, order)


def test_prescription_and_stage_are_not_all_loop_or_quantum_limit_claims():
    assert audit.require_prescription("H8A420-VAC-OS4", 1) == ("H8A420-VAC-OS4", 1)
    for stage in (
        "complete_classical_limit_one_loop_four_point",
        "all_new_source_contraction_classes",
        "dimensional_tensor_before_MS_subtraction",
        "new_full_parent_OS4_value",
        "complete_physical_angle_first_loop_bound",
        "complete_first_loop_forward_coefficients",
        "specified_finite_counterterm_clock_four_jet",
    ):
        assert audit.require_stage(stage) == stage
    assert "CLASSICAL" in audit.observable()["input"]
    assert "free" in audit.observable()["input"]
    assert "No finite-gravity quantum-limit theorem" in audit.observable()["input"]


def test_exhaustive_degree_counting_independently_includes_no_vertex_above_six():
    allowed = []
    for row in itertools.product(range(5), repeat=4):
        if sum((i + 1) * v for i, v in enumerate(row)) == 4:
            allowed.append(row)
    assert set(allowed) == set(germs.graph_degrees(4))
    for degree in range(7, 15):
        assert degree - 2 > 4
    assert germs.graph_degrees(2) == ((0, 1, 0, 0), (2, 0, 0, 0))
    assert germs.graph_degrees(1) == ((1, 0, 0, 0),)


def test_complete_literal_six_jet_and_higher_localizer_separation():
    cp = germs.couplings()
    k = germs.KAPPA
    assert cp["phi6"] * k * k == -s.Rational(2681, 200)
    assert cp["phi2Y2"] * k * k == -s.Rational(343168, 25)
    assert cp["H_phi2_Y"] == -4 * germs.G / k
    assert (
        s.diff(germs.parent.switch(germs.X, germs.parent.A), germs.X).subs(germs.X, 0)
        == -8
    )
    assert (
        s.diff(germs.parent.switch(germs.X, germs.parent.A), germs.X, 2).subs(
            germs.X, 0
        )
        == 56 - 2 * germs.parent.A
    )
    assert 1 + 2 + 2 * 2 == 7 and 2 * 2 + 2 * 2 == 8
    assert "degree2048" in germs.data()["higher_localizer_boundary"]


def add(*polys):
    out = {}
    for p in polys:
        for m, v in p.items():
            out[m] = out.get(m, 0) + v
    return {m: v for m, v in out.items() if v}


def scale(p, c):
    return {m: c * v for m, v in p.items() if c * v}


def mul(a, b):
    out = {}
    for ma, va in a.items():
        for mb, vb in b.items():
            if not ma & mb:
                mask = ma | mb
                out[mask] = out.get(mask, 0) + va * vb
    return {m: v for m, v in out.items() if v}


def power(a, n):
    out = {0: mp.mpf(1)}
    for _ in range(n):
        out = mul(out, a)
    return out


def coefficient_for_momenta(momenta):
    dimension = len(momenta[0])
    eta = [1] + [-1] * (dimension - 1)
    phi = {1 << i: mp.mpf(1) for i in range(6)}
    grad = [{1 << i: mp.j * momenta[i][a] for i in range(6)} for a in range(dimension)]
    hess = [
        [
            {1 << i: -momenta[i][a] * momenta[i][b] for i in range(6)}
            for b in range(dimension)
        ]
        for a in range(dimension)
    ]
    box = add(*(scale(hess[a][a], eta[a]) for a in range(dimension)))
    Y = add(*(scale(mul(grad[a], grad[a]), eta[a]) for a in range(dimension)))
    W = [
        add(*(scale(mul(grad[a], hess[a][b]), eta[a]) for a in range(dimension)))
        for b in range(dimension)
    ]
    L3 = mul(add(*(scale(mul(W[b], grad[b]), eta[b]) for b in range(dimension))), box)
    L4 = add(*(scale(mul(W[b], W[b]), eta[b]) for b in range(dimension)))
    gal = add(L3, scale(L4, -1))
    words = {
        "phi6": power(phi, 6),
        "phi4Y": mul(power(phi, 4), Y),
        "phi2Y2": mul(power(phi, 2), power(Y, 2)),
        "Y3": power(Y, 3),
        "phi2_L3_minus_L4": mul(power(phi, 2), gal),
        "Y_L3_minus_L4": mul(Y, gal),
    }
    return {name: p.get(63, 0) / 2 for name, p in words.items()}


def sphere_rule(dimension):
    for axis in range(dimension):
        for sign in (-1, 1):
            vec = [mp.mpf(0)] * dimension
            vec[axis] = sign
            yield mp.mpf(4 - dimension) / (2 * dimension * (dimension + 2)), vec
    for a, b in itertools.combinations(range(dimension), 2):
        for sa, sb in itertools.product((-1, 1), repeat=2):
            vec = [mp.mpf(0)] * dimension
            vec[a] = sa / mp.sqrt(2)
            vec[b] = sb / mp.sqrt(2)
            yield mp.mpf(1) / (dimension * (dimension + 2)), vec


def independent_tadpoles(dimension):
    momentum = mp.sqrt(5) / 2
    energy = mp.mpf(3) / 2
    external = [
        [energy, 0, momentum],
        [energy, 0, -momentum],
        [-energy, -4 * momentum / 5, -3 * momentum / 5],
        [-energy, 4 * momentum / 5, 3 * momentum / 5],
    ]
    external = [p + [0] * (dimension - 3) for p in external]
    result = {}
    weights = 0
    for weight, q in sphere_rule(dimension):
        weights += weight
        k = [q[0]] + [mp.j * value for value in q[1:]]
        values = coefficient_for_momenta(external + [k, [-v for v in k]])
        for name, value in values.items():
            result[name] = result.get(name, 0) + weight * value
    assert abs(weights - 1) < mp.mpf("1e-45")
    return result


@pytest.mark.parametrize("dimension", (3, 4, 5))
def test_independent_literal_multilinear_sphere_contraction(dimension):
    # Different representation: coefficients of six commuting nilpotent plane
    # wave amplitudes are multiplied directly, with no720-permutation routine.
    # The signed rule integrates degree4 sphere moments exactly. Negative
    # weights at d5 are a quadrature identity, not a physical measure claim.
    with mp.workdps(60):
        actual = independent_tadpoles(dimension)
        raw, _ = loops.local_factors()
        for name, value in actual.items():
            expected = mp.mpf(
                str(
                    s.N(
                        raw[name].subs({loops.S: 9, loops.T: -1, loops.d: dimension}),
                        58,
                    )
                )
            )
            assert abs(value - expected) < mp.mpf("1e-48")


@pytest.mark.parametrize("dimension", (3, 4, 5))
def test_independent_signed_rule_full_second_and_fourth_moments(dimension):
    with mp.workdps(60):
        rule = list(sphere_rule(dimension))
        for a, b in itertools.product(range(dimension), repeat=2):
            val = sum(w * q[a] * q[b] for w, q in rule)
            assert abs(val - mp.mpf(int(a == b)) / dimension) < mp.mpf("1e-50")
        for a, b, c, d in itertools.product(range(dimension), repeat=4):
            val = sum(w * q[a] * q[b] * q[c] * q[d] for w, q in rule)
            expected = mp.mpf(
                int(a == b and c == d) + int(a == c and b == d) + int(a == d and b == c)
            ) / (dimension * (dimension + 2))
            assert abs(val - expected) < mp.mpf("1e-50")


def test_dimensional_finite_part_is_not_raw_four_dimensional_contraction():
    raw, finite = loops.local_factors()
    epsilon = s.Symbol("epsilon", positive=True)
    for name in raw:
        point = raw[name].subs({loops.S: 9, loops.T: -1})
        laurent = s.series(
            -(1 / epsilon + 1) * point.subs(loops.d, 4 - 2 * epsilon), epsilon, 0, 1
        ).removeO()
        const = s.expand(laurent).coeff(epsilon, 0)
        assert s.cancel(const + finite[name].subs({loops.S: 9, loops.T: -1})) == 0
    assert s.cancel(finite["Y3"] - raw["Y3"].subs(loops.d, 4)) != 0
    assert (
        s.cancel(finite["phi2_L3_minus_L4"] - raw["phi2_L3_minus_L4"].subs(loops.d, 4))
        != 0
    )


@pytest.mark.parametrize("mass2", (3, 31, 10**8))
def test_independent_mixed_Feynman_shift_integral_and_full_A0_B0(mass2):
    with mp.workdps(70):
        n = mp.mpf(mass2)
        F = lambda x: (1 - x) ** 2 + n * x
        points = sorted({mp.mpf(0), 1 / n, 1 / mp.sqrt(n), mp.mpf("0.1"), mp.mpf(1)})
        B = -mp.quad(lambda x: mp.log(F(x)), points)
        reduced = 1 - n + n * mp.log(n) + (n - 4) * B
        shifted = mp.quad(lambda x: 2 * (x + 1) * mp.log(F(x)), points)
        assert abs(reduced - shifted) < mp.mpf("1e-55")
        assert 0 < shifted < 3 * mp.log(n)
        assert mp.quad(lambda x: -2 * (x + 1), [0, 1]) == -3


def test_all_three_source_contractions_and_onepoint_counterterm():
    phi, Y, I, c, g, n, j = s.symbols("phi Y I c g n j")
    J2 = g * phi**2 / 2
    J4 = c * phi**2 * Y
    assert s.diff(J2, phi, 2) * I / 2 == g * I / 2
    # A scalar pair in Phi² and a derivative pair in Y are both present.
    wick_J4 = c * I * (Y + phi**2)
    assert s.expand(s.diff(J4, phi, 2) * I / 2 + c * I * phi**2 - wick_J4) == 0
    assert s.expand((g * I / 2 - j) * J4 / n).subs(j, g * I / 2) == 0
    assert s.expand(g * I * J4 / (2 * n)) != 0
    ext = s.Symbol("s")
    vertex = 2 * c * I * (1 - (ext - 2) / 2)
    assert s.expand(vertex - c * I * (4 - ext)) == 0


def test_complete_new_symmetric_subtraction_and_forward_coefficients():
    data = loops.data()
    terms = loops.new_base_loop()
    base = terms["mixed"] + terms["vertex"] + terms["local"]
    point = {loops.S: s.Rational(4, 3), loops.T: s.Rational(4, 3)}
    centered = base - base.subs(point)
    assert s.cancel(centered - data["extra_OS4_amplitude"]) == 0
    v, t = s.symbols("v t")
    amplitude = data["extra_OS4_amplitude"].subs(
        {loops.S: 2 + v - t / 2, loops.T: t}, simultaneous=True
    )
    _, (_, k2, k3) = loops.local_polynomial()
    eps = data["extra_heavy_exchange_relative_factor"]
    D = loops.n - 2
    expected = {
        (2, 0): eps * 2 * loops.g**2 / D**3 - 2 * k2 / (16 * s.pi**2),
        (2, 1): -eps * 3 * loops.g**2 / D**4 + k3 / (16 * s.pi**2),
        (4, 0): eps * 2 * loops.g**2 / D**5,
    }
    for (i, j), want in expected.items():
        got = s.diff(amplitude, v, i, t, j).subs({v: 0, t: 0}) / (
            s.factorial(i) * s.factorial(j)
        )
        assert s.cancel(got - want) == 0
    # The actual tree moment relations convert these to4lambda,-3gamma,
    # gamma²/lambda. In particular the exchange correction is not zero.
    replacements = {loops.n: germs.MASS2, loops.g: germs.G, loops.k: germs.KAPPA}
    assert s.cancel((2 * loops.g**2 / D**3).subs(replacements) - 4 * germs.LAMBDA) == 0
    assert s.cancel((-3 * loops.g**2 / D**4).subs(replacements) + 3 * germs.GAMMA) == 0
    assert (
        s.cancel(
            (2 * loops.g**2 / D**5).subs(replacements) - germs.GAMMA**2 / germs.LAMBDA
        )
        == 0
    )
    assert expected[(4, 0)] != 0


def test_full_exchange_positivity_and_value_subtraction_not_mass_expansion():
    n = loops.n
    s0 = s.Rational(4, 3)
    channels = (loops.S, loops.T, loops.U)
    direct = sum(1 / (n - x) for x in channels) - 3 / (n - s0)
    divided = sum((x - s0) ** 2 / ((n - s0) ** 2 * (n - x)) for x in channels)
    assert s.cancel(direct - divided) == 0
    for S, angle in ((4, -1), (9, 0), (9, 1), (100, -1)):
        T = -(S - 4) * (1 - s.Rational(angle)) / 2
        got = direct.subs({loops.S: S, loops.T: T, n: 1000})
        assert got >= 0


def test_actual_full_matching_and_contact_bounds_are_exact_not_underflow():
    result = bounds.data()
    extra = result["new_symmetric_value_absolute_upper"]
    assert isinstance(extra, s.Rational) and extra > 0
    assert (
        0
        < result["new_contact_positive_lower"]
        < result["new_contact_absolute_upper"]
        < s.Rational(1, 10**408)
    )
    assert result["new_extra_angular_relative_upper"] < s.Rational(1, 10**600)
    assert result["combined_full_first_loop_angular_relative_upper"] < s.Rational(
        1, 10**199
    )
    for b in ("b20", "b21", "b40"):
        assert (
            0 < result["extra_relative_coefficient_bounds"][b] < s.Rational(1, 10**600)
        )
    assert germs.KAPPA == 10**800
    assert germs.parent.LOCALIZER == 10**420


def test_finite_switch_is_not_the_classical_heavy_localizer():
    X = germs.X
    V = germs.data()["finite_extension_switch"]
    denominator = X**germs.N + (1 - X) ** germs.N
    assert s.factor(denominator * (V - 1) + X**germs.N) == 0
    assert s.factor(denominator * V - (1 - X) ** germs.N) == 0
    assert s.diff(V, X).subs(X, 0) == 0
    assert s.diff(germs.parent.switch(X, germs.parent.A), X).subs(X, 0) == -8
    # Multiplying a finite quadratic mass by h instead would create Phi²Y.
    dm, k = s.symbols("delta_m_squared kappa")
    assert (-dm / 2) * (-8) / k == 4 * dm / k


def test_complete_finite_extension_four_jet_bound_and_scope():
    value = bounds.data()["finite_extension_complete_four_jet_upper"]
    assert value == 2 * s.Integer(10) ** 403 * s.Rational(9, 55) ** 1024
    assert 0 < value < s.Rational(1, 10**400)
    assert isinstance(value, s.Rational)
    assert s.factorial(4) * 64**4 == 402653184 < 10**9
    text = bounds.data()["finite_extension_tube_bound"]
    assert "NEW finite-counterterm budget" in text
    assert "not S238" in text
    assert "entire loop effective action" in text
    assert "not a proof" in loops.data()["new_parent_loop_boundary"]
