"""Full-profile, mixed-halfedge and Gaussian/Grassmann diagnostics."""

from fractions import Fraction
from itertools import product

import mpmath as mp
import pytest
import sympy as s
from p8_exceptional_vacuum import analytic
from p8_polynomial_vacuum import model
from p8_vacuum_fermion_local_matching import calibration as old_fermion
from p8_vacuum_full_analytic_target_match import audit as previous
from p8_vacuum_protected_yukawa_profile import (
    audit,
    calibration,
    clock,
    order,
    profile,
    transport,
)


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_named_exact_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_each_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_rejected_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def mp_profile(z, R):
    return z / (1 + (z / R) ** 8) ** (mp.mpf(1) / 8)


def number(value):
    q = s.Rational(value)
    return mp.mpf(int(q.p)) / int(q.q)


@pytest.mark.parametrize(
    "x",
    [
        s.Rational(1, 10**12),
        s.Rational(1, 10),
        s.Rational(1, 2),
        1,
        2,
        10,
        10**6,
        10**12,
    ],
)
@pytest.mark.parametrize("sign", [-1, 1])
def test_full_profile_global_real_caps_and_derivatives(x, sign):
    with mp.workdps(160):
        z = sign * number(x)
        f = mp_profile(z, 1)
        d1 = mp.diff(lambda v: mp_profile(v, 1), z)
        d2 = mp.diff(lambda v: mp_profile(v, 1), z, 2)
        assert abs(f) < 1
        assert 0 < d1 <= 1
        assert abs(d2) <= 9
        assert abs(d1 - (1 + z**8) ** (-mp.mpf(9) / 8)) < mp.mpf("1e-145")
        assert abs(d2 + 9 * z**7 / (1 + z**8) ** (mp.mpf(17) / 8)) < mp.mpf("1e-145")
        assert mp_profile(-z, 1) == -f


def test_full_complex_profile_Cauchy_extraction_through_degree_twenty_five():
    with mp.workdps(150):
        radius = mp.mpf(1) / 8
        count = 96
        angles = [2 * mp.pi * j / count for j in range(count)]
        values = [mp_profile(radius * mp.exp(1j * a), 1) for a in angles]
        for degree in range(26):
            observed = mp.fsum(
                v * mp.exp(-1j * degree * a) for v, a in zip(values, angles)
            ) / (count * radius**degree)
            expected = (
                number(s.binomial(-s.Rational(1, 8), (degree - 1) // 8))
                if degree >= 1 and (degree - 1) % 8 == 0
                else mp.mpf(0)
            )
            assert abs(observed - expected) < mp.mpf("1e-80")


@pytest.mark.parametrize("j", range(8))
def test_eight_finite_complex_branchpoints_not_entire(j):
    z = s.exp(s.I * s.pi * s.Rational(2 * j + 1, 8))
    assert s.simplify(z**8 + 1) == 0
    assert s.simplify(z * s.conjugate(z)) == 1
    assert s.im(z) != 0


@pytest.mark.parametrize("power", [-10, 0, 100, 300, 400, 500])
@pytest.mark.parametrize("sign", [-1, 1])
def test_actual_mass_floor_at_widely_separated_real_arguments(power, sign):
    d = calibration.data()
    with mp.workdps(100):
        R = number(d["fixed_profile_radius"])
        m = number(d["same_fermion_mass"])
        y = mp.sqrt(number(d["actual_Yukawa_squared_rational_upper"]))
        z = sign * mp.mpf(10) ** power
        f = mp_profile(z, R)
        for mass in (m + y * f, m - y * f):
            assert (
                number(
                    d["strict_pointwise_paired_mass_enclosure"][
                        "both_mass_strict_lower"
                    ]
                )
                < mass
            )
            assert mass < number(
                d["strict_pointwise_paired_mass_enclosure"]["both_mass_strict_upper"]
            )


def test_original_linear_crossings_and_profile_protection_are_distinct():
    m, y, R = s.Integer(10), s.Integer(2), s.Integer(1)
    assert m - y * (m / y) == 0
    assert m + y * (-m / y) == 0
    f = profile.expression(m / y, R)
    assert bool(0 < m - y * f < m)
    assert profile.gap(m, y * y, R, s.Rational(1, 4))[
        "both_mass_strict_lower"
    ] == s.Rational(15, 2)


@pytest.mark.parametrize("E", [0, 2, 4, 8, 10])
def test_independent_zero_dimensional_fermion_determinant_loop_grade(E):
    # A Grassmann pair gives D/m^2=1-Y f(phi)^2/m^2.
    # Its first change is +Y phi^10/(4m^2 R^8).
    # In W=h log Z the coefficient at first order in Y has
    # Gaussian mean J/K and variance h/K.
    J, h, K, m, R, Y = s.symbols("J h K m R Y", positive=True)
    moment = sum(
        s.binomial(10, 2 * r)
        * s.factorial2(2 * r - 1)
        * (h / K) ** r
        * (J / K) ** (10 - 2 * r)
        for r in range(6)
    )
    deltaW = s.expand(h * Y * moment / (4 * m * m * R**8))
    coefficient = s.expand(deltaW).coeff(J, E)
    powers = s.Poly(coefficient, h).monoms()
    assert len(powers) == 1
    observed = powers[0][0]
    assert observed == 6 - E // 2
    assert observed == order.new_vertex_floor(E)
    assert coefficient != 0
    if E <= 4:
        assert all(coefficient.coeff(h, loop) == 0 for loop in range(3))


@pytest.mark.parametrize("h", [s.Rational(1, 100), s.Rational(1, 50)])
def test_full_profile_Gaussian_integral_has_predicted_first_difference(h):
    with mp.workdps(70):
        variance = number(h)
        R, m, Y = mp.mpf(3), mp.mpf(10), mp.mpf(1)

        def integrand(x):
            phi = mp.sqrt(variance) * x
            return mp.exp(-x * x / 2) * (phi * phi - mp_profile(phi, R) ** 2)

        observed = (
            Y
            / m**2
            * mp.quad(integrand, [-mp.inf, -4, 0, 4, mp.inf])
            / mp.sqrt(2 * mp.pi)
        )
        leading = 945 * Y * variance**5 / (4 * m * m * R**8)
        assert observed > 0
        assert abs(observed / leading - 1) < mp.mpf("1e-5")


@pytest.mark.parametrize("q2", [0, 2, 7])
def test_paired_full_logdet_higher_field_change_not_an_all_fields_identity(q2):
    with mp.workdps(110):
        m, y, R = mp.mpf(3), mp.mpf(2) / 5, mp.mpf(2)
        phi = mp.mpf("1e-4")
        f = mp_profile(phi, R)
        new = (q2 + (m + y * f) ** 2) * (q2 + (m - y * f) ** 2)
        old = (q2 + (m + y * phi) ** 2) * (q2 + (m - y * phi) ** 2)
        observed = mp.log(new / old) / phi**10
        expected = -y * y * (q2 - m * m) / (2 * R**8 * (q2 + m * m) ** 2)
        assert expected != 0
        assert abs(observed / expected - 1) < mp.mpf("1e-8")


def test_independent_mixed_field_halfedge_enumeration_and_possible_higher_points():
    # Vertices: one new Phi^9 fermion pair, ordinary Phi fermion
    # pair, Phi^4, H Phi^2, and gauge fermion pair. This is a
    # necessary halfedge enumeration; connectivity is not claimed.
    found = {0: [], 1: [], 2: [], 4: [], 8: [], 10: []}
    for Ephi, EH in ((0, 0), (0, 1), (2, 0), (4, 0), (8, 0), (10, 0)):
        for ny, n4, nH, ng in product(range(4), repeat=4):
            phi = 9 + ny + 4 * n4 + 2 * nH - Ephi
            H = nH - EH
            gauge = ng
            fermion = 2 * (1 + ny + ng)
            if min(phi, H, gauge) < 0 or any(n % 2 for n in (phi, H, gauge, fermion)):
                continue
            internal = (phi + H + gauge + fermion) // 2
            vertices = 1 + ny + n4 + nH + ng
            loops = internal - vertices + 1
            if loops < 0:
                continue
            E = Ephi + EH
            found[E].append(loops)
            assert loops >= order.new_vertex_floor(E)
            valences = (11,) + (3,) * ny + (4,) * n4 + (3,) * nH + (3,) * ng
            assert order.grade(E, valences, (0,) * len(valences)) == loops
    for E, values in found.items():
        assert values
        if E <= 4:
            assert min(values) > 2
    assert min(found[8]) == 2
    assert min(found[10]) == 1


@pytest.mark.parametrize(
    "vertices,edges,external,orders",
    [(2, 2, 4, 0), (3, 4, 2, 1), (1, 1, 9, 0), (2, 3, 8, 0), (4, 6, 0, 2)],
)
def test_subgraph_contraction_cannot_reduce_counterterm_loop_weight(
    vertices, edges, external, orders
):
    summed_vertex_weight = s.Rational(2 * edges + external - 2 * vertices, 2) + orders
    loop_grade = edges - vertices + 1 + orders
    local_weight = s.Rational(external - 2, 2) + loop_grade
    assert summed_vertex_weight == local_weight


@pytest.mark.parametrize("k", range(1, 6))
def test_new_EFT_coefficients_have_the_correct_D_dimensional_lift(k):
    eps = s.Symbol("epsilon")
    operator_dim = (8 * k + 1) * (1 - eps) + 2 * (s.Rational(3, 2) - eps)
    coefficient_dim = (4 - 2 * eps) - operator_dim
    assert s.expand(coefficient_dim + 8 * k - (8 * k + 1) * eps) == 0
    assert order.new_vertex_floor(4, k) == 4 * k


@pytest.mark.parametrize("time", [0, s.Rational(1, 7), s.Rational(2, 3), 2])
def test_conditional_mass_derivative_bounds_for_nonlinear_argument(time):
    with mp.workdps(80):
        m, y, R, omega = mp.mpf(10), mp.mpf(1) / 10, mp.mpf(2), mp.mpf(3) / 2
        z = lambda t: R * mp.sin(omega * t)
        mass = lambda t: m + y * mp_profile(z(t), R)
        t = number(time)
        caps = clock.ratios(
            10, s.Rational(1, 10), 2, 3, s.Rational(9, 2), s.Rational(1, 2)
        )
        assert abs(mp.diff(mass, t)) / mass(t) ** 2 < number(
            caps["first_mass_derivative_over_mass_squared_upper"]
        )
        assert abs(mp.diff(mass, t, 2)) / mass(t) ** 3 < number(
            caps["second_mass_derivative_over_mass_cubed_upper"]
        )


def test_literal_flat_clock_map_from_direct_time_differentiation():
    t, kappa, lam, gamma, c = s.symbols("t kappa lambda gamma c", positive=True)
    phi = s.sqrt(kappa) * t
    X = s.diff(phi, t) ** 2
    E = phi + s.diff(phi, t, 2)
    K = lambda expr: s.diff(expr, t, 2) + 2 * expr
    Z = s.diff(phi, t) ** 2 * s.diff(phi, t, 2)
    correction = (
        -c * phi**3 / 6
        + lam * (2 * phi * X + phi**2 * E)
        + gamma * (2 * Z + phi * X + phi**3 / 3)
        - gamma * phi * K(X)
        - gamma * phi * K(phi * E) / 2
    )
    actual = s.expand(phi + correction)
    expected = (1 + 2 * lam * kappa - 2 * gamma * kappa) * phi + (
        lam - 2 * gamma / 3 - c / 6
    ) * phi**3
    assert s.expand(actual - expected) == 0
    derivative = s.diff(actual, t).subs(t, 0)
    assert (
        s.simplify(derivative / s.sqrt(kappa))
        == 1 + 2 * lam * kappa - 2 * gamma * kappa
    )


def test_actual_positive_Yukawa_lower_and_failed_naive_mass_screen():
    d = calibration.data()
    f = old_fermion.data()
    assert f["rational_Yukawa_squared_lower"] > s.Rational(4, 10**206)
    assert bool(f["actual_Yukawa_squared"] > f["rational_Yukawa_squared_lower"])
    p = model.data()["actual_parameters"]
    slope = 1 + 2 * p["lambda"] * analytic.KAPPA - 2 * p["gamma"] * analytic.KAPPA
    assert d["flat_linear_clock_literal_map_slope_at_zero"] == slope
    assert d["flat_linear_clock_literal_map_mass_derivative_ratio_lower"] > 10**97
    assert d["diagnostic_linear_argument_mass_derivative_ratios"][
        "first_mass_derivative_over_mass_squared_upper"
    ] < s.Rational(1, 10**102)
    assert "not a proof of particle production" in d["naive_map_screen"]


def test_exact_rational_profile_gap_and_clock_inputs():
    d = profile.gap(Fraction(10), s.Integer(1), Fraction(1), s.Rational(1, 2))
    assert d["both_mass_strict_lower"] == 5
    assert d["both_mass_strict_upper"] == 15
    assert clock.ratios(10, Fraction(1, 10), 1, 1, 0, s.Rational(1, 2))[
        "first_mass_derivative_over_mass_squared_upper"
    ] == s.Rational(1, 250)


def test_observable_transport_is_named_and_frontier_only_appended():
    assert all(
        row["through_loop"] < order.new_vertex_floor(row["source_endpoints"])
        for row in transport.rows()
    )
    assert audit.frontier() == previous.frontier()
    assert audit.matching()[:-1] == previous.matching()
    assert "NEW_ANALYTIC_EFT_BRANCH" in audit.matching()[-1]["status"]
    assert audit.controls()["original_P8_not_closed"] is True
    assert (
        "not the original globally renormalizable" in transport.data()["not_inferred"]
    )
