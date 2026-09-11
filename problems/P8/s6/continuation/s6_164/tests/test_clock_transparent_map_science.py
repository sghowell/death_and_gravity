"""Independent interpolation, finite-Fourier action and clock diagnostics."""

from fractions import Fraction
from itertools import product

import mpmath as mp
import pytest
import sympy as s
from p8_exceptional_vacuum import analytic
from p8_offshell_vacuum import norms as original_norms
from p8_vacuum_clock_transparent_map import (
    audit,
    calibration,
    counting,
    gate,
    norms,
    transport,
)
from p8_vacuum_protected_yukawa_profile import audit as previous


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_named_exact_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_each_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_rejected_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "x",
    [
        s.Rational(-1, 4),
        0,
        s.Rational(1, 10),
        s.Rational(1, 2),
        s.Rational(9, 10),
        1,
        s.Rational(5, 4),
    ],
)
def test_independent_beta_integral_and_reflection(x):
    with mp.workdps(100):
        z = mp.mpf(int(s.Rational(x).p)) / int(s.Rational(x).q)
        observed = 51480 * mp.quad(lambda t: t**7 * (1 - t) ** 7, [0, z])
        expected = mp.mpf(str(gate.data()["step"].subs(gate.X, x).evalf(105)))
        assert abs(observed - expected) < mp.mpf("1e-90")
        reflected = 51480 * mp.quad(lambda t: t**7 * (1 - t) ** 7, [0, 1 - z])
        assert abs(observed + reflected - 1) < mp.mpf("1e-90")


def test_literal_fifteenth_degree_polynomial_is_not_a_rational_switch():
    T = gate.data()["step"]
    assert s.denom(T) == 1
    assert s.degree(T, gate.X) == 15
    assert set(s.Poly(T, gate.X).monoms()) == {(j,) for j in range(8, 16)}
    assert gate.data()["gate"].subs(gate.X, 2) != 0
    assert gate.data()["gate"].subs(gate.X, 2) > 1


def test_independent_mixed_clock_variation_and_nonzero_eighth_order():
    eps, a, b, r0, r1, r2 = s.symbols("epsilon a b r0 r1 r2")
    q = gate.data()["gate"]
    changed = s.expand(
        q.subs(gate.X, 1 + eps * a + eps**2 * b) * (r0 + eps * r1 + eps**2 * r2)
    )
    for degree in range(8):
        assert s.expand(changed).coeff(eps, degree) == 0
    assert s.factor(s.expand(changed).coeff(eps, 8) - 6435 * a**8 * r0) == 0
    assert changed.subs({a: 1, b: 0, r0: 1, r1: 0, r2: 0}).coeff(eps, 8) == 6435


def test_direct_vacuum_substitution_retains_degree_three_but_changes_nineteen():
    t, P, Y, R, kappa = s.symbols("t P Y R kappa")
    F = t * P + gate.data()["gate"].subs(gate.X, t * t * Y / kappa) * t**3 * R
    difference = s.expand(F - t * P - t**3 * R)
    for degree in range(19):
        assert difference.coeff(t, degree) == 0
    assert difference.coeff(t, 19) == -6435 * R * Y**8 / kappa**8
    assert s.degree(difference, t) == 33
    assert s.expand(F).coeff(t, 3) == R
    assert s.expand(F).coeff(t, 1) == P


def convolve(a, b):
    result = {}
    for i, ci in a.items():
        for j, cj in b.items():
            result[i + j] = result.get(i + j, Fraction(0)) + ci * cj
    return {i: c for i, c in result.items() if c}


def add(a, b):
    result = dict(a)
    for i, c in b.items():
        result[i] = result.get(i, Fraction(0)) + c
    return {i: c for i, c in result.items() if c}


def scale(a, c):
    return {i: c * v for i, v in a.items() if c * v}


def power(a, n):
    value = {0: Fraction(1)}
    for _ in range(n):
        value = convolve(value, a)
    return value


def finite_maps(amplitude=Fraction(1, 3), Rcoefficient=Fraction(1, 100), kappa=100):
    # A periodic one-dimensional diagnostic of the same product and
    # multiplier inequalities, not a claim this cosine is Schwartz.
    psi = {-1: amplitude / 2, 1: amplitude / 2}
    X = {}
    for i, ci in psi.items():
        for j, cj in psi.items():
            X[i + j] = X.get(i + j, Fraction(0)) - i * j * ci * cj / kappa
    cubic = scale(power(psi, 3), Rcoefficient)
    step = {}
    for j, c in gate.data()["step_nonzero_coefficients"].items():
        step = add(step, scale(power(X, j), Fraction(int(c))))
    old = add(psi, cubic)
    change = scale(convolve(step, cubic), -1)
    return psi, old, add(old, change), change


def stationary_action(field, L, g, M):
    quadratic = sum(
        Fraction(k * k - 1, 2) * c * field.get(-k, Fraction(0))
        for k, c in field.items()
    )
    source = convolve(field, field)
    local = -L * convolve(source, source).get(0, Fraction(0)) / 24
    heavy = (
        g
        * sum(c * source.get(-k, Fraction(0)) / (M - k * k) for k, c in source.items())
        / 8
    )
    return quadratic + local + heavy


@pytest.mark.parametrize("amplitude", [Fraction(1, 4), Fraction(1, 3), Fraction(1, 2)])
def test_exact_finite_Fourier_full_stationary_action_comparison(amplitude):
    C, kappa = Fraction(1, 100), 100
    psi, old, new, change = finite_maps(amplitude, C, kappa)
    assert max(map(abs, old)) == 3
    assert max(map(abs, change)) == 33
    assert max(map(abs, convolve(new, new))) == 66
    U2 = sum(c * c for c in psi.values())
    delta = norms.gate_difference_cap(s.Rational(4, kappa)) * s.Rational(C)
    assert sum(abs(c) for c in change.values()) <= delta
    assert sum(c * c for c in change.values()) <= delta * delta * s.Rational(U2)
    L, g, M = Fraction(1, 5), Fraction(1, 7), 5000
    observed = stationary_action(new, L, g, M) - stationary_action(old, L, g, M)
    bound = norms.action_difference(C, delta, L, g, M)[
        "complete_action_difference_coefficient"
    ]
    assert observed != 0
    assert abs(s.Rational(observed)) <= bound * s.Rational(U2)


def test_hidden_new_heavy_pole_detected_if_old_source_radius_is_reused():
    _, old, new, _ = finite_maps()
    assert max(map(abs, convolve(old, old))) == 6
    assert convolve(new, new)[66] != 0
    assert 66**2 > 6**2
    assert stationary_action(old, Fraction(1, 5), Fraction(1, 7), 66**2) is not None
    with pytest.raises(ZeroDivisionError):
        stationary_action(new, Fraction(1, 5), Fraction(1, 7), 66**2)


@pytest.mark.parametrize(
    "C,delta", list(product((0, s.Rational(1, 10), s.Rational(1, 4)), repeat=2))
)
def test_positive_higher_delta_powers_are_bounded_not_dropped(C, delta):
    exact = norms.action_difference(C, delta, s.Rational(2, 3), s.Rational(3, 7), 5000)
    linear = norms.linear_enclosure(C, delta, s.Rational(2, 3), s.Rational(3, 7), 5000)
    for key in (
        "free_action_difference_coefficient",
        "local_quartic_difference_coefficient",
        "full_heavy_resolvent_difference_coefficient",
        "complete_action_difference_coefficient",
    ):
        assert exact[key] <= linear[key]
    if delta:
        assert exact["free_action_difference_coefficient"] > 100 * (1 + C) * delta


def test_actual_finite_coefficient_cap_and_declared_rational_allowances():
    d = calibration.data()
    Xcap = s.Rational(4, analytic.KAPPA)
    actual = norms.gate_difference_cap(Xcap)
    assert actual < d["declared_polynomial_gate_difference_cap"]
    assert (
        actual * original_norms.data()["actual_C_R"]
        < d["declared_map_difference_sup_and_L2_factor"]
    )
    assert d["complete_new_map_action_difference"][
        "complete_action_difference_coefficient"
    ] < s.Rational(1, 10**6790)
    assert d["clock_transparent_full_analytic_target_match_coefficient"] < s.Rational(
        1, 10**800
    )
    assert (
        d["clock_transparent_full_analytic_target_match_coefficient"]
        > d["previous_full_analytic_target_match_coefficient"]
    )


def test_covariant_unit_clock_gate_independent_of_FLRW_scale_factor():
    t, kappa = s.symbols("t kappa", positive=True)
    a = s.Function("a")(t)
    inverse = s.diag(1, -1 / a**2, -1 / a**2, -1 / a**2)
    gradient = s.Matrix([s.sqrt(kappa), 0, 0, 0])
    X = (gradient.T * inverse * gradient)[0] / kappa
    assert s.simplify(X) == 1
    arbitrary_old_R = s.Function("old_R")(t)
    F = s.sqrt(kappa) * t + gate.data()["gate"].subs(gate.X, X) * arbitrary_old_R
    assert F == s.sqrt(kappa) * t
    assert s.diff(F, t) == s.sqrt(kappa)
    assert s.diff(F, t, 2) == 0
    assert counting.data()["clock_variations_identical_through_order"] == 7


def test_clock_source_dictionary_does_not_claim_background_or_state_solution():
    d = transport.data()
    assert d["on_unit_clock_mass_derivative_ratio_bounds"][
        "first_mass_derivative_over_mass_squared_upper"
    ] < s.Rational(1, 10**102)
    assert d["on_unit_clock_mass_derivative_ratio_bounds"][
        "second_mass_derivative_over_mass_cubed_upper"
    ] < s.Rational(1, 10**201)
    assert "does not thereby acquire the target bounce" in d["background_limit"]
    assert "No scalar/gravity parent matching" in d["background_limit"]
    assert "ordinary-Psi" in d["vacuum_prescription"]


@pytest.mark.parametrize("j", range(8, 16))
def test_extra_dimensional_weight_tracks_new_scalar_degree(j):
    eps = s.Symbol("epsilon")
    scalar_degree = 3 + 2 * j
    old_map_epsilon_weight = 2 * eps
    gate_epsilon_weight = 2 * j * eps
    assert (
        s.expand(
            old_map_epsilon_weight + gate_epsilon_weight - (scalar_degree - 1) * eps
        )
        == 0
    )


def test_recoverable_original_scope_and_exact_scalar_helpers():
    assert norms.gate_difference_cap(Fraction(0)) == 0
    assert (
        norms.action_difference(0, 0, 1, 1, 5000)[
            "complete_action_difference_coefficient"
        ]
        == 0
    )
    assert audit.frontier() == previous.frontier()
    assert audit.matching()[:-1] == previous.matching()
    assert audit.controls()["original_P8_not_closed"] is True
