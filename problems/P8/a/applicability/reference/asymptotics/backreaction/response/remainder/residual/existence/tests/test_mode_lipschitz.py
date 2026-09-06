from fractions import Fraction
from math import factorial

import pytest
import sympy as sp
from p8a_existence import mode_lipschitz as modes


def test_mode_response_exact_identity_dictionary():
    checks = modes.identities()
    assert len(checks) == 11
    assert all(value == 0 for value in checks.values()), checks


@pytest.mark.parametrize("value", [True, False, 0.0, 1.0, float("inf"),
                                  float("nan"), sp.Float("0.1"), sp.oo,
                                  -sp.oo, sp.nan, sp.Symbol("x"), 1j])
def test_numeric_interfaces_reject_inexact_boolean_nonfinite_and_symbolic(value):
    with pytest.raises((TypeError, ValueError)):
        modes.rational(value)
    with pytest.raises((TypeError, ValueError)):
        modes.full_history_bounds(value, 1)
    with pytest.raises((TypeError, ValueError)):
        modes.shared_history_bounds(0, 1, value)


@pytest.mark.parametrize("value", [0, 1, -2, "3/7", sp.Rational(2, 9), "0.125"])
def test_exact_rational_inputs_are_retained(value):
    assert modes.rational(value) == sp.Rational(value)


@pytest.mark.parametrize("args", [(2, 1, 1), ("7/2", 1, 1), (3, 0, 1),
                                 (3, -1, 1), (3, 1, 0), (3, 1, -1)])
def test_unintegrable_order_and_empty_frequency_time_domains_rejected(args):
    with pytest.raises(ValueError):
        modes.dyson_coefficients(*args)


def test_rational_majorant_domain_is_explicit_and_closed_at_z_quarter():
    assert modes.full_history_bounds("1/4", 1)["exponential_upper"] == 2
    assert modes.full_history_bounds(0, 3)["derivative"] == 0
    with pytest.raises(ValueError):
        modes.full_history_bounds("251/1000", 1)
    with pytest.raises(ValueError):
        modes.full_history_bounds(-1, 1)
    with pytest.raises(ValueError):
        modes.exponential_upper(1)
    with pytest.raises(ValueError):
        modes.exponential_upper(-1)
    with pytest.raises(ValueError):
        modes.log_ratio_upper("1/2")
    with pytest.raises(ValueError):
        modes.log_ratio_upper(2, 0)
    with pytest.raises(ValueError):
        modes.shared_history_bounds(0, 1, 2)


def test_dyson_radial_integrals_against_separate_fraction_arithmetic():
    duration, split = Fraction(3, 2), Fraction(5, 7)
    derivative_cap = Fraction(1, 100)
    amplitude = derivative_cap * duration
    for order in range(3, 13):
        infrared = split**2 * (2 * duration**2) ** order / (2 * factorial(order))
        ultraviolet = (
            (2 * duration) ** order * split ** (2 - order)
            / (factorial(order) * (order - 2))
        )
        data = modes.derivative_order_bound(
            order, str(derivative_cap), str(duration), str(split)
        )
        assert data["infrared"] == sp.Rational(infrared.numerator, infrared.denominator)
        assert data["ultraviolet"] == sp.Rational(ultraviolet.numerator, ultraviolet.denominator)
        coefficient = order**2 * (infrared + ultraviolet) * amplitude ** (order - 1)
        assert data["lipschitz"] == sp.Rational(coefficient.numerator, coefficient.denominator)


def test_tail_bound_dominates_exact_partial_sums_and_keeps_cubic_origin():
    data = modes.full_history_bounds("1/100", 2)
    partial = sum(
        modes.derivative_order_bound(n, "1/100", 2, "1/2")["lipschitz"]
        for n in range(3, 21)
    )
    assert 0 < partial < data["higher_derivative"]
    assert data["quadratic_derivative"] == sp.Rational(4, 25)
    assert data["higher_derivative"] > 18 * data["strength"] ** 2
    assert modes.full_history_bounds(0, 2)["higher_derivative"] == 0


def test_marked_factor_localization_and_future_value_gain():
    full = modes.full_history_bounds("1/1000", 3)
    for future in [sp.Rational(1, 10), sp.Rational(1, 4), sp.Integer(1), sp.Integer(3)]:
        data = modes.shared_history_bounds("1/1000", 3, future)
        partial = sum(
            modes.shared_order_bound(n, "1/1000", 3, future, "1/3")["lipschitz"]
            for n in range(3, 17)
        )
        assert 0 < partial < data["higher_derivative"]
        assert data["higher_derivative"] == future * full["higher_derivative"] / 3
        assert data["value"] == future * data["derivative"]
        assert data["derivative"] <= full["derivative"]


def test_marked_simplex_and_two_chain_binomial_factors():
    location, duration = sp.symbols("s T", positive=True)
    for order in range(1, 10):
        marked = sum(
            location**j * (duration - location) ** (order - 1 - j)
            / (sp.factorial(j) * sp.factorial(order - 1 - j))
            for j in range(order)
        )
        assert sp.factor(marked - duration ** (order - 1) / sp.factorial(order - 1)) == 0
        product_sum = sum(
            sp.S.One / (sp.factorial(j) * sp.factorial(order - j))
            for j in range(order + 1)
        )
        insertion_sum = 2 * sum(
            sp.S.One / (sp.factorial(j) * sp.factorial(order - 1 - j))
            for j in range(order)
        )
        assert product_sum == 2**order / sp.factorial(order)
        assert insertion_sum == order * product_sum


def test_abel_derivative_error_has_the_required_uniform_integrable_limit():
    radius, cutoff, duration = sp.symbols("r c T", positive=True)
    error_mass = sp.integrate(cutoff**2 / (radius**2 + cutoff**2), (radius, 0, duration))
    assert sp.simplify(error_mass - cutoff * sp.atan(duration / cutoff)) == 0
    assert sp.limit(error_mass, cutoff, 0, dir="+") == 0
    singular_mass = sp.integrate(sp.log(1 + cutoff**2 / radius**2), (radius, 0, sp.oo))
    assert singular_mass == sp.pi * cutoff


def test_log_majorant_uses_an_explicit_positive_tail():
    assert modes.log_ratio_upper(1) == 0
    for ratio in [sp.Rational(3, 2), sp.Integer(2), sp.Integer(12), sp.Integer(100)]:
        argument = (ratio - 1) / (ratio + 1)
        lower = sum(2 * argument ** (2 * j + 1) / (2 * j + 1) for j in range(8))
        assert lower < modes.log_ratio_upper(ratio, 16) < modes.log_ratio_upper(ratio, 8)
        assert modes.log_ratio_upper(ratio, 8) < ratio - 1


def test_original_prepared_derivative_cap_has_margin_and_rational_enclosures():
    data = modes.calibration()
    expected_m = Fraction(2, 10**14) * Fraction(61013499, 8192)
    strength = expected_m * 27
    exponential = 1 / (1 - 2 * strength)
    full = 2 * strength + 18 * strength**2 * exponential
    assert data["full"]["derivative"] == sp.Rational(full.numerator, full.denominator)
    assert data["full"]["derivative_cap"] > data["delta"] * data["prepared_U1"]
    assert data["full_rounding_gap"] > 0
    assert data["shared_rounding_gap"] > 0
    assert data["full"]["strength"] < modes.Z_MAX
    assert data["shared"]["derivative"] < data["full"]["derivative"]


def test_first_dyson_derivative_transfer_including_nonzero_start_control():
    time, inner, frequency = sp.symbols("t s k", positive=True)
    start_value, slope = sp.symbols("u0 u1", real=True)
    kernel = sp.sin(frequency * (time - inner)) / frequency
    potential = start_value + slope * inner
    first = -sp.integrate(kernel * potential * sp.exp(-sp.I * frequency * inner),
                          (inner, 0, time))
    inserted = -sp.integrate(kernel * slope * sp.exp(-sp.I * frequency * inner),
                             (inner, 0, time))
    defect = sp.simplify(sp.expand_complex(
        sp.diff(first, time) + sp.I * frequency * first - inserted
    ))
    endpoint = -start_value * sp.sin(frequency * time) / frequency
    assert sp.simplify(defect - endpoint) == 0
    assert sp.simplify(defect.subs(start_value, 0)) == 0
    assert sp.simplify(defect.subs({start_value: 1, frequency: 1})) != 0


def test_quadratic_kernel_for_a_genuinely_zero_initial_linear_potential():
    time, inner = sp.symbols("t r", positive=True)
    slope = sp.symbols("m", real=True)
    # Integrate the derivative kernel after exchanging the two time integrals.
    derivative = sp.integrate(
        slope**2 * inner * (time**2 - inner**2) / (2 * (time - inner)),
        (inner, 0, time),
    )
    assert sp.simplify(derivative - 5 * slope**2 * time**3 / 12) == 0
    response = 5 * slope**2 * time**4 / 48
    assert sp.simplify(sp.diff(response, time) - derivative) == 0
    # A wrong factor two in F=2k|v|^2 cannot pass this comparison.
    assert sp.simplify(sp.diff(2 * response, time) - derivative) != 0


def test_scaling_keeps_the_lipschitz_constant_dimensionless():
    reference = modes.full_history_bounds("1/1000", 3)
    for rescaling in [sp.Rational(1, 3), sp.Integer(2), sp.Integer(7)]:
        changed = modes.full_history_bounds(
            sp.Rational(1, 1000) / rescaling**3, 3 * rescaling
        )
        assert changed["derivative"] == reference["derivative"]
        assert changed["value"] == rescaling * reference["value"]
        old_shared = modes.shared_history_bounds("1/1000", 3, "1/4")
        new_shared = modes.shared_history_bounds(
            sp.Rational(1, 1000) / rescaling**3,
            3 * rescaling,
            rescaling / 4,
        )
        assert new_shared["derivative"] == old_shared["derivative"]


def test_no_unproved_scheme_or_full_stress_parameters_are_accepted():
    with pytest.raises(TypeError):
        modes.full_history_bounds(0, 1, gamma=1)
    with pytest.raises(TypeError):
        modes.shared_history_bounds(0, 1, 1, initial_state="arbitrary")
