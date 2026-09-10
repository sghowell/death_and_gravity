"""Independent Fraction anchors, derivative bounds and reference-map tests."""

from fractions import Fraction
from math import factorial

import pytest
import sympy as sp
from p8_vacuum_fermion_local_matching import (
    anchors,
    audit,
    calibration,
    potential,
    reference,
    twopoint,
)


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_named_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_every_explicit_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[x[0] for x in audit.bad_cases()]
)
def test_invalid_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("n", range(17))
def test_moments_independently_with_integer_factorials(n):
    assert anchors.moment(n) == sp.Rational(
        Fraction(factorial(n) ** 2, factorial(2 * n + 1))
    )


@pytest.mark.parametrize("mass", (36, 10**6, 10**200))
@pytest.mark.parametrize("terms", range(1, 9))
def test_anchor_intervals_independent_fraction_and_nesting(mass, terms):
    a = {
        n: Fraction(3 * factorial(n) ** 2, n * (2 * n + 3) * factorial(2 * n + 1))
        for n in range(1, terms + 1)
    }
    S = sum(a[n] / mass ** (2 * n) for n in a)
    P = sum((n + 1) * a[n] / mass ** (2 * n) for n in a)
    R = sum(n * a[n] / mass ** (2 * n) for n in a)
    ratio = Fraction(1, 4 * mass * mass)
    tail = Fraction(3, 5) * ratio ** (terms + 1) / (1 - ratio)
    expected = {
        "mass_increment_over_2NY_div_Q": (
            Fraction(2, 3) - S - tail,
            Fraction(2, 3) - S,
        ),
        "slope_over_2NY_div_Q": (Fraction(2, 3) - P - 2 * tail, Fraction(2, 3) - P),
        "OS_zero_remainder_over_2NY_div_Q": (-R - tail, -R),
    }
    d = anchors.enclosure(mass, terms)
    for prefix, (low, high) in expected.items():
        assert d[prefix + "_lower"] == sp.Rational(low)
        assert d[prefix + "_upper"] == sp.Rational(high)
        assert low < high
        if terms > 1:
            previous = anchors.enclosure(mass, terms - 1)
            assert previous[prefix + "_lower"] < d[prefix + "_lower"]
            assert d[prefix + "_upper"] < previous[prefix + "_upper"]
    assert d["common_positive_series_tail_upper"] == sp.Rational(tail)


@pytest.mark.parametrize("degree", range(6, 35, 2))
def test_positive_higher_field_coefficients_and_monotonicity(degree):
    d = potential.data()
    formula = d["positive_higher_even_field_coefficient"]
    symbol = next(iter(formula.free_symbols))
    actual = formula.subs(symbol, degree)
    independent = Fraction(288 * factorial(degree - 5), factorial(degree))
    assert actual == sp.Rational(independent) > 0
    assert formula.subs(symbol, degree + 2) < actual
    assert actual <= sp.Rational(2, 5)


@pytest.mark.parametrize("A", (sp.Rational(0), sp.Rational(1, 8), sp.Rational(1, 4)))
@pytest.mark.parametrize("s", (-1, 3, 1 + 2 * sp.I))
def test_complex_parameter_denominator_and_absolute_numerator(A, s):
    m = sp.Integer(36)
    numerator = m * m * (4 * A * A - 2 * A) + s * A * A
    absolute_upper = m * m * (4 * A * A + 2 * A) + 3 * A * A
    assert sp.simplify(numerator * sp.conjugate(numerator)) <= absolute_upper**2
    denominator = m * m - A * s
    assert sp.re(denominator) >= m * m - sp.Rational(3, 4) > 0


@pytest.mark.parametrize("mass", (36, 100, 10**200))
def test_zero_remainder_sign_and_uniform_bound_independently(mass):
    d = anchors.enclosure(mass)
    B_over_C = (Fraction(7, 15) * mass * mass + Fraction(1, 10)) / (
        2 * (mass * mass - Fraction(3, 4)) ** 2
    )
    low = d["OS_zero_remainder_over_2NY_div_Q_lower"]
    high = d["OS_zero_remainder_over_2NY_div_Q_upper"]
    assert -sp.Rational(B_over_C) < low < high < 0
    assert (
        0
        < d["slope_over_2NY_div_Q_lower"]
        < d["slope_over_2NY_div_Q_upper"]
        < sp.Rational(2, 3)
    )


def test_potential_low_coefficients_from_independent_log_series():
    r = sp.symbols("r")
    logs = sum((-1) ** (n + 1) * r**n / sp.Integer(n) for n in range(1, 9))
    single = (1 + r) ** 4 * (2 * logs - sp.Rational(3, 2))
    result = sp.expand(-3 * (single + single.subs(r, -r)))
    assert [result.coeff(r, n) for n in (0, 2, 4, 6, 8)] == [
        9,
        12,
        -16,
        sp.Rational(2, 5),
        sp.Rational(3, 70),
    ]


def test_inert_vacuum_energy_not_lost():
    d = potential.data()
    phi = d["reference_field"]
    m = d["fermion_mass"]
    scale = next(
        s
        for s in d["inert_one_loop_MSbar_potential"].free_symbols
        if str(s) == "positive_reference_scale"
    )
    active = d["active_one_loop_MSbar_potential"].subs({phi: 0, scale: m})
    inert = d["inert_one_loop_MSbar_potential"].subs(scale, m)
    assert sp.simplify(inert - 6 * active) == 0
    assert sp.simplify(active + inert - d["all_flavor_vacuum_energy_at_scale_m"]) == 0


def test_actual_mass_reference_is_not_the_physical_pole():
    d = calibration.data()
    assert d["mass_reference_parameter_bounds"][1] < 0
    assert 0 < d["canonical_kinetic_normalization_bounds"][0]
    assert 0 < d["canonical_OS_remainder_coefficient_upper"] < sp.Rational(1, 10**606)
    assert d["canonical_potential_mass_squared_lower"] > 0
    assert d["unnormalized_selected_functional_quartic_lower"] > 0
    assert "tree-scalar plus one-loop-fermion functional" in d["scope"]
    assert "not the complete quantum candidate" in d["scope"]


def test_reference_field_and_potential_couplings_not_old_numbers():
    d = reference.data()
    fp = next(s for s in d["kinetic_normalization"].free_symbols)
    G = next(
        s
        for s in d["canonical_heavy_cubic"].free_symbols
        if str(s) == "cubic_reference"
    )
    assert sp.simplify(
        (d["canonical_heavy_cubic"] / G).subs(fp, sp.Rational(1, 4))
    ) == sp.Rational(4, 3)
    assert "not an all-loop" in d["scope"]
    assert d["canonical_heavy_cubic"].is_number is False
    assert "Old scalar loops" in potential.data()["scope"]
    assert (
        "entire" in twopoint.data()["scope"].lower()
        or "full tadpole" in twopoint.data()["scope"]
    )


def test_exact_counts_and_nonclosure():
    assert len(audit.residuals()) == 67
    assert len(audit.gates()) == 32
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 30
    assert audit.controls()["original_P8_not_closed"] is True
