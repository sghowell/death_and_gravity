"""Independent Fraction, kernel, boundary and insufficiency regressions."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_vacuum_gauge_light_cut_subtraction import (
    analytic,
    audit,
    calibration,
    domain,
    logarithm,
    subtraction,
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


@pytest.mark.parametrize("n", range(1, 9))
def test_log_interval_independent_fraction_and_nesting(n):
    r = Fraction(1, 3)
    low = 2 * sum(r ** (2 * j + 1) / Fraction(2 * j + 1) for j in range(n))
    tail = 2 * r ** (2 * n + 1) / (Fraction(2 * n + 1) * (1 - r * r))
    d = logarithm.enclosure(n)
    assert d["lower"] == sp.Rational(low)
    assert d["upper"] == sp.Rational(low + tail)
    if n > 1:
        previous = logarithm.enclosure(n - 1)
        assert previous["lower"] < d["lower"] < d["upper"] < previous["upper"]


@pytest.mark.parametrize("mass", (36, 100, 10**6, 10**100, 10**200))
def test_finite_mass_coefficient_interval_independent_fraction(mass):
    r = Fraction(1, 3)
    low_log = 2 * sum(r ** (2 * j + 1) / Fraction(2 * j + 1) for j in range(8))
    tail = 2 * r**17 / (17 * (1 - r * r))
    eps = Fraction(4000000, mass)
    eta = Fraction(5, 2) * eps + eps * eps / 8
    error = 288 * eta
    low = 2 * low_log - Fraction(21, 4) - error
    high = 2 * (low_log + tail) - Fraction(21, 4) + error
    d = calibration.point(mass)
    assert d["cut_b2_over_K_lower"] == sp.Rational(low)
    assert d["cut_b2_over_K_upper"] == sp.Rational(high)
    assert d["strictly_negative_subtraction_coefficient"] == (high < 0)


@pytest.mark.parametrize("n", range(11))
def test_removable_quotient_independent_polynomial_division(n):
    d = analytic.data()["independent_polynomial_divided_difference_checks"][n]
    x, s = sp.symbols("cut_variable channel_s")
    quotient, remainder = sp.div(x**n - s**n, x - s, x)
    assert remainder == 0
    assert sp.expand(quotient - d["removable_quotient"]) == 0
    assert sp.integrate(quotient, (x, 0, 6)) == d["integral"]


def test_leading_center_coefficient_from_local_power_series():
    v = sp.symbols("center_displacement")
    s, u = 2 + v, 2 - v
    single = lambda w: 18 + 6 * w + w * w * (sp.log(6 - w) - sp.log(w))
    germ = single(s) + single(u) + sp.I * sp.pi * (s * s - u * u)
    coefficient = sp.expand(sp.series(germ, v, 0, 3).removeO()).coeff(v, 2)
    assert sp.simplify(coefficient - 2 * sp.log(2) + sp.Rational(21, 4)) == 0


def test_wrong_same_sign_crossed_germ_changes_center_coefficient():
    s = sp.symbols("channel_s")
    u = 4 - s
    correct_imag = sp.I * sp.pi * (s * s - u * u)
    wrong_imag = sp.I * sp.pi * (s * s + u * u)
    assert sp.diff(correct_imag, s, 2) == 0
    assert sp.diff(wrong_imag, s, 2) / 2 == 2 * sp.I * sp.pi


def test_known_cut_does_not_fix_regular_matching_coefficient():
    s, c = sp.symbols("channel_s arbitrary_regular_coefficient")
    regular = c * (s - 2) ** 2
    assert sp.expand(regular.subs(s, 4 - s) - regular) == 0
    assert sp.diff(regular, s, 2) / 2 == c
    # No retuning of the actual parent action is performed by this diagnostic.
    assert "not permission" in subtraction.data()["scope"]


def test_large_complex_disc_is_not_only_real_interval_control():
    d = domain.data()
    a = analytic.data()
    assert d["complex_s_disc_radius"] == a["rho_control_radius"] == 5
    assert a["quotient_control_radius"] == 4
    assert a["coefficient_Cauchy_radius"] == 1
    assert d["strict_cumulative_component_norm_upper"] == 18
    assert d["minimum_fermion_mass"] == 36


def test_actual_strict_cut_band_and_uncomputed_full_amplitude():
    d = calibration.data()
    point = d["actual"]
    assert -4 < point["cut_b2_over_K_lower"] < point["cut_b2_over_K_upper"] < -3
    assert point["finite_mass_cut_b2_error_over_K_upper"] < sp.Rational(1, 10**190)
    assert d["absolute_subtraction_b2_rational_upper"] < sp.Rational(4, 10**1620)
    assert "not a new-model amplitude budget" in d["scope"]


def test_valid_small_mass_bound_can_be_inconclusive():
    point = calibration.point(36)
    assert point["strictly_negative_subtraction_coefficient"] is False
    assert point["strict_minus_four_to_minus_three_K_band"] is False


def test_exact_counts_and_nonclosure():
    assert len(audit.residuals()) == 80
    assert len(audit.gates()) == 29
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 24
    assert audit.controls()["original_P8_not_closed"] is True
