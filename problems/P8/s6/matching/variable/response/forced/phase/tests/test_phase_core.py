"""Exact physical moment, continuous error and strict-domain tests."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_forced_phase import core, independent, verify


def test_exact_physical_weight_and_moment_identities():
    checks = core.checks()
    assert len(checks["residuals"]) == 8
    assert all(value == 0 for value in checks["residuals"].values())
    assert all(value > 0 for value in checks["margins"].values())


@pytest.mark.parametrize("denominator", [1000, 10000, 10**6])
@pytest.mark.parametrize("deficit", [Fraction(0), Fraction(1, 200), Fraction(1, 100)])
def test_separate_Fraction_error_chain(denominator, deficit):
    a = core.calibration(sp.Rational(1, denominator), sp.Rational(deficit))
    b = independent.calibration(Fraction(1, denominator), deficit)
    assert a.keys() == b.keys()
    assert all(a[key] == sp.Rational(b[key]) for key in a)


def test_all_radius_proof_uses_binomial_coefficients_and_monotone_factors():
    checked = independent.checks()
    coeffs = checked["inverse_ninth_nonnegative_coefficients"]
    t = sp.Symbol("t", nonnegative=True)
    actual = sp.Poly(1-(1-9*t)*(1+t)**9, t)
    assert [sp.Rational(value) for value in coeffs] == [actual.nth(i) for i in range(11)]
    assert all(value >= 0 for value in coeffs)
    assert checked["rectangular_moment_lower"] == Fraction(5, 28)
    c = core.calibration()
    assert c["amplitude_over_r_squared_lower"] > sp.Rational(2, 25)
    assert c["phase_separation_over_r_squared_lower"] > sp.Rational(3, 20)


def test_smoothing_has_room_inside_the_admissible_L1_deficit():
    c = core.calibration()
    assert c["smooth_pulse_deficit_over_r_upper"] == 4*c["smooth_transition_width_over_r"]
    assert c["smooth_pulse_deficit_over_r_upper"] < c["deficit_over_r"]


@pytest.mark.parametrize("momentum", [0, 1, 4])
def test_conditional_bound_has_actual_spatial_band_and_amplitude_scaling(momentum):
    bound = core.separation_bound(momentum_squared=momentum)
    assert bound == sp.Rational(3, 20*10**6)
    assert core.separation_bound(amplitude=3, momentum_squared=momentum) == 3*bound


@pytest.mark.parametrize("bad", [True, 0.001, sp.Float(".001"), 0, -1, sp.Rational(1, 100), sp.oo, sp.nan])
def test_outside_radius_or_inexact_input_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        core.calibration(bad)


@pytest.mark.parametrize("bad", [True, 0.0, -1, sp.Rational(1, 99)])
def test_profile_deficit_hypothesis_not_dropped(bad):
    with pytest.raises((TypeError, ValueError)):
        core.calibration(deficit=bad)


def test_invalid_domain_and_scope_controls_are_exact():
    controls = verify.controls()
    assert controls["rejected_inputs"] == 20
    assert controls["actual_subsequence_limits_not_assumed"] is True
    assert controls["conditional_numeric_bound_does_not_certify_a_profile"] is True


@pytest.mark.parametrize("bad", [True, 0.001, "1/1000", 0, -1, Fraction(1, 100)])
def test_separate_engine_rejects_invalid_radius(bad):
    with pytest.raises((TypeError, ValueError)):
        independent.calibration(bad)
