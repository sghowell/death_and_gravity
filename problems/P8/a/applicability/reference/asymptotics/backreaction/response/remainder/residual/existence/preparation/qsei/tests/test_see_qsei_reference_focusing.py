"""Actual-source credit and the explicit remaining focusing obstruction."""

import pytest
import sympy as sp
from p8a_see_qsei import focusing, reference, sampling


def test_new_reference_credit_uses_actual_fixed_point_distance():
    data = reference.calibration()
    assert data["actual_fixed_point_distance"] < data["distance_cap"]
    assert data["numerator_loss_upper"] < data["new_numerator_lower"]
    assert data["positive_EED_dimensionless"] == sp.Rational(1, 4199040)
    # The larger full ball radius does not prove this positive-credit margin.
    full_ball_loss = data["numerator_loss_upper"]*sp.Rational(1, 10**6)/data["distance_cap"]
    assert full_ball_loss > data["new_numerator_lower"]


def test_reference_identity_and_renormalization_scale_dimensions():
    assert all(value == 0 for value in reference.identities().values())
    base = reference.positive_credit()
    assert reference.positive_credit(normalization=2) == base/16
    assert reference.positive_credit(eta_star=2) == base/256
    assert reference.positive_credit(hbar=3) == 3*base


@pytest.mark.parametrize("kwargs", [{"normalization": 0}, {"eta_star": -1}, {"hbar": -1},
                                    {"normalization": 1.0}, {"eta_star": True}])
def test_positive_credit_rejects_unphysical_or_inexact_scales(kwargs):
    with pytest.raises((TypeError, ValueError)):
        reference.positive_credit(**kwargs)


def test_focusing_index_is_strictly_out_of_reach_inside_new_slab():
    data = focusing.calibration()
    span = data["proper_span_upper"]
    assert data["index_times_T0_lower"] == focusing.index_lower(span)
    assert focusing.index_lower(span/2) > focusing.index_lower(span) > sp.Rational(3, 4)
    assert data["normal_jacobian_lower"] == sp.Rational(8, 27)
    assert data["Q2_over_available_duration_squared_lower"] == 160000000
    assert all(value > 0 for value in data["strict_margins"].values())


@pytest.mark.parametrize("duration", [0, -1, True, 0.1, sp.oo, sp.Rational(1, 10**9)])
def test_focusing_api_does_not_extrapolate_past_proved_domain(duration):
    with pytest.raises((TypeError, ValueError)):
        focusing.index_lower(duration)


def test_short_index_identities_and_trial_function_control():
    assert all(value == 0 for value in focusing.identities().values())
    tau = sampling.calibration()["proper_span_upper"]
    x = sp.Symbol("x", nonnegative=True)
    linear = 1-x/tau
    # Worst permitted Ricci, evaluated on a real admissible trial, is
    # still above the uniform lower bound and the contraction threshold.
    direct = sp.integrate(3*sp.diff(linear, x)**2-linear**2/4, (x, 0, tau))
    assert direct == 3/tau-tau/12
    assert direct > focusing.index_lower(tau) > sp.Rational(3, 4)
