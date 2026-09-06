"""Rational bracket and norm estimates independent of rounded pole values."""

import json

import pytest
import sympy as sp
from p8a_continuation import bounds, independent
from p8a_preparation import verify as prior


def test_fraction_engine_replays_every_rational_constant():
    assert independent.serialize(bounds.calibration()) == independent.replay(json.loads(prior.REPORT.read_text()))


def test_exact_brackets_and_large_norms():
    data = bounds.calibration()
    assert data["positive_pole_lower"] == 10**6
    assert data["positive_pole_upper"] == 2*10**6
    assert data["norm_lower_at_10_to_minus_5"] > 1000
    assert data["norm_lower_at_10_to_minus_4"] > 10**28
    assert data["weighted_norm_upper"] == sp.Rational(240, 769)
    assert all(value > 0 for value in data["strict_margins"].values())


def test_coarse_pole_brackets_are_not_the_isolated_logarithmic_pole():
    data = bounds.calibration()
    assert data["positive_pole_lower"] > 100000*data["isolated_log_pole_upper"]
    assert data["A11_existing_length"] < sp.Rational(1, 10**5)


def test_rational_exponential_lower_and_duration_domain():
    assert bounds.integer_duration_norm_lower(100) > bounds.integer_duration_norm_lower(10)
    # A negative coarse lower estimate for short duration is not clipped
    # into a false claim of a negative operator norm.
    assert bounds.integer_duration_norm_lower(1) < 0
    for value in (0, -1, True, 1.0, sp.Rational(3, 2), sp.oo):
        with pytest.raises((TypeError, ValueError)):
            bounds.integer_duration_norm_lower(value)


def test_damped_pair_budget_and_positive_residue_are_both_needed():
    data = bounds.calibration()
    assert data["positive_step_residue_lower"] == sp.Rational(1, 17)
    assert data["damped_pair_step_absolute_cap"] == sp.Rational(8, 3)
    direct = ((sp.Rational(8, 3)**10-1)*data["positive_step_residue_lower"]
              -data["damped_pair_step_absolute_cap"])
    assert direct == data["norm_lower_at_10_to_minus_5"]
    assert direct < direct+data["damped_pair_step_absolute_cap"]


def test_weighted_gain_does_not_control_the_unweighted_norm():
    data = bounds.calibration()
    assert data["weighted_sigma"]*sp.Rational(1, 10**4) == 200
    assert data["weighted_norm_upper"]*2**200 > 10**50
