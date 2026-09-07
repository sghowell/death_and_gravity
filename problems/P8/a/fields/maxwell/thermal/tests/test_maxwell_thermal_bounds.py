"""Continuous rational C3 history embedding, including all endpoints."""

import sympy as sp
from p8a_maxwell_cosmology import geometry
from p8a_maxwell_thermal import bounds, dynamics


def test_every_scheme_response_and_reciprocal_identity():
    assert all(sp.simplify(value) == 0 for value in bounds.identities().values())


def test_actual_history_is_strictly_inside_the_pinned_anchored_radiation_tube():
    data = bounds.history(dynamics.DELTA_MAX)
    assert data["C3_error_upper"] == [sp.Rational(31, 28125000000), sp.Rational(93, 781250000),
                                      sp.Rational(3131, 292968750), sp.Rational(149048, 146484375)]
    assert min(data["strict_A18_tube_margins"]) > 0
    assert all(error < cap for error, cap in zip(data["C3_error_upper"], geometry.ERROR_CAPS, strict=True))
    assert data["anchored_at_observer"] and data["actual_state_and_SEE_history_realized"]
    assert not data["A18_future_bounds_proved_beyond_this_branch_endpoint"]


def test_smaller_delta_uniformly_improves_each_continuous_error_bound():
    large = bounds.history(dynamics.DELTA_MAX)
    small = bounds.history(dynamics.DELTA_MAX/10)
    assert small["C3_error_upper"] == [error/10 for error in large["C3_error_upper"]]


def test_numerator_controls_are_strict_and_no_denominator_is_zero():
    data = bounds.calibration()
    assert min(data["numerator_lower_values_at_one_quarter"]) > 0
    assert data["strict_named_coupling_margin"] > 0
    assert data["strict_backward_y_above_one"] > 0
