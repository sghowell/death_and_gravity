"""Smooth future-complete geometric countercontrol is not a quantum solution."""

import sympy as sp
from p8a_maxwell_focusing import controls


def test_every_polynomial_mollification_input_is_exact():
    assert all(sp.expand(value) == 0 for value in controls.identities().values())


def test_positive_convolution_preserves_strict_geometric_caps_and_history():
    data = controls.calibration()
    assert all(margin > 0 for margin in data["cap_margins"])
    assert data["history_margin"] > 0
    assert data["second_switch_derivative_square_margin"] > 0
    assert data["future_static_by_time_over_tau"] == sp.Rational(1001, 1000)


def test_geometry_countercontrol_does_not_become_a_Maxwell_SEE_witness():
    data = controls.calibration()
    assert data["future_timelike_and_null_complete_geometry"] is True
    assert data["actual_allowed_Maxwell_SEE_solution"] is False
    assert data["required_sigma_exceeds_one_by"] > 6
    assert data["initial_pointwise_SEC_holds"] is False
