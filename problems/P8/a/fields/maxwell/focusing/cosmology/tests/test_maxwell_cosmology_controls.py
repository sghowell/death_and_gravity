"""Nonvacuity is proved with a full interval of complete smooth geometries."""

import pytest
import sympy as sp
from p8a_maxwell_cosmology import controls, geometry


def test_raw_smooth_joins_and_derivative_factor_identities():
    assert all(sp.simplify(value) == 0 for value in controls.identities().values())


def test_complete_smooth_family_has_uniform_jet_and_actual_tube_margins():
    data = controls.complete_geometry()
    assert data["raw_global_jet_caps_through_four"] == [4, 102, 8352, 630528, 887218176]
    assert min(data["strict_past_tube_margins"]) > 0
    assert min(data["future_cap_margins_from_global_bounds"]) >= 0
    assert data["future_timelike_and_null_complete"]
    assert not data["actual_quantum_SEE_solution_asserted"]
    assert data["covers_every_reference_power_in_interval"]


def test_original_tight_future_candidate_is_rejected_without_quantum_input():
    data = controls.rejected_tight_future()
    assert data["strict_incompatibility_gap"] == sp.Rational(4707907, 62554080)
    assert not data["QSEI_used_in_this_exclusion"]
    assert geometry.FUTURE_CAPS != geometry.EXPANDED_CAPS


def test_smoothing_width_must_actually_preserve_all_four_history_errors():
    with pytest.raises(ValueError):
        controls.complete_geometry(sp.Rational(1, 1000))


def test_all_four_global_reference_caps_follow_from_denominator_and_power_bounds():
    for j, cap in enumerate(controls.raw_global_jet_caps()):
        assert cap >= sp.factorial(j)*2**(2*j+1)


def test_named_reference_tube_is_not_an_initial_SEC_violation_example():
    data = geometry.calibration()
    assert data["this_tube_Ricci_upper_times_tau_squared"] < 0
    assert not data["this_tube_demonstrates_initial_SEC_violation"]
