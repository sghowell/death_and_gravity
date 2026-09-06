import copy
import json

import pytest
from p8_composite_asymptotics import verify


def test_certificate_replays_sources_and_full_immutable_lineage():
    actual = verify.build_report()
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual)
    assert actual["prior_S6_9_composite_light_sha256"] == verify.PRIOR_SHA
    assert actual["claim"] == "P8-S6.10.COMPOSITE"


@pytest.mark.parametrize("key", ["domain", "actual_limiting_homogeneous_response",
                               "fixed_interval_convergence", "singular_limit_and_amplitude_guard",
                               "not_established"])
def test_changed_scope_or_finite_Y_instability_inflation_fails_replay(key):
    actual = verify.build_report()
    altered = copy.deepcopy(actual)
    altered[key] = "uniform finite-amplitude nonlinear instability at all finite Y"
    with pytest.raises(ValueError, match="differs from exact replay"):
        verify.validate_report(altered, actual)


def test_normalization_state_direction_and_amplitude_controls_are_replayed():
    actual = verify.build_report()
    controls = actual["omission_amplitude_and_direction_controls"]
    assert controls["true_initial_heavy_acceleration"] == "55/2"
    assert controls["omitted_normalization_false_acceleration"] == "-2"
    assert controls["backward_z_below_one_margin"] == "13/90"
    assert controls["small_metric_amplitude_absolute_projection_limit"] == "0"
    assert controls["zero_heavy_data_unforced_solution"] == "0"
