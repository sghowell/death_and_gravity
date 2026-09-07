import copy
import json

import pytest
from p8_trimetric_cones import verify


def test_frozen_report_matches_exact_read_only_replay():
    actual = verify.build_report()
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual)
    assert actual["claim"] == "P8-S6.14.AUXILIARY"
    assert sum(map(len, actual["exact_residuals"].values())) == 28
    assert "tests/test_trimetric_cones_covariant_audit.py" in actual["source_sha256"]


def test_replay_rejects_a_full_parent_to_all_light_EFT_overclaim():
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed["general_full_parent_cone_theorem"]["scope"] = "All light-only EFTs excluded"
    with pytest.raises(ValueError):
        verify.validate_report(changed, actual)


def test_checked_controls_keep_actual_and_algebraic_only_examples_distinct():
    c = verify.controls()
    assert c["background_hypotheses"]["mixed_sign_full_flat_relative_spring"] == "4"
    assert c["background_hypotheses"]["mixed_sign_auxiliary_c_e"] == "1/2"
    assert c["tensor_source_clock_and_symmetry"]["dropping_probe_contact_at_A2_q1_j1"] == "1"
    bridge = verify.clock_bridge()
    assert len(bridge) == 2
    assert all(value == "0" for fixture in bridge for value in fixture["primary_minus_independent"].values())


def test_prior_pin_replays_action_original_curvature_and_frozen_matter():
    assert verify.prior_checks() == {"S6_13_actual_auxiliary_parent": verify.PRIOR_SHA,
                                      "original_B_curvature_and_physical_cone_contract": verify.P8_FORMULATION_SHA,
                                      "conditional_original_CD_M1_free_chi": verify.CD_SHA}
