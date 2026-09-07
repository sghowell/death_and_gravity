import copy
import json

import pytest
from p8_trimetric import verify


def test_frozen_auxiliary_report_matches_read_only_exact_replay():
    actual = verify.build_report()
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual)
    assert actual["claim"] == "P8-S6.13.AUXILIARY"
    assert sum(map(len, actual["exact_residuals"].values())) == 49
    assert actual["independent_Fraction_replay"]["full_matrix_first_jet_directions_checked"] == 64
    assert "tests/test_trimetric_covariant_audit.py" in actual["source_sha256"]


def test_report_rejects_extending_the_no_flat_verdict_across_beta4_countercontrol():
    actual = verify.build_report()
    bad = copy.deepcopy(actual)
    bad["separately_named_beta4_countercontrol"]["flat_checks"] = "No regular extended flat vacuum"
    with pytest.raises(ValueError):
        verify.validate_report(bad, actual)


def test_report_pins_context_and_the_additional_flat_vacuum_contract():
    assert verify.prior_checks() == {"S6_12_source_preserving_branch": verify.PRIOR_SHA,
                                     "adopted_S6_contract": verify.CONTRACT_SHA}


def test_checked_controls_distinguish_full_truncated_and_changed_actions():
    controls = verify.controls()
    original = controls["chosen_parent_no_flat_and_cancellation"]
    assert original["dropping_e_f_equations_false_vacuum_u_Euler_00"] == "0"
    assert original["original_regular_g_flat_Euler_00"] == "-2"
    assert original["actual_action_at_formal_cancel_over_B_h0_volume"] == "9/32"
    assert controls["separately_named_beta4_extension"]["positive_vacuum_mass_squared_G_F_q_1"] == "2"
