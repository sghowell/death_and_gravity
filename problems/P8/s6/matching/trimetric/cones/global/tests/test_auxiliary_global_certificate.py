import copy
import json

import pytest
from p8_trimetric_global import verify


def test_frozen_report_matches_exact_recursive_read_only_replay():
    actual = verify.build_report()
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual)
    assert actual["claim"] == "P8-S6.15.AUXILIARY"
    assert sum(map(len, actual["exact_residuals"].values())) == 58
    assert "tests/test_auxiliary_global_independent_audit.py" in actual["source_sha256"]
    assert len(actual["source_sha256"]) == 16


def test_exact_dynamics_and_actual_solution_bridges_use_two_separate_engines():
    phase, rolling = verify.dynamics_bridge(), verify.rolling_bridge()
    assert len(phase) == 5
    assert len(rolling) == 4
    assert {entry["q"] for entry in rolling} == {"1", "-1"}
    assert all(value == "0" for entry in phase+rolling for value in entry["primary_minus_independent"].values())


def test_controls_reject_false_vacuum_and_keep_singular_background_scope():
    checked = verify.controls()
    guards = checked["actual_solution_and_domain_guards"]
    assert guards["positive_cap_does_not_supply_full_vacuum_E_u00"] == "-2"
    assert guards["singular_TT_hessian"] == "0"
    assert guards["all_zero_links_arbitrary_metric_H_prime_at_zero"] == "2"


def test_replay_rejects_scope_promotion_to_all_parent_actions():
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed["primary_full_background_theorem"]["conclusion"] = "Every parent action is excluded"
    with pytest.raises(ValueError):
        verify.validate_report(changed, actual)


def test_prior_pins_keep_actual_action_target_and_adopted_S6_separate():
    assert verify.prior_checks() == {
        "S6_14_actual_auxiliary_action_and_cones": verify.PRIOR_SHA,
        "original_B_curvature_and_physical_metric": verify.P8_FORMULATION_SHA,
        "conditional_original_CD_M1_endpoint_target": verify.CD_SHA,
        "adopted_S6_controlled_matching_contract": verify.S6_FORMULATION_SHA,
    }
