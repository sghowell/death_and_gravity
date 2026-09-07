import copy
import json

import pytest
from p8_hr_tree import verify


def test_frozen_certificate_replays_without_writes():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_manifest_and_explicit_ancestor_pins_are_replayed():
    report = verify.build_report()
    assert len(report["prior_context_sha256"]) == 5
    assert report["prior_context_sha256"]["S6_16_constant_star_action_and_controls"] == verify.STAR_SHA
    assert report["prior_context_sha256"]["S6_5_bimetric_monotonicity_context"] == verify.BIMETRIC_SHA
    for relative, expected in report["source_sha256"].items():
        assert verify.sha(verify.ROOT/relative) == expected
    assert "tests/test_hr_tree_independent_audit.py" in report["source_sha256"]
    assert sum(map(len, report["exact_residuals"].values())) == 60


def test_promoting_the_result_beyond_its_physical_target_domain_is_rejected():
    actual = verify.build_report()
    false = copy.deepcopy(actual)
    false["action_and_domain"]["Einstein_coefficients"] = "Any physical target, including disconnected zero-EH target"
    with pytest.raises(ValueError):
        verify.validate_report(false, actual)
    false = copy.deepcopy(actual)
    false["status"] = "P8_COMPLETELY_CLOSED"
    with pytest.raises(ValueError):
        verify.validate_report(false, actual)


def test_independent_bridge_counts_and_exact_outputs():
    bridges = verify.independent_bridges()
    assert len(bridges["literal_edge_stresses"]) == 5
    assert sum(len(item["primary_minus_independent"]) for item in bridges["literal_edge_stresses"]) == 35
    assert len(bridges["Gaussian_vs_leaf_incidence"]) == 4
    assert sum(len(item["primary_minus_independent_fluxes"]) for item in bridges["Gaussian_vs_leaf_incidence"]) == 12
    assert len(bridges["fixed_component_first_jets"]) == 5
    assert sum(len(item["primary_minus_independent"]) for item in bridges["fixed_component_first_jets"]) == 40


def test_omission_controls_and_strict_endpoint_equality_guard():
    controls = verify.checked_controls()
    assert controls["wrong_weight_and_normalization_defects"] == {
        "wrong_one_lapse_Bianchi_weight": "-12", "wrong_two_lapse_null_weight": "4", "raw_J_normalization_omission": "1"}
    assert controls["actual_algebraic_active_component"] == [0, 2]
    assert controls["actual_zero_target_Hprime_at_zero"] == 2
    assert controls["actual_nonNEC_null_at_bounce"] == -4
    assert controls["invalid_domain_calls_rejected"] == 15
    assert controls["CD_endpoint_equality_control"]["status"] == "INCONCLUSIVE"
