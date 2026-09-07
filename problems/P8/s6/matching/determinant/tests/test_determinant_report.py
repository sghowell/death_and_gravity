import copy
import json

import pytest
from p8_determinant import verify


@pytest.fixture(scope="module")
def actual():
    return verify.build_report()


def test_determinant_frozen_report_replays_every_hash_and_ancestor(actual):
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual)


def test_determinant_source_pins_exclude_mutable_tree_and_include_reserved_audit(actual):
    pins = actual["prior_context_sha256"]
    assert pins["S6_16_frozen_STAR_and_arithmetic_lineage"] == verify.PRIOR_SHA
    assert pins["S6_13_literal_auxiliary_action"] == verify.ACTION_SHA
    assert pins["adopted_S6_matching_contract"] == verify.CONTRACT_SHA
    sources = actual["source_sha256"]
    assert "tests/test_determinant_independent_audit.py" in sources
    assert not any("tree" in path.lower() for path in sources)
    assert all(verify.sha(verify.ROOT/path) == digest for path, digest in sources.items())


def test_determinant_report_retains_all_scope_boundaries(actual):
    assert "S_a" in actual["fixed_action_domain"]["regular_sum"]
    assert "S_N" in actual["fixed_action_domain"]["regular_sum"]
    assert "not" in actual["auxiliary_source_preserving_identity"]["chart"]
    assert actual["checked_omission_controls"]["singular_time_sum"]["actual_solution"] is False
    assert actual["checked_omission_controls"]["invalid_domain_calls_rejected"] == 8
    assert actual["independent_Fraction_replay"]["full_coframe_first_jet_directions"] == 192
    assert "OPEN" in actual["status"]


@pytest.mark.parametrize("key", ["theorem", "not_established", "fixed_action_domain", "source_sha256"])
def test_determinant_claim_or_hypothesis_omission_is_detected(actual, key):
    changed = copy.deepcopy(actual)
    del changed[key]
    with pytest.raises(ValueError):
        verify.validate_report(changed, actual)


def test_determinant_fake_singular_solution_promotion_is_detected(actual):
    changed = copy.deepcopy(actual)
    changed["checked_omission_controls"]["singular_time_sum"]["actual_solution"] = True
    with pytest.raises(ValueError):
        verify.validate_report(changed, actual)


def test_determinant_source_hash_and_identity_corruption_are_detected_read_only(actual):
    for key in ("source_sha256", "exact_residuals"):
        changed = copy.deepcopy(actual)
        changed[key] = {}
        with pytest.raises(ValueError):
            verify.validate_report(changed, actual)
