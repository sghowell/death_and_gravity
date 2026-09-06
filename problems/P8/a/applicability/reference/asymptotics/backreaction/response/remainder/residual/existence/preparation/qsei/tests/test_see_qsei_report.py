"""Read-only ancestor replay and independent arithmetic/source seal."""

import copy
import json

import pytest
from p8a_see_qsei import independent, verify


@pytest.fixture(scope="module")
def actual_report():
    return verify.build_report()


def test_stored_certificate_replays_exactly(actual_report):
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual_report)


def test_all_sources_and_reserved_independent_audit_are_pinned(actual_report):
    hashes = actual_report["source_sha256"]
    assert "tests/test_see_qsei_covariant_audit.py" in hashes
    for relative, digest in hashes.items():
        assert verify.sha(verify.ROOT/relative) == digest
    assert actual_report["prior_sha256"] == verify.PINS


def test_fraction_replay_rejects_state_distance_corruption():
    prior = json.loads(verify.prior.REPORT.read_text())
    corrupt = copy.deepcopy(prior)
    corrupt["derived_constants"]["full_map"]["actual_gate"]["fixed_point_distance"] = "1/1000000"
    with pytest.raises(ValueError, match="margin"):
        independent.replay(corrupt)


def test_mutated_bound_or_closed_scope_does_not_replay(actual_report):
    for key, value in [("status", "P8_CLOSED"), ("prior_sha256", {})]:
        corrupt = copy.deepcopy(actual_report)
        corrupt[key] = value
        with pytest.raises(ValueError, match="certificate"):
            verify.validate_report(corrupt, actual_report)
    corrupt = copy.deepcopy(actual_report)
    corrupt["derived_constants"]["proper_H2_sampling"]["rounded_coefficient"] = "1"
    with pytest.raises(ValueError, match="certificate"):
        verify.validate_report(corrupt, actual_report)


def test_nonclosure_and_target_state_scope_are_explicit(actual_report):
    theorem = actual_report["theorem"]
    assert "all Hadamard" in theorem["target_states"]
    assert theorem["coefficient_is_globally_optimal"] is False
    assert theorem["higher_potential_jet_norms_assumed"] is False
    assert theorem["prior_A9_numerical_bound_transferred"] is False
    assert actual_report["focusing_result"]["sufficient_index_trigger_available"] is False
    assert "completion of original P8a or P8" in actual_report["not_established"]
