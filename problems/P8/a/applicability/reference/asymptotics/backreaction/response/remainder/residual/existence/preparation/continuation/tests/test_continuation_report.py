"""The immutable source seal and explicitly limited claim scope."""

import copy
import json

import pytest
from p8a_continuation import independent, verify


@pytest.fixture(scope="module")
def actual_report():
    return verify.build_report()


def test_frozen_report_matches_full_read_only_replay(actual_report):
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual_report)


def test_source_hashes_include_independent_covariant_audit(actual_report):
    hashes = actual_report["source_sha256"]
    assert "tests/test_continuation_covariant_audit.py" in hashes
    assert actual_report["prior_sha256"] == {"A11": verify.PIN}
    for path, digest in hashes.items():
        assert verify.sha(verify.ROOT/path) == digest


def test_fraction_gate_rejects_changed_physical_amplitude():
    prior = json.loads(verify.prior.REPORT.read_text())
    corrupt = copy.deepcopy(prior)
    corrupt["derived_constants"]["full_map"]["delta"] = "1/1000"
    with pytest.raises(ValueError, match="margin"):
        independent.replay(corrupt)


def test_claim_cannot_be_mutated_into_actual_instability_or_closure(actual_report):
    for key, value in [("status", "P8_CLOSED"), ("prior_sha256", {})]:
        corrupt = copy.deepcopy(actual_report)
        corrupt[key] = value
        with pytest.raises(ValueError, match="certificate"):
            verify.validate_report(corrupt, actual_report)
    corrupt = copy.deepcopy(actual_report)
    corrupt["theorem"]["full_rolling_Frechet_derivative"] = True
    with pytest.raises(ValueError, match="certificate"):
        verify.validate_report(corrupt, actual_report)


def test_preserved_remainder_and_nonclosure_are_explicit(actual_report):
    theorem = actual_report["theorem"]
    assert len(theorem["whole_remaining_map_retained"]) == 5
    assert theorem["full_rolling_Frechet_derivative"] is False
    assert theorem["model_changed"] is False
    assert theorem["condition_imposed_on_actual_state"] is False
    assert "nonlinear runaway, instability, generic no-go or a necessary maximum step length" in actual_report["not_established"]
