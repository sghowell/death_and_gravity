"""Read-only lineage, source seal and unmodified-source claim controls."""

import copy
import json

import pytest
from p8a_extension import independent, verify


@pytest.fixture(scope="module")
def actual_report():
    return verify.build_report()


def test_stored_extension_report_replays_exactly(actual_report):
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual_report)


def test_pins_include_the_reserved_independent_covariant_audit(actual_report):
    hashes = actual_report["source_sha256"]
    assert "tests/test_extension_covariant_audit.py" in hashes
    assert actual_report["prior_sha256"] == verify.PINS
    for relative, digest in hashes.items():
        assert verify.sha(verify.ROOT/relative) == digest


def test_independent_replay_rejects_changed_inverse_or_amplitude():
    old = json.loads(verify.original.REPORT.read_text())
    prior = json.loads(verify.prior.REPORT.read_text())
    for key, value in [("weighted_norm_upper", "100"), ("delta", "1/1000")]:
        corrupt = copy.deepcopy(prior)
        corrupt["derived_constants"]["rational_gate"][key] = value
        with pytest.raises(ValueError):
            independent.replay(old, corrupt)


def test_changed_source_or_automatic_QSEI_claim_cannot_replay(actual_report):
    corrupt = copy.deepcopy(actual_report)
    corrupt["negative_controls"]["cutoff_rescaled_to_L"] = True
    with pytest.raises(ValueError, match="certificate"):
        verify.validate_report(corrupt, actual_report)
    corrupt = copy.deepcopy(actual_report)
    corrupt["negative_controls"]["old_QSEI_or_positive_reference_credit_automatically_transferred"] = True
    with pytest.raises(ValueError, match="certificate"):
        verify.validate_report(corrupt, actual_report)


def test_actual_extension_and_remaining_scope_are_both_explicit(actual_report):
    theorem = actual_report["theorem"]
    assert theorem["actual_full_density_and_trace_SEE"] is True
    assert theorem["agrees_with_A11_on_original_slab"] is True
    assert theorem["smooth_on_the_same_longer_slab"] is True
    assert "completion of P8a or P8" in actual_report["not_established"]
