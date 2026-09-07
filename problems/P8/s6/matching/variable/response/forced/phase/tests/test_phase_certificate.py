"""Exact report/source/claim integrity with a fresh recursive parent replay."""

import copy
import json

import pytest
from p8_forced_phase import verify


@pytest.fixture(scope="module")
def report():
    return verify.build_report()


def test_all_sources_including_separate_audit_are_pinned(report):
    assert report["source_sha256"] == {str(p.relative_to(verify.ROOT)): verify.sha(p) for p in verify.source_files()}
    assert "tests/test_phase_independent_audit.py" in report["source_sha256"]
    assert report["prior_sha256"]["S6_26_actual_source_and_recursive_prepared_response_parent_contract"] == verify.FORCED_SHA


def test_full_error_and_precise_nonconvergence_scope(report):
    assert "1200r^4" in report["proof_chain"]["full_error"]
    assert "3r^2/20" in report["proof_chain"]["actual_separation"]
    assert "delta-independent" in report["unchanged_physical_problem"]["pulse"]
    assert any("delta-dependent" in value for value in report["not_established"])
    assert any("No original" in value for value in report["not_established"])
    assert len(report["independent_Fraction_bridges"]) == 9


@pytest.mark.parametrize("key", ["status", "proof_chain", "source_sha256", "not_established"])
def test_omitted_or_changed_scope_rejected(report, key):
    modified = copy.deepcopy(report)
    modified[key] = "DROPPED"
    with pytest.raises(ValueError):
        verify.validate_report(modified, report)


def test_unreplayed_extra_claim_rejected(report):
    modified = copy.deepcopy(report)
    modified["all_low_energy_reductions_excluded"] = True
    with pytest.raises(ValueError):
        verify.validate_report(modified, report)


def test_frozen_report_is_read_only(report):
    before = verify.REPORT.read_bytes()
    verify.validate_report(json.loads(before), report)
    assert verify.REPORT.read_bytes() == before
