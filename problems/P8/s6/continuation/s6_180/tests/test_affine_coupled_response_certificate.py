"""Read-only exact report, explicit source manifest and field mutation tests."""

import copy
import json

import pytest
from p8_vacuum_affine_coupled_response import verify

KEYS = [
    "schema",
    "claim",
    "date",
    "status",
    "prior_sha256",
    "source_sha256",
    "formulation",
    "written_proofs",
    "exact_residuals",
    "named_exact_check_count",
    "checked_scalar_entries",
    "proof_checks",
    "actual_current_parent_tree_and_uncancelled_time_channel",
    "matched_retained_inverse_and_physical_reconstruction",
    "observable_and_scope",
    "primitive_and_matching_frontier",
    "controls",
    "verdict",
    "not_established",
    "verification_boundary",
]


def test_complete_read_only_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_exact_source_and_field_manifest():
    d = verify.build_report()
    assert set(d) == set(KEYS)
    assert d["source_sha256"] == {
        str(p.relative_to(verify.ROOT)): verify.sha(p) for p in verify.source_files()
    }


@pytest.mark.parametrize("key", KEYS)
def test_every_report_field_mutation_rejected(key):
    d = verify.build_report()
    bad = copy.deepcopy(d)
    bad[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, d)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_or_extra_field_rejected(mutation):
    d = verify.build_report()
    bad = copy.deepcopy(d)
    if mutation == "missing":
        del bad["matched_retained_inverse_and_physical_reconstruction"]
    else:
        bad["self_consistent_quantum_background_complete"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, d)


def test_exact_counts_and_scope():
    d = verify.build_report()
    assert len(d["source_sha256"]) == 18
    assert d["named_exact_check_count"] == 109
    assert d["checked_scalar_entries"] == 140
    assert len(d["proof_checks"]) == 24
    assert len(d["controls"]) == 9
    assert d["controls"]["rejected_inputs"] == 101
    assert "NOT_SMALL_FULL_INVERSE_QUANTUM_BACKGROUND" in d["status"]
