"""Read-only complete ancestry, exact source hashes and every-field mutations."""

import copy
import json

import pytest
from p8_vacuum_gauge_yukawa_screen import verify

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
    "new_candidate_and_one_loop_marginal_flow",
    "leading_heavy_fermion_threshold",
    "leading_matched_operator_gauge_cut",
    "prospective_boundary_calibration_not_full_matching",
    "controls",
    "verdict",
    "not_established",
    "verification_boundary",
]


def test_complete_read_only_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_every_field_and_exact_source_manifest():
    d = verify.build_report()
    assert set(d) == set(KEYS)
    assert d["source_sha256"] == {
        str(p.relative_to(verify.ROOT)): verify.sha(p) for p in verify.source_files()
    }


@pytest.mark.parametrize("key", KEYS)
def test_every_report_field_mutation_rejected(key):
    d = verify.build_report()
    changed = copy.deepcopy(d)
    changed[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, d)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_or_extra_field_rejected(mutation):
    d = verify.build_report()
    changed = copy.deepcopy(d)
    if mutation == "missing":
        del changed["leading_matched_operator_gauge_cut"]
    else:
        changed["original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, d)


def test_exact_counts_and_scope():
    d = verify.build_report()
    assert d["named_exact_check_count"] == 132
    assert d["checked_scalar_entries"] == 248
    assert len(d["proof_checks"]) == 38
    assert len(d["source_sha256"]) == 18
    assert d["controls"]["rejected_inputs"] == 16
    assert "NOT_FULL_FINITE_MASS_MATCHING" in d["status"]
