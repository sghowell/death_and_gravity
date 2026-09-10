"""Complete ancestry replay, pinned sources and every-field mutation controls."""

import copy
import json

import pytest
from p8_vacuum_two_loop_double_bubble import verify

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
    "exact_disjoint_refinement_selection",
    "every_actual_restricted_forest_and_local_subtraction",
    "complete_finite_bubble_and_heavy_triangle_factors",
    "actual_integrated_group_and_combined_partial_bound",
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
def test_missing_or_extra_report_field_rejected(mutation):
    d = verify.build_report()
    changed = copy.deepcopy(d)
    if mutation == "missing":
        del changed["actual_integrated_group_and_combined_partial_bound"]
    else:
        changed["original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, d)


def test_exact_counts_and_original_scope():
    d = verify.build_report()
    assert d["named_exact_check_count"] == d["checked_scalar_entries"] == 142
    assert len(d["proof_checks"]) == 34
    assert len(d["source_sha256"]) == 18
    assert d["controls"]["rejected_inputs"] == 26
    assert (
        "NOT_REMAINING_SIXTEEN_WINEGLASS_FINITE_POTENTIAL_INSERTIONS_COMPLETE_TWO_LOOP_POLE_LSZ"
        in d["status"]
    )
