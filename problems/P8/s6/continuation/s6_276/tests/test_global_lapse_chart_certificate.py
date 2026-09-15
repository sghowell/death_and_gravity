"""Whole immutable report contract, not only selected JSON fields."""

import copy
import json

import pytest
from p8_vacuum_affine_global_lapse_chart_obstruction import verify

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
    "whole_joint_auxiliary_atlas_and_exact_source_fold",
    "whole_constraint_realization_and_connected_sheet_obstruction",
    "observable_and_scope",
    "primitive_and_matching_frontier",
    "controls",
    "verdict",
    "not_established",
    "verification_boundary",
]


def test_entire_native_report_read_only_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_entire_source_manifest_and_ordered_twenty_field_contract():
    data = verify.build_report()
    assert list(data) == KEYS
    assert data["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path)
        for path in verify.source_files()
    }


@pytest.mark.parametrize("key", KEYS)
def test_mutation_of_every_report_field_is_rejected(key):
    data = verify.build_report()
    bad = copy.deepcopy(data)
    bad[key] = "not_the_certified_value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_or_extra_field_is_rejected(mutation):
    data = verify.build_report()
    bad = copy.deepcopy(data)
    if mutation == "missing":
        del bad["whole_joint_auxiliary_atlas_and_exact_source_fold"]
    else:
        bad["original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_exact_counts_and_scope_are_complete():
    data = verify.build_report()
    assert len(data["source_sha256"]) == 17
    assert (data["named_exact_check_count"], data["checked_scalar_entries"]) == (
        97,
        158,
    )
    assert len(data["proof_checks"]) == 48
    assert len(data["controls"]) == 8 and data["controls"]["rejected_inputs"] == 72
    assert "NOT_QUANTUM_OR_ORIGINAL_V_G_B_P8_NO_GO" in data["status"]
    assert data["prior_sha256"]["S275_corrected_finite_scope_not_refuted"] is True
    assert (
        data["prior_sha256"]["S269_to_S274_physical_qualifications_unchanged"] is True
    )
    assert len(data["observable_and_scope"]["historical_qualification"]) == 6
