"""Whole twenty-field report equality, complete immutable source manifest."""

import copy
import json

import pytest
from p8_vacuum_affine_triangle_curvature_coefficient import verify

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
    "whole_original_lift_and_triangle_coefficient",
    "whole_complete_dressing_and_original_calibrations",
    "observable_and_scope",
    "primitive_and_matching_frontier",
    "controls",
    "verdict",
    "not_established",
    "verification_boundary",
]


def test_entire_native_read_only_report():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_entire_manifest_and_twenty_fields():
    data = verify.build_report()
    assert list(data) == KEYS
    assert data["source_sha256"] == {
        str(p.relative_to(verify.ROOT)): verify.sha(p) for p in verify.source_files()
    }


@pytest.mark.parametrize("key", KEYS)
def test_every_field_mutation_rejected(key):
    data = verify.build_report()
    bad = copy.deepcopy(data)
    bad[key] = "not_the_certified_value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_extra_fields_rejected(mutation):
    data = verify.build_report()
    bad = copy.deepcopy(data)
    if mutation == "missing":
        del bad["whole_complete_dressing_and_original_calibrations"]
    else:
        bad["original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_counts_and_historical_scope():
    from p8_vacuum_affine_triangle_curvature_coefficient import audit

    data = verify.build_report()
    assert len(data["source_sha256"]) == 21
    assert data["named_exact_check_count"] == len(audit.residuals()) == 514
    assert data["checked_scalar_entries"] == audit.scalar_entry_count() == 536
    assert len(data["proof_checks"]) == len(audit.gates()) == 42
    assert len(data["controls"]) == 8 and data["controls"]["rejected_inputs"] == 670
    assert "NOT_FULL_HARD_V_G_B_P8_CLOSURE" in data["status"]
    assert data["prior_sha256"]["S279_remains_rejected_and_archived"] is True
    assert data["prior_sha256"]["S275_S276_S277_and_scoped_P8a_unchanged"] is True
    assert len(data["observable_and_scope"]["historical_qualification"]) == 6
