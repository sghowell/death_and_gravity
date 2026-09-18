"""Whole twenty-field equality and frozen-source contract."""

import copy
import json

import pytest
from p8_vacuum_affine_leading_cloud_logarithmic_coefficient import verify

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
    "whole_original_source_and_defined_leading_cloud",
    "whole_conditioned_complete_coefficient_moments",
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
        del bad["whole_conditioned_complete_coefficient_moments"]
    else:
        bad["original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_counts_and_historical_scope():
    from p8_vacuum_affine_leading_cloud_logarithmic_coefficient import audit

    data = verify.build_report()
    assert len(data["source_sha256"]) == 17
    assert data["named_exact_check_count"] == len(audit.residuals()) == 257
    assert data["checked_scalar_entries"] == audit.scalar_entry_count() == 257
    assert len(data["proof_checks"]) == len(audit.gates()) == 30
    assert len(data["controls"]) == 8 and data["controls"]["rejected_inputs"] == 93
    assert "NOT_INTERACTING_ALL_N_V_G_B_P8_CLOSURE" in data["status"]
    assert data["prior_sha256"]["S279_remains_rejected_and_archived"] is True
    assert data["prior_sha256"]["S275_S276_S277_and_scoped_P8a_unchanged"] is True
    assert len(data["observable_and_scope"]["historical_qualification"]) == 6
