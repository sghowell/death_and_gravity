"""Whole twenty-field native equality and immutable source contract."""

import copy
import json

import pytest
from p8_vacuum_affine_calorimetric_soft_resummation import verify

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
    "whole_original_source_and_total_energy_sum",
    "whole_Gamma_conversion_and_detector_limits",
    "observable_and_scope",
    "primitive_and_matching_frontier",
    "controls",
    "verdict",
    "not_established",
    "verification_boundary",
]


def test_entire_native_report_read_only_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_entire_source_manifest_and_twenty_field_contract():
    data = verify.build_report()
    assert list(data) == KEYS
    assert data["source_sha256"] == {
        str(p.relative_to(verify.ROOT)): verify.sha(p) for p in verify.source_files()
    }


@pytest.mark.parametrize("key", KEYS)
def test_every_report_field_mutation_rejected(key):
    data = verify.build_report()
    bad = copy.deepcopy(data)
    bad[key] = "not_the_certified_value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_and_extra_fields_rejected(mutation):
    data = verify.build_report()
    bad = copy.deepcopy(data)
    if mutation == "missing":
        del bad["whole_Gamma_conversion_and_detector_limits"]
    else:
        bad["original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_exact_counts_and_historical_scope():
    data = verify.build_report()
    assert len(data["source_sha256"]) == 18
    assert (data["named_exact_check_count"], data["checked_scalar_entries"]) == (
        223,
        223,
    )
    assert len(data["proof_checks"]) == 37
    assert len(data["controls"]) == 8 and data["controls"]["rejected_inputs"] == 73
    assert "NOT_FULL_V_G_B_P8_CLOSURE" in data["status"]
    assert data["prior_sha256"]["S279_remains_rejected_and_archived"] is True
    assert data["prior_sha256"]["S275_S276_S277_and_scoped_P8a_unchanged"] is True
    assert len(data["observable_and_scope"]["historical_qualification"]) == 6
