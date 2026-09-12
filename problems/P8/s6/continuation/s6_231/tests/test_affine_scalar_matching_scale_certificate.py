"""Immutable scalar matching report and every-field mutation rejection."""

import copy
import json

import pytest
from p8_vacuum_affine_scalar_matching_scale import verify

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
    "complete_original_tree_identical_partial_waves_and_first_elastic_coefficient",
    "necessary_unitarity_correction_and_conditional_physical_matching_tradeoff",
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
    data = verify.build_report()
    assert set(data) == set(KEYS)
    assert data["source_sha256"] == {
        str(p.relative_to(verify.ROOT)): verify.sha(p) for p in verify.source_files()
    }


@pytest.mark.parametrize("key", KEYS)
def test_every_report_field_mutation_rejected(key):
    data = verify.build_report()
    bad = copy.deepcopy(data)
    bad[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_and_extra_fields_rejected(mutation):
    data = verify.build_report()
    bad = copy.deepcopy(data)
    if mutation == "missing":
        del bad[
            "necessary_unitarity_correction_and_conditional_physical_matching_tradeoff"
        ]
    else:
        bad["physical_cutoff"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_exact_counts_and_conditional_boundary():
    data = verify.build_report()
    assert len(data["source_sha256"]) == 17
    assert (
        data["named_exact_check_count"] == 34 and data["checked_scalar_entries"] == 34
    )
    assert len(data["proof_checks"]) == 21
    assert len(data["controls"]) == 9 and data["controls"]["rejected_inputs"] == 223
    assert "CONDITIONAL_MATCHING_DISJUNCTION" in data["status"]
    assert "NOT_PHYSICAL_CUTOFF_QUANTUM_DECOUPLING" in data["status"]
    assert (
        data["prior_sha256"][
            "original_state_prescription_and_all_frozen_scientific_bytes_unchanged"
        ]
        is True
    )
