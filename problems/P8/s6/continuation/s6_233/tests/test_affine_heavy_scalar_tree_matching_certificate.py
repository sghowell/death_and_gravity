"""Immutable separate two-scalar tree matching report and every-field mutation rejection."""

import copy
import json

import pytest
from p8_vacuum_affine_heavy_scalar_tree_matching import verify

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
    "separate_local_two_scalar_classical_model_and_potential",
    "complete_original_tree_remainder_and_full_first_elastic_comparison",
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
        del bad["complete_original_tree_remainder_and_full_first_elastic_comparison"]
    else:
        bad["physical_cutoff"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_exact_counts_and_conditional_boundary():
    data = verify.build_report()
    assert len(data["source_sha256"]) == 17
    assert (
        data["named_exact_check_count"] == 42 and data["checked_scalar_entries"] == 42
    )
    assert len(data["proof_checks"]) == 23
    assert len(data["controls"]) == 9 and data["controls"]["rejected_inputs"] == 213
    assert (
        "SEPARATE_V2S_T1_CLASSICAL_MODEL_AND_COMPLETE_TREE_MATCHING" in data["status"]
    )
    assert "NOT_FULL_QUANTUM_ORIGINAL_AFFINE_PARENT" in data["status"]
    assert (
        data["prior_sha256"][
            "original_state_prescription_and_all_frozen_scientific_bytes_unchanged"
        ]
        is True
    )
