"""Immutable current clock and complete finite heavy local response report with every-field rejection."""

import copy
import json

import pytest
from p8_vacuum_affine_heavy_clock_quadratic import verify

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
    "complete_current_clock_quadratic_comparison",
    "complete_physical_scalar_and_tensor_local_response",
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
        del bad["complete_physical_scalar_and_tensor_local_response"]
    else:
        bad["quantum_limit_commutes"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_exact_counts_and_conditional_boundary():
    data = verify.build_report()
    assert len(data["source_sha256"]) == 17
    assert data["named_exact_check_count"] == 49
    assert data["checked_scalar_entries"] == 64
    assert len(data["proof_checks"]) == 41
    assert len(data["controls"]) == 9 and data["controls"]["rejected_inputs"] == 284
    assert (
        "CURRENT_CLASSICAL_COEFFICIENT_COMPARISON_AND_COMPLETE_FINITE_LOCAL_HEAVY_RESPONSE"
        in data["status"]
    )
    assert (
        "NOT_FULL_QUANTUM_RESPONSE_INVERSE_INTERACTING_BOUNCE_UV_REGGE_ORIGINAL_V_G_B_OR_P8"
        in data["status"]
    )
    assert (
        data["prior_sha256"][
            "all_frozen_original_and_separate_model_scientific_bytes_unchanged"
        ]
        is True
    )
