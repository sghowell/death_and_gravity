"""Immutable native report, exact manifest and every-field mutation rejection."""

import copy
import json

import pytest
from p8_vacuum_affine_full_finite_local_reference import verify

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
    "actual_complete_local_Hessian_factorization_and_coefficients",
    "source_time_reciprocity_conjugated_remainder_and_full_finite_inverse",
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
        del bad["source_time_reciprocity_conjugated_remainder_and_full_finite_inverse"]
    else:
        bad["actual_full_quantum_inverse"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_exact_counts_and_uncompleted_actual_quantum_boundary():
    data = verify.build_report()
    assert len(data["source_sha256"]) == 17
    assert (
        data["named_exact_check_count"] == 31 and data["checked_scalar_entries"] == 33
    )
    assert len(data["proof_checks"]) == 24
    assert len(data["controls"]) == 9 and data["controls"]["rejected_inputs"] == 184
    assert "NOT_FULL_NONLOCAL_CURVED_QUANTUM_SYSTEM" in data["status"]
    assert (
        data["prior_sha256"][
            "original_state_prescription_and_all_frozen_scientific_bytes_unchanged"
        ]
        is True
    )
