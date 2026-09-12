"""Immutable native report, exact manifest and every-field mutation rejection."""

import copy
import json

import pytest
from p8_vacuum_affine_homogeneous_shear_inverse import verify

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
    "actual_homogeneous_leading_pair_full_dimensional_matching_and_original_local_tree",
    "actual_complete_weak_log_remainder_and_ordered_prepared_inverse",
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
        del bad["actual_complete_weak_log_remainder_and_ordered_prepared_inverse"]
    else:
        bad["actual_full_quantum_inverse"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_exact_counts_and_uncompleted_actual_quantum_boundary():
    data = verify.build_report()
    assert len(data["source_sha256"]) == 17
    assert (
        data["named_exact_check_count"] == 28 and data["checked_scalar_entries"] == 38
    )
    assert len(data["proof_checks"]) == 24
    assert len(data["controls"]) == 9 and data["controls"]["rejected_inputs"] == 173
    assert "NOT_NONZERO_TRANSFER_FULL_S222" in data["status"]
    assert (
        data["prior_sha256"][
            "original_state_prescription_and_all_frozen_scientific_bytes_unchanged"
        ]
        is True
    )
