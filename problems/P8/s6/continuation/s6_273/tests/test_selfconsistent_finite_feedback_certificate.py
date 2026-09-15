"""Immutable whole-source and complete-report finite hybrid contract."""

import copy
import json

import pytest
from p8_vacuum_affine_selfconsistent_finite_feedback import verify

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
    "whole_source_canonical_force_and_actual_reference_symmetry",
    "whole_full_spatial_average_operator_bounds_and_coupled_turnaround",
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
        str(path.relative_to(verify.ROOT)): verify.sha(path)
        for path in verify.source_files()
    }


@pytest.mark.parametrize("key", KEYS)
def test_every_whole_report_field_mutation_is_rejected(key):
    data = verify.build_report()
    bad = copy.deepcopy(data)
    bad[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_and_extra_fields_are_rejected(mutation):
    data = verify.build_report()
    bad = copy.deepcopy(data)
    if mutation == "missing":
        del bad["whole_full_spatial_average_operator_bounds_and_coupled_turnaround"]
    else:
        bad["original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_full_counts_and_original_research_boundaries():
    data = verify.build_report()
    assert len(data["source_sha256"]) == 21
    assert (data["named_exact_check_count"], data["checked_scalar_entries"]) == (
        54,
        1331,
    )
    assert len(data["proof_checks"]) == 66
    assert len(data["controls"]) == 11 and data["controls"]["rejected_inputs"] == 396
    assert "NOT_HOMOGENEOUS_QUANTIZATION_OR_ORIGINAL_V_G_B_P8" in data["status"]
    assert (
        data["prior_sha256"][
            "S261_refutation_and_S265_integrability_boundary_unchanged"
        ]
        is True
    )
    assert (
        data["prior_sha256"][
            "new_model_is_classical_homogeneous_quantum_nonzero_hybrid"
        ]
        is True
    )
