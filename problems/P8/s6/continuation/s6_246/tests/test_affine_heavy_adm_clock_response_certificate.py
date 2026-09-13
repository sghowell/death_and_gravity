"""Immutable complete actual ADM/common-clock report and every-field mutation checks."""

import copy
import json

import pytest
from p8_vacuum_affine_heavy_adm_clock_response import verify

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
    "complete_actual_scalar_metric_and_clock",
    "complete_reference_graph_and_restored_projection",
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
        del bad["complete_reference_graph_and_restored_projection"]
    else:
        bad["full_quantum_inverse_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_exact_counts_and_complete_actual_response_boundary():
    data = verify.build_report()
    assert len(data["source_sha256"]) == 17
    assert (data["named_exact_check_count"], data["checked_scalar_entries"]) == (
        65,
        200,
    )
    assert len(data["proof_checks"]) == 50
    assert len(data["controls"]) == 9 and data["controls"]["rejected_inputs"] == 172
    assert (
        "COMPLETE_ACTUAL_GAUSSIAN_ADM_COMMON_CLOCK_REFERENCE_RESPONSE_WITH_WARD_COMPLETED_ORIGINAL_PROJECTION"
        in data["status"]
    )
    assert (
        "NOT_BARE_CUTOFF_LOCAL_SHAPE_IDENTIFICATION_QUANTUM_INVERSE_NONLINEAR_UV_REGGE_ORIGINAL_V_G_B_OR_P8"
        in data["status"]
    )
    assert (
        data["prior_sha256"][
            "all_frozen_original_and_separate_model_scientific_bytes_unchanged"
        ]
        is True
    )
