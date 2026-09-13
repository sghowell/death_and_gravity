"""Frozen entire homogeneous physical-background Gaussian vertex report."""

import copy
import json

import pytest
from p8_vacuum_affine_physical_background_vertices import verify

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
    "entire_current_parent_germs_coupled_reduction_and_joint_measure",
    "whole_physical_vertices_embedding_contacts_and_fixed_product_response",
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
        del bad["whole_physical_vertices_embedding_contacts_and_fixed_product_response"]
    else:
        bad["original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_exact_counts_and_unchanged_physical_boundary():
    data = verify.build_report()
    assert len(data["source_sha256"]) == 17
    assert (data["named_exact_check_count"], data["checked_scalar_entries"]) == (
        104,
        956,
    )
    assert len(data["proof_checks"]) == 39
    assert len(data["controls"]) == 9 and data["controls"]["rejected_inputs"] == 182
    assert (
        "EXACT_COMPLETE_HOMOGENEOUS_GAUGE_FIXED_GAUSSIAN_VERTICES_JOINT_MEASURE_AND_FIXED_PRODUCT_RESPONSE"
        in data["status"]
    )
    assert (
        "NOT_FULL_COVARIANT_VERTEX_MEASURE_CURVED_SUBTRACTION_PHYSICAL_MEAN_LOOP_BOUND_NONLINEAR_UV_REGGE_ORIGINAL_V_G_B_OR_P8"
        in data["status"]
    )
    assert (
        data["prior_sha256"]["no_interacting_state_or_QG3_retuning_substituted"] is True
    )
