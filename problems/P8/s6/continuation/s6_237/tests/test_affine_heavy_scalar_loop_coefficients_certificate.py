"""Immutable separate first-loop low coefficient matching report and every-field mutation rejection."""

import copy
import json

import pytest
from p8_vacuum_affine_heavy_scalar_loop_coefficients import verify

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
    "complete_four_order_loop_and_logarithmic_moment_reduction",
    "full_joint_bidisk_and_actual_coefficient_matching_bounds",
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
        del bad["full_joint_bidisk_and_actual_coefficient_matching_bounds"]
    else:
        bad["physical_cutoff"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_exact_counts_and_conditional_boundary():
    data = verify.build_report()
    assert len(data["source_sha256"]) == 17
    assert (
        data["named_exact_check_count"] == 60 and data["checked_scalar_entries"] == 60
    )
    assert len(data["proof_checks"]) == 36
    assert len(data["controls"]) == 9 and data["controls"]["rejected_inputs"] == 273
    assert "SEPARATE_COMPLETE_FIRST_LOOP_LOW_COEFFICIENT_MATCHING" in data["status"]
    assert (
        "NOT_OMITTED_LOOPS_EXACT_PHYSICAL_COEFFICIENTS_FULL_UV_ORIGINAL_AFFINE_PARENT"
        in data["status"]
    )
    assert (
        data["prior_sha256"][
            "all_frozen_original_and_separate_model_scientific_bytes_unchanged"
        ]
        is True
    )
