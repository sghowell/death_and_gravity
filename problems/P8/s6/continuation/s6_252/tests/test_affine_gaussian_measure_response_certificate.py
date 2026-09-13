"""Immutable Gaussian measure, CTP and continuum-boundary report."""

import copy
import json

import pytest
from p8_vacuum_affine_gaussian_measure_response import verify

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
    "whole_current_Gaussian_measure_and_fixed_state_kernels",
    "complete_CTP_vertices_and_continuum_Wick_extension_boundary",
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
        del bad["complete_CTP_vertices_and_continuum_Wick_extension_boundary"]
    else:
        bad["original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_exact_counts_and_unchanged_physical_boundary():
    data = verify.build_report()
    assert len(data["source_sha256"]) == 17
    assert (data["named_exact_check_count"], data["checked_scalar_entries"]) == (
        150,
        586,
    )
    assert len(data["proof_checks"]) == 26
    assert len(data["controls"]) == 9 and data["controls"]["rejected_inputs"] == 185
    assert (
        "EXACT_FULL_QUADRATIC_GAUSSIAN_MEASURE_AND_FINITE_CTP_WITH_WRITTEN_CONTINUUM_WICK_NOISE_AND_EXTENSION_FAMILY"
        in data["status"]
    )
    assert (
        "NOT_NONLINEAR_COVARIANT_MEASURE_SELECTED_CURVED_COUNTERFUNCTIONAL_PHYSICAL_LOOP_MEAN_BOUND_NONLINEAR_UV_REGGE_ORIGINAL_V_G_B_OR_P8"
        in data["status"]
    )
    assert (
        data["prior_sha256"]["no_full_nonlocal_or_interacting_state_substituted"]
        is True
    )
