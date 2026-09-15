"""Immutable complete-source and twenty-field finite comparison contract."""

import copy
import json

import pytest
from p8_vacuum_affine_hybrid_core_regulator_comparison import verify

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
    "whole_source_full_quadratic_and_complete_nonlinear_remainder",
    "whole_exact_operators_and_coupled_core_regulator_comparisons",
    "observable_and_scope",
    "primitive_and_matching_frontier",
    "controls",
    "verdict",
    "not_established",
    "verification_boundary",
]


def test_complete_read_only_report_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_exact_complete_source_and_field_manifest():
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
def test_missing_and_extra_report_fields_are_rejected(mutation):
    data = verify.build_report()
    bad = copy.deepcopy(data)
    if mutation == "missing":
        del bad["whole_exact_operators_and_coupled_core_regulator_comparisons"]
    else:
        bad["original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_exact_counts_and_retained_original_research_boundaries():
    data = verify.build_report()
    assert len(data["source_sha256"]) == 19
    assert (data["named_exact_check_count"], data["checked_scalar_entries"]) == (
        61,
        289,
    )
    assert len(data["proof_checks"]) == 57
    assert len(data["controls"]) == 13 and data["controls"]["rejected_inputs"] == 496
    assert "NOT_EXACT_SUPPORT_REGULATOR_REMOVAL_OR_ORIGINAL_V_G_B_P8" in data["status"]
    assert (
        data["prior_sha256"][
            "same_four_finite_hybrid_solutions_not_replacement_dynamics"
        ]
        is True
    )
    assert (
        data["prior_sha256"][
            "S261_refutation_and_S265_integrability_boundary_unchanged"
        ]
        is True
    )
