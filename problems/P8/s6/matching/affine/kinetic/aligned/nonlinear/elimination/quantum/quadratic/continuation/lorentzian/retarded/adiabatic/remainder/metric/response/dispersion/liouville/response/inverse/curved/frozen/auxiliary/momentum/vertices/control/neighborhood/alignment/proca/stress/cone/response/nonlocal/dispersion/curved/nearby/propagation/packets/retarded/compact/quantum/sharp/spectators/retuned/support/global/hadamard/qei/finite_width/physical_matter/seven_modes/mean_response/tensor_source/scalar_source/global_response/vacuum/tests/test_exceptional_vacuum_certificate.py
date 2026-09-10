"""Complete native ancestry, source manifest and every report-field mutation."""

import copy
import json

import pytest
from p8_exceptional_vacuum import verify

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
    "smooth_classical_family",
    "naive_transition_control",
    "analytic_finite_jet_family",
    "actual_vacuum_jets",
    "clock_tube_four_jet_bounds",
    "uniform_co_scaled_family",
    "canonical_decoupling_limit",
    "actual_labelled_scalar_contact",
    "healthy_heavy_tree_matching",
    "exact_calibrations",
    "controls",
    "verdict",
    "not_established",
    "verification_boundary",
]


def test_readonly_complete_report_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_own_source_manifest_and_complete_report_fields():
    actual = verify.build_report()
    assert actual["source_sha256"] == {
        str(p.relative_to(verify.ROOT)): verify.sha(p) for p in verify.source_files()
    }
    assert set(actual) == set(KEYS)


@pytest.mark.parametrize("key", KEYS)
def test_each_report_field_mutation_rejected(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_or_extra_report_field_rejected(mutation):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    if mutation == "missing":
        del changed["healthy_heavy_tree_matching"]
    else:
        changed["full_V_G_B_and_original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_exact_counts_and_original_scope():
    actual = verify.build_report()
    assert actual["named_exact_check_count"] == 103
    assert actual["checked_scalar_entries"] == 103
    assert len(actual["proof_checks"]) == 53
    assert actual["controls"]["rejected_inputs"] == 61
    assert len(actual["source_sha256"]) == 21
    assert actual["claim"] == "P8-S6.108.EXCEPTIONAL_VACUUM_TREE_MATCHING"
