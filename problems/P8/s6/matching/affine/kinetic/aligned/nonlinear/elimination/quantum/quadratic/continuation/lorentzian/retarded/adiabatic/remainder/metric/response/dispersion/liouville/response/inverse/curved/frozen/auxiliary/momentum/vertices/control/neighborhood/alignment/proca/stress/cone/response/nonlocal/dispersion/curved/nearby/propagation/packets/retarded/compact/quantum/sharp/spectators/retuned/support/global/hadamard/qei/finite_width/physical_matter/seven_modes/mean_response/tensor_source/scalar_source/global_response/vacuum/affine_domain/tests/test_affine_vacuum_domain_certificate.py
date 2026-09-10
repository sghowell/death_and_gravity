"""Complete read-only ancestry, source hashes and every report-field mutation."""

import copy
import json

import pytest
from p8_affine_vacuum_domain import verify

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
    "new_analytic_reduced_family",
    "old_affine_extension_obstructions",
    "actual_principal_affine_lift",
    "regular_lower_dictionary",
    "general_gradient_quotient",
    "uniform_clock_tube_and_domain_bounds",
    "actual_vacuum_jets",
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
        del changed["regular_lower_dictionary"]
    else:
        changed["full_V_G_B_and_original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_exact_counts_and_original_scope():
    actual = verify.build_report()
    assert actual["named_exact_check_count"] == 70
    assert actual["checked_scalar_entries"] == 4930
    assert len(actual["proof_checks"]) == 51
    assert actual["controls"]["rejected_inputs"] == 54
    assert len(actual["source_sha256"]) == 20
    assert actual["claim"] == "P8-S6.109.REGULAR_ANALYTIC_CLASSICAL_AFFINE_VACUUM_LIFT"
