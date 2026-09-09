"""Full ancestry, every source byte and each report-field mutation."""

import copy
import json

import pytest
from p8_proca_mean_response import verify

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
    "physical_relative_Proca_stress",
    "actual_leading_mean_response",
    "physical_Ward_identity",
    "fixed_band_and_interval_bounds",
    "global_weighted_response_and_frame",
    "independent_center_stress_jets",
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
        del changed["actual_leading_mean_response"]
    else:
        changed["absolute_SEE_and_UV_completion_proved"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_exact_counts_and_original_scope():
    actual = verify.build_report()
    assert actual["named_exact_check_count"] == 65
    assert actual["checked_scalar_entries"] == 451
    assert len(actual["proof_checks"]) == 47
    assert actual["controls"]["rejected_inputs"] == 16
    assert len(actual["source_sha256"]) == 20
    assert actual["claim"] == "P8-S6.104.GLOBAL_RELATIVE_PROCA_MEAN_RESPONSE"
