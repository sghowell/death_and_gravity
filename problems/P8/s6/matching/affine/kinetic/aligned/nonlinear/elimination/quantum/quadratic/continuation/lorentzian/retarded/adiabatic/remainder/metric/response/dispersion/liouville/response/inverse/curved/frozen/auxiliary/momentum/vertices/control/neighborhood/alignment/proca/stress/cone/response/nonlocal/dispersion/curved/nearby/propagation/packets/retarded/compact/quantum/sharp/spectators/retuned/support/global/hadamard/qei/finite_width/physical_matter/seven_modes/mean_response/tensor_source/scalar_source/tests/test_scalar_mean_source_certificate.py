"""Complete ancestry, own source bytes and every report-field mutation."""

import copy
import json

import pytest
from p8_scalar_mean_source import verify

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
    "actual_scalar_quadratic_background_action",
    "actual_natural_scalar_phase",
    "actual_leading_scalar_sources",
    "actual_leading_scalar_mean_response",
    "actual_physical_matter_observable",
    "independent_scalar_center_jets",
    "whole_interval_state_and_response_bounds",
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
        str(path.relative_to(verify.ROOT)): verify.sha(path)
        for path in verify.source_files()
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
        del changed["actual_leading_scalar_sources"]
    else:
        changed["absolute_SEE_and_UV_completion_proved"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_exact_counts_and_original_scope():
    actual = verify.build_report()
    assert actual["named_exact_check_count"] == 33
    assert actual["checked_scalar_entries"] == 166
    assert len(actual["proof_checks"]) == 43
    assert actual["controls"]["rejected_inputs"] == 32
    assert len(actual["source_sha256"]) == 21
    assert actual["claim"] == "P8-S6.106.ACTUAL_COUPLED_SCALAR_MEAN_RESPONSE"
