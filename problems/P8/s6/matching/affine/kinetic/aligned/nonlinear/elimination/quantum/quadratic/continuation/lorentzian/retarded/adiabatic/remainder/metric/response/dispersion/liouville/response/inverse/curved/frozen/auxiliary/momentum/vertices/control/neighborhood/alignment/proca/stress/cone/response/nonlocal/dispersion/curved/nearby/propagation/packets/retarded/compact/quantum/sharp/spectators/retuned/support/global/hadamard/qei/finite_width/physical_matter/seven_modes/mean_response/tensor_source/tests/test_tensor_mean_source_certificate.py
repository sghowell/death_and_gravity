"""Full ancestry and own-source/report mutation replay."""

import copy
import json

import pytest
from p8_tensor_mean_source import verify

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
    "actual_tensor_action_and_lapse_source",
    "positive_actual_tensor_state_difference",
    "actual_leading_tensor_mean_response",
    "actual_clock_and_joint_source_Ward_identity",
    "nonlinear_center_physical_density_chart",
    "global_response_and_clock_transfer_bounds",
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
        del changed["actual_tensor_action_and_lapse_source"]
    else:
        changed["absolute_SEE_and_UV_completion_proved"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_exact_counts_and_original_scope():
    actual = verify.build_report()
    assert actual["named_exact_check_count"] == 53
    assert actual["checked_scalar_entries"] == 89
    assert len(actual["proof_checks"]) == 37
    assert actual["controls"]["rejected_inputs"] == 16
    assert len(actual["source_sha256"]) == 20
    assert actual["claim"] == "P8-S6.105.ACTUAL_TENSOR_CLOCK_MEAN_RESPONSE"
