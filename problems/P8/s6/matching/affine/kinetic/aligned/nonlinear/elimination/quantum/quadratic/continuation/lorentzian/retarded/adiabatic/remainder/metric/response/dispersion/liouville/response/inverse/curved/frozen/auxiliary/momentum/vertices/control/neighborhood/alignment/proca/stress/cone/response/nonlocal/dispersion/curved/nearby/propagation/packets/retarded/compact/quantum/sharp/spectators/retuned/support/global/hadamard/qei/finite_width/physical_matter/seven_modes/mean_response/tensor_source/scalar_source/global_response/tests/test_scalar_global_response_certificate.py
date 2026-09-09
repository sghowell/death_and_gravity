"""Complete ancestry, source manifest and every report-field mutation."""

import copy
import json

import pytest
from p8_scalar_global_response import verify

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
    "actual_weighted_scalar_phase",
    "compensated_actual_mean_kernels",
    "whole_half_line_forcing_envelopes",
    "global_scalar_mean_and_observable_bounds",
    "independent_variable_and_asymptotic_controls",
    "joint_seven_mode_relative_response",
    "explicit_amplitude_calibrations",
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
        del changed["compensated_actual_mean_kernels"]
    else:
        changed["absolute_SEE_and_original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_exact_counts_and_original_scope():
    actual = verify.build_report()
    assert actual["named_exact_check_count"] == 47
    assert actual["checked_scalar_entries"] == 203
    assert len(actual["proof_checks"]) == 48
    assert actual["controls"]["rejected_inputs"] == 51
    assert len(actual["source_sha256"]) == 20
    assert actual["claim"] == "P8-S6.107.GLOBAL_SEVEN_MODE_RELATIVE_MEAN_CONTROL"
