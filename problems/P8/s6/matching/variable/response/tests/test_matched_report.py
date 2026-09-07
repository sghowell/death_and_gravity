"""Read-only evidence inventory, parent replay and corruption controls."""

import copy
import json

import pytest
import sympy as sp
from p8_variable_response import verify


def test_frozen_matched_report_replays_exactly():
    actual = verify.build_report()
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual)


def test_manifest_is_complete_and_owns_no_ancestor():
    report = verify.build_report()
    assert report["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()
    }
    assert "tests/test_response_independent_audit.py" in report["source_sha256"]
    assert all(not name.startswith("../") for name in report["source_sha256"])
    assert all("certificates/" not in name for name in report["source_sha256"])


def test_parent_report_is_replayed_not_just_a_trusted_number():
    pins = verify.prior_checks()
    assert pins["S6_20_actual_variable_background_and_complete_canonical_TT"] == verify.VARIABLE_SHA
    assert pins["adopted_S6_contract_unchanged"] == verify.CONTRACT_SHA


@pytest.mark.parametrize("field", ["source_sha256", "norm_calibration", "exact_connection",
                                  "physical_source_and_data", "not_established"])
def test_any_claim_or_evidence_mutation_is_rejected(field):
    actual = verify.build_report()
    tampered = copy.deepcopy(actual)
    tampered[field] = {}
    with pytest.raises(ValueError, match="read-only exact replay"):
        verify.validate_report(tampered, actual)


def test_report_does_not_promote_source_independent_or_low_band_verdict():
    report = verify.build_report()
    assert "OPEN" in report["status"]
    assert report["checked_controls"]["prepared_light_error_retains_source_L1"] is True
    assert report["checked_controls"]["full_mass_remainder_joint_analyticity"] is False
    assert report["norm_calibration"]["upper_box_error_not_small"] == "40"
    assert any("not a low temporal band" in entry for entry in report["not_established"])


@pytest.mark.parametrize("value", [0.1, sp.Float("0.1"), sp.oo, sp.nan])
def test_no_rounded_or_nonfinite_certificate_number(value):
    with pytest.raises(TypeError):
        verify.serialize(value)
