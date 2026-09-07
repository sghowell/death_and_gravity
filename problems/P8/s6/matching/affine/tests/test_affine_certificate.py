"""Read-only replay, source manifest and nonzero/inexact negative controls."""
import copy
import json
from fractions import Fraction

import pytest
import sympy as sp
from p8_affine import verify


@pytest.fixture(scope="module")
def report():
    return verify.build_report()


def test_exact_frozen_replay(report):
    verify.validate_report(json.loads(verify.REPORT.read_text()), report)


def test_complete_current_source_manifest(report):
    assert len(report["source_sha256"]) == 20
    assert report["source_sha256"] == {str(path.relative_to(verify.ROOT)): verify.sha(path)
                                       for path in verify.source_files()}


def test_named_checks_not_inflated_by_matrix_entries(report):
    assert report["named_exact_check_count"] == 114
    assert report["checked_scalar_entries"] == 1617
    assert len(report["proof_checks"]["independent_domain_and_ODE"]) == 16
    assert len(report["proof_checks"]["premise_discharge"]) == 14


def test_controls_and_nonunit_physical_contract(report):
    assert report["controls"]["rejected_inputs"] == 46
    assert report["physical_nonunit_calibration"]["p_physical"] == "147/100"
    assert report["physical_nonunit_calibration"]["Delta_physical"] == "3"
    assert "ORIGINAL_P8_OPEN" in report["status"]
    assert report["full_connection_rank_and_inverse"]["UV_or_kinetic_health_claim"] is False


@pytest.mark.parametrize("key", ["status", "source_sha256", "proof_checks", "controls",
                                 "solution_lifting", "not_established"])
def test_tampered_report_rejected(report, key):
    changed = copy.deepcopy(report)
    changed[key] = "tampered"
    with pytest.raises(ValueError):
        verify.validate_report(changed, report)


@pytest.mark.parametrize("value", [1.0, sp.Float(1), sp.oo, -sp.oo, sp.zoo, sp.nan])
def test_inexact_or_nonfinite_certificate_entry_rejected(value):
    with pytest.raises((TypeError, ValueError)):
        verify.serialize({"nested": [value]})


def test_exact_serialization_and_matrix_metadata():
    assert verify.serialize([Fraction(2, 3), sp.Rational(4, 5), sp.Matrix([[0]])]) == ["2/3", "4/5", [["0"]]]
    assert verify.certify_residuals({"matrix": sp.zeros(2, 3)}) == {
        "matrix": {"shape": [2, 3], "all_entries_exactly_zero": True}}


@pytest.mark.parametrize("value", [sp.Integer(1), sp.Matrix([[0, 1]]), Fraction(1, 2)])
def test_nonzero_identity_rejected(value):
    with pytest.raises(ValueError):
        verify.certify_residuals({"not_zero": value})
