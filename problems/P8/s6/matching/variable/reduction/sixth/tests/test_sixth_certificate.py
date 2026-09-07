"""Ordinary report replay and promotion-boundary controls."""
import copy
import json

import pytest
import sympy as sp
from p8_sixth_reduction import verify


def test_frozen_report_rebuilds_ordinary_ancestry():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_every_source_and_both_scientific_audits_are_pinned():
    actual = verify.build_report()
    sources = {str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}
    assert actual["source_sha256"] == sources
    assert len(sources) == 14
    assert "tests/test_sixth_independent_audit.py" in sources
    assert "tests/test_sixth_center_independent_audit.py" in sources


def test_exact_counts_domain_and_derivative_guards():
    actual = verify.build_report()
    assert [len(group) for group in actual["exact_residuals"].values()] == [30, 40]
    assert actual["exact_residual_count"] == 70
    assert actual["controls"]["rejected_inputs"] == 26
    assert "ORIGINAL_P8_OPEN" in actual["status"]


@pytest.mark.parametrize("key", ("status", "prior_sha256", "source_sha256", "exact_residuals",
                                  "actual_center_sixth_E_physical", "actual_center_combined_E_dimensionless",
                                  "individual_strict_floor", "mandatory_combined_low_symbol_polynomial",
                                  "not_established"))
def test_changed_evidence_is_not_accepted(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "changed"
    with pytest.raises(ValueError):
        verify.validate_report(changed, actual)


def test_literal_physical_numeric_calibration_and_nonclaims():
    values = verify.numeric_center(4, 3, 2)
    assert values["sixth_E"][2] == -312
    assert values["sixth_E"][4] == -264
    assert values["sixth_E"][6] == -6
    assert values["fourth_E2"] == 30
    assert values["combined_E2"] == -282
    assert not values["probe_is_actual_FLRW_background"]
    assert not values["covariant_Xi_extracted"]
    assert verify.numeric_center(sp.Rational(201, 100))["sixth_E"][2] < -sp.Rational(793, 400)


def test_individual_floor_is_not_a_combined_remainder_bound():
    actual = verify.build_report()
    assert not actual["individual_strict_floor"]["is_a_combined_or_full_remainder_bound"]
    assert actual["controls"]["fourth_order_cancellation_retained"]
    assert actual["mandatory_combined_low_symbol_polynomial"][-1] == "0"
    assert actual["continuous_individual_low_symbol_polynomial"][-1] == "-2"
