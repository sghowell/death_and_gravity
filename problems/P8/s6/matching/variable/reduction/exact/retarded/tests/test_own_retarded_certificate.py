"""Ordinary frozen ancestry replay and exact response-contract controls."""
import copy
import json
from fractions import Fraction

import pytest
import sympy as sp
from p8_own_retarded import verify


def test_frozen_report_rebuilds_ordinary_ancestry():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_sources_and_three_scientific_audits_are_pinned():
    report = verify.build_report()
    sources = {str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}
    assert len(sources) == 15
    assert report["source_sha256"] == sources
    for name in ("test_own_retarded_action_audit.py", "test_own_retarded_bounds.py", "test_own_retarded_green.py"):
        assert "tests/"+name in sources


def test_exact_identities_proof_checks_and_controls():
    report = verify.build_report()
    assert [len(group) for group in report["exact_residuals"].values()] == [10, 14, 12]
    assert report["exact_residual_count"] == 36
    assert [len(group) for group in report["proof_checks"].values()] == [18, 22]
    assert all(value is True for group in report["proof_checks"].values() for value in group.values())
    assert report["controls"]["rejected_inputs"] == 34
    assert "ORIGINAL_P8_OPEN" in report["status"]


@pytest.mark.parametrize("key", ("status", "prior_sha256", "source_sha256", "exact_residuals",
                                  "domain", "input_state_contract", "finite_delta_coefficient_box",
                                  "Green_and_pulse_calibration", "literal_original_g_readout", "not_established"))
def test_changed_evidence_or_contract_is_rejected(key):
    report = verify.build_report()
    changed = copy.deepcopy(report)
    changed[key] = "changed"
    with pytest.raises(ValueError):
        verify.validate_report(changed, report)


def test_exact_nonunit_physical_readout_and_literal_variation_factor():
    d = verify.numeric_response(amplitude=Fraction(1, 100), parent_planck_squared=3, time_scale=2)
    assert d["Q0_lower"] == Fraction(63, 102400)
    assert d["Qx0_lower"] == Fraction(21, 8000)
    assert d["QT0_lower"] == sp.Rational(21, 640)
    assert d["physical_window_duration"] == sp.Rational(1, 50)
    assert d["germ_local_minimax_Q0_error_lower"] == Fraction(63, 204800)
    assert d["physical_original_g_equation_magnitude_lower"] == 3*128*Fraction(63, 102400)/(4*Fraction(1, 625)*Fraction(1251, 625))
    assert d["literal_g_action_derivative_magnitude_lower"] == d["physical_original_g_equation_magnitude_lower"]/2


def test_zero_history_comparator_and_retarded_state_are_part_of_contract():
    report = verify.build_report()
    assert "zero history OR" in report["input_state_contract"]
    assert "positive lower bounds concern the pulse subclass" in report["input_state_contract"]
    assert report["controls"]["zero_hidden_data_at_left_endpoint_explicit"]
    assert report["controls"]["flux_jump_is_one_not_R_derivative_one"]
    assert report["germ_local_minimax_Q0_error_per_amplitude"] == "63/2048"


def test_no_fixed_band_matter_source_or_variational_inverse_promotion():
    report = verify.build_report()
    c = report["controls"]
    assert c["shrinking_duration_is_not_a_fixed_low_frequency_class"]
    assert c["prescribed_metric_is_not_a_conserved_matter_source"]
    assert c["causal_equation_not_silently_a_single_copy_variational_action"]
    d = verify.numeric_response()
    assert not d["fixed_physical_low_frequency_band_claim"]
    assert not d["canonical_matter_source_or_undriven_physical_g_solution"]
