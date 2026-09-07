"""Ordinary recursive report replay, physical units and evidence controls."""
import copy
import json
from fractions import Fraction

import pytest
from p8_own_forced import verify


def test_frozen_report_rebuilds_ordinary_ancestry():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_sources_and_independent_scientific_audits_are_pinned():
    report = verify.build_report()
    files = {str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}
    assert len(files) == 17
    assert report["source_sha256"] == files
    for name in ("audit", "interfaces", "loading", "phase"):
        assert "tests/test_own_forced_"+name+".py" in files


def test_exact_proof_and_arithmetic_control_counts():
    report = verify.build_report()
    assert report["exact_residual_count"] == 64
    assert [len(group) for group in report["exact_residuals"].values()] == [8, 25, 13, 18]
    assert [len(group) for group in report["proof_checks"].values()] == [22, 17, 11, 14]
    assert all(value is True for group in report["proof_checks"].values() for value in group.values())
    assert report["controls"]["rejected_inputs"] == 45
    assert "ORIGINAL_P8_OPEN" in report["status"]


@pytest.mark.parametrize("key", ("status", "prior_sha256", "source_sha256", "exact_residuals",
                                  "domain", "input_state_contract", "finite_parameter_loading",
                                  "phase_comparison_and_scalar_readout", "explicit_fixed_physical_pulse",
                                  "original_Q_limsup_minus_liminf_lower_per_eta", "not_established"))
def test_altered_evidence_or_source_state_contract_is_rejected(key):
    report = verify.build_report()
    changed = copy.deepcopy(report)
    changed[key] = "changed"
    with pytest.raises(ValueError):
        verify.validate_report(changed, report)


def test_exact_nonunit_physical_calibration():
    data = verify.physical_contract(amplitude=Fraction(1, 100), time_scale=2)
    assert data["initial_physical_time"] == Fraction(-1, 50)
    assert data["loading_physical_time"] == Fraction(-1, 100)
    assert data["observation_physical_time"] == Fraction(1, 100)
    assert data["total_physical_duration"] == Fraction(3, 100)
    assert data["physical_plateau_width"] == Fraction(1, 640)
    assert data["loaded_Y_lower"] == Fraction(9, 128000)
    assert data["loaded_Y_T_lower"] == Fraction(1, 48)
    assert data["limsup_minus_liminf_original_Q_lower"] == Fraction(1, 20000)
    assert data["stronger_original_Q_gap_lower"] == Fraction(237, 4256000)


def test_limsup_gap_is_not_a_positive_response_floor_at_every_delta():
    data = verify.physical_contract()
    assert data["gap_not_a_pointwise_positive_Q_floor_for_every_delta"]
    assert data["same_q_and_tau_for_every_delta"]
    assert not data["physical_low_frequency_band_or_cutoff_claimed"]
    assert not data["canonical_matter_or_undriven_g_solution_claimed"]


def test_central_error_and_effective_description_boundaries_are_explicit():
    controls = verify.build_report()["controls"]
    assert controls["both_actual_error_neighborhoods_retained"]
    assert controls["strict_limiting_loading_uses_stronger_uniform_margin"]
    assert controls["actual_remainder_subsequential_limits_not_assumed"]
    assert controls["delta_dependent_phase_retaining_effective_descriptions_not_excluded"]
    assert controls["retarded_inverse_not_silently_a_single_copy_variational_action"]
