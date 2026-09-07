"""Ordinary recursive replay and exact physical-input/report controls."""
import copy
import json
from fractions import Fraction

import pytest
import sympy as sp
from p8_exact_stationary.intervals import Interval
from p8_own_scattering import verify


def test_frozen_report_rebuilds_ordinary_ancestry():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_sources_and_separately_authored_audits_are_pinned():
    report = verify.build_report()
    sources = {str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}
    assert len(sources) == 23
    assert report["source_sha256"] == sources
    for name in ("arb_audit", "connection", "connection_audit", "interfaces", "jets", "potential"):
        assert "tests/test_own_scattering_"+name+".py" in sources


def test_exact_algebra_and_whole_domain_proof_counts():
    report = verify.build_report()
    assert [len(group) for group in report["exact_residuals"].values()] == [9, 45, 18]
    assert report["exact_residual_count"] == 72
    assert [len(group) for group in report["proof_checks"].values()] == [18, 16, 9]
    assert all(value is True for group in report["proof_checks"].values() for value in group.values())
    assert report["controls"]["rejected_inputs"] == 32
    assert "ORIGINAL_P8_OPEN" in report["status"]


@pytest.mark.parametrize("key", ("status", "prior_sha256", "source_sha256", "exact_residuals",
                                  "potential_domain", "transfer_domain", "probe_and_state_contract",
                                  "fraction_partial_jet_potential_proof", "finite_window_margins",
                                  "independent_Arb_total_jet_cover", "not_established"))
def test_changed_evidence_or_contract_is_rejected(key):
    report = verify.build_report()
    changed = copy.deepcopy(report)
    changed[key] = "changed"
    with pytest.raises(ValueError):
        verify.validate_report(changed, report)


def test_physical_window_and_potential_units_are_not_inner_units():
    data = verify.physical_contract(time_scale=2)
    assert data["physical_half_width"] == Fraction(1, 50)
    assert data["physical_duration"] == Fraction(1, 25)
    assert data["physical_potential_remainder_upper"] == 11
    assert data["homogeneous_mixing_lower"] == Fraction(1, 80)
    assert not data["fixed_temporal_band_or_EFT_cutoff"]
    assert not data["zero_data_matter_response"]


def test_complex_wave_output_is_exact_algebra_not_real_physical_input():
    assert verify.serialize(sp.I) == "I"
    assert verify.serialize(sp.Matrix([[1, sp.I]])) == [["1", "I"]]
    assert verify.serialize(Interval(Fraction(1, 3), Fraction(2, 3))) == {"lower": "1/3", "upper": "2/3"}
    with pytest.raises(TypeError):
        verify.physical_contract(delta=sp.I)
    for value in (1.0, sp.Float("0.1"), sp.oo, -sp.oo, sp.zoo, sp.nan):
        with pytest.raises(TypeError):
            verify.serialize(value)


def test_actual_mixing_not_a_vacuum_source_or_EFT_certificate():
    report = verify.build_report()
    checks = report["controls"]
    assert checks["conditional_remainder_premise_discharged_in_this_report"]
    assert checks["fictitious_free_exterior_not_physical_vacuum"]
    assert checks["fixed_physical_window_not_fixed_low_frequency_EFT_band"]
    assert checks["homogeneous_data_not_zero_data_source_loading"]
    assert checks["retarded_single_copy_variational_action_not_assumed"]
