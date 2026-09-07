"""Read-only ordinary ancestry/report replay and claim-boundary controls."""
import copy
import json

import pytest
import sympy as sp
from p8_exact_stationary import verify


def test_frozen_report_rebuilds_ordinary_ancestry():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_every_source_and_all_three_scientific_audits_are_pinned():
    report = verify.build_report()
    sources = {str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}
    assert report["source_sha256"] == sources
    assert len(sources) == 16
    for name in ("test_exact_stationary.py", "test_exact_intervals.py", "test_exact_independent_audit.py"):
        assert "tests/"+name in sources


def test_exact_identities_interval_proofs_and_rejection_counts():
    report = verify.build_report()
    assert [len(group) for group in report["exact_residuals"].values()] == [37, 14]
    assert report["exact_residual_count"] == 51
    assert len(report["fraction_whole_box"]["checks"]) == 31
    assert len(report["independent_arb_whole_box"]["strict_rational_enclosures"]) == 8
    assert report["controls"]["rejected_inputs"] == 35
    assert "ORIGINAL_P8_OPEN" in report["status"]


@pytest.mark.parametrize("key", ("status", "prior_sha256", "source_sha256", "exact_residuals",
                                  "actual_domain", "fraction_whole_box", "independent_arb_whole_box",
                                  "center_jets", "inner_equation_coefficients", "not_established"))
def test_changed_evidence_or_scope_is_rejected(key):
    report = verify.build_report()
    changed = copy.deepcopy(report)
    changed[key] = "changed"
    with pytest.raises(ValueError):
        verify.validate_report(changed, report)


def test_literal_physical_center_units_and_domain():
    data = verify.numeric_center(parent_planck_squared=3, time_scale=2)
    assert data["c"] == sp.Rational(201, 100)
    assert data["N"] == sp.Rational(40401, 20000)
    assert data["f00"] == sp.Rational(40401, 20000)**2
    assert data["b_TT"] == -sp.Rational(201*401, 20000)
    assert data["K_f"] == sp.Rational(480000, 40401)
    assert data["own_f_mass_squared"] == 402
    assert not data["full_parent_solution"]
    assert not data["physical_cutoff_claimed"]


def test_original_action_and_off_shell_probe_are_retained():
    report = verify.build_report()
    assert report["prior_sha256"]["does_not_use_S6_31_deformed_interactions"]
    assert report["prior_sha256"]["same_off_shell_unit_volume_probe_not_actual_FLRW"]
    assert not report["fraction_whole_box"]["tensor_Green_inverse_selected"]
    assert report["controls"]["lapse_uses_total_implicit_derivative"]
    assert report["controls"]["literal_c_equals_2_excluded"]


def test_proof_booleans_are_not_zero_symbolic_residuals():
    report = verify.build_report()
    assert all(value is True for value in report["fraction_whole_box"]["checks"].values())
    assert all(value == "0" for group in report["exact_residuals"].values() for value in group.values())
