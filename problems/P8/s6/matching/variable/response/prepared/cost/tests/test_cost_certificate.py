"""Full fresh ancestry, source hashes, exact serialization and report guards."""
import copy
import json

import pytest
import sympy as sp
from p8_preparation_cost import verify


@pytest.fixture(scope="module")
def rebuilt():
    return verify.build_report()


def test_frozen_report_matches_fresh_build(rebuilt):
    verify.validate_report(json.loads(verify.REPORT.read_text()), rebuilt)


def test_all_live_source_hashes_match_the_manifest(rebuilt):
    files = {str(path.relative_to(verify.ROOT)): path for path in verify.source_files()}
    assert set(files) == set(rebuilt["source_sha256"])
    assert all(verify.sha(path) == rebuilt["source_sha256"][name] for name, path in files.items())


def test_report_pins_the_actual_prepared_full_response_ancestry(rebuilt):
    assert verify.PREPARED_SHA in rebuilt["prior_sha256"].values()
    assert rebuilt["claim"] == "P8-S6.29.COST"
    assert "ORIGINAL_S6_P8_OPEN" in rebuilt["status"]


def test_full_identities_and_independent_comparisons(rebuilt):
    counts = [len(values) for values in rebuilt["exact_residuals"].values()]
    assert counts == [44, 60, 7, 8]
    assert sum(int(row["comparisons"]) for row in rebuilt["literal_rational_to_Arb_coefficient_bridges"]) == 256
    assert rebuilt["independent_Fraction_source_Hermite_replay"]["exact_scalar_function_or_polynomial_comparisons"] == 291
    assert rebuilt["independent_rational_Gramian_LDL"]["strict_positive_pivots"] == 4


def test_evidence_and_domain_controls_are_not_omitted(rebuilt):
    assert rebuilt["controls"]["rejected_inputs"] == 34
    assert rebuilt["controls"]["H0_2_minimum_at_boundary_budget_not_assumed"]
    assert rebuilt["controls"]["computed_band_moments_and_optimizer_not_claimed"]


def test_changed_report_does_not_replay(rebuilt):
    changed = copy.deepcopy(rebuilt)
    changed["continuous_K_energy_and_singular_value_bounds"]["uniform_Frobenius_gramian_lower"] = "1"
    with pytest.raises(ValueError):
        verify.validate_report(changed, rebuilt)


@pytest.mark.parametrize("invalid", [0.1, sp.Float("0.1"), sp.oo, sp.nan, sp.I, sp.Symbol("x")])
def test_exact_report_data_cannot_hide_inexact_or_free_expressions(invalid):
    with pytest.raises(TypeError):
        verify.serialize(invalid)


def test_exact_pi_in_the_Fourier_trace_is_not_a_decimal():
    assert verify.serialize(1/sp.pi) == "1/pi"
