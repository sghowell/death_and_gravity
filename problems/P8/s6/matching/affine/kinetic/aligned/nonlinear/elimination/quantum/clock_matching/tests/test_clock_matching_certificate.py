"""Read-only fixed-prescription matching and exact source manifest."""
import copy
import json

import pytest
from p8_vector_clock_matching import verify

KEYS = ("schema", "claim", "date", "status", "source_sha256", "prior_sha256",
        "formulation", "written_proofs", "exact_residuals", "named_exact_check_count",
        "checked_scalar_entries", "proof_checks", "literal_component", "named_prescription",
        "counterterm_density", "first_variation_argument", "dimensional_scalar_curvature_tensor",
        "mass_pole_extension_and_evanescent_energy", "matched_local_first_variation_coefficients",
        "actual_clock_local_coefficients_at_mu_mass", "continuous_actual_clock_local_envelopes",
        "complex_dimension_continuation", "complex_dimension_coefficient_and_tail_envelopes",
        "uniform_integrable_regulator_limit", "matched_finite_part", "physical_scale_example",
        "physical_bound_boundary", "local_boundary_terms", "controls", "verdict",
        "not_established", "verification_boundary")


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_local_sources_and_report_keys_are_pinned():
    report = verify.build_report()
    assert report["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}
    assert set(KEYS) == set(report)


@pytest.mark.parametrize("key", KEYS)
def test_report_mutation_is_rejected(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_or_extra_report_fields_are_rejected(mutation):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    if mutation == "missing":
        del changed["uniform_integrable_regulator_limit"]
    else:
        changed["full_quantum_bounce_and_VGB_established"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_counts_and_rejection_boundary():
    d = verify.build_report()
    assert d["named_exact_check_count"] == 43
    assert d["checked_scalar_entries"] == 45
    assert len(d["proof_checks"]) == 24
    assert d["controls"]["rejected_inputs"] == 59
    assert d["controls"]["complex_dimension_not_conjugated_in_analytic_pair"] is True
    assert len(d["source_sha256"]) == 14
