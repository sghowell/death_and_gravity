"""Read-only higher-regularity certificate and full source-manifest checks."""
import copy
import json

import pytest
from p8_vector_regularity import verify

KEYS = ("schema", "claim", "date", "status", "prior_sha256", "source_sha256",
        "formulation", "written_proofs", "exact_residuals", "named_exact_check_count",
        "checked_scalar_entries", "proof_checks", "fixed_input_boundary",
        "frequency_residual_envelopes", "reference_residual_definition",
        "initial_preparation_envelopes", "frequency_band_boundary",
        "oscillatory_evolution_boundary", "actual_readout_row_rule",
        "projected_reference_boundary", "readout_derivative_envelopes",
        "lower_reference_failure_control", "radial_integral_envelope",
        "radial_bound_formula", "local_derivative_coefficients_and_envelopes",
        "physical_scale_example", "continuous_regularity_conclusion", "controls",
        "verdict", "not_established", "verification_boundary")


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
        del changed["projected_reference_boundary"]
    else:
        changed["fixed_background_bounds_are_full_quantum_feedback"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_counts_and_rejection_boundary():
    data = verify.build_report()
    assert data["named_exact_check_count"] == 86
    assert data["checked_scalar_entries"] == 381
    assert len(data["proof_checks"]) == 98
    assert data["controls"]["rejected_inputs"] == 120
    assert len(data["source_sha256"]) == 14
