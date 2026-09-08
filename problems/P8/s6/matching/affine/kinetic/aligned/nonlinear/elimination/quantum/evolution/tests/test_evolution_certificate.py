"""Read-only C1 vector and clock-source report integrity."""
import copy
import json

import pytest
from p8_vector_evolution import verify

KEYS = ("schema", "claim", "date", "status", "prior_sha256", "source_sha256",
        "formulation", "written_proofs", "exact_residuals", "named_exact_check_count",
        "checked_scalar_entries", "proof_checks", "literal_component", "polynomial_majorant_method",
        "actual_residual_and_tail_derivative_envelopes", "phase_separated_evolution",
        "oscillatory_mixing_bound", "physical_readout_derivative", "differentiated_finite_integral_bound",
        "ordinary_conservation_controls", "direct_homogeneous_clock_variation",
        "clock_source_regulator_limit", "actual_local_derivative_and_clock_source_coefficients",
        "continuous_local_derivative_and_clock_source_envelopes", "physical_scale_example",
        "physical_units_and_scope", "controls", "verdict", "not_established", "verification_boundary")


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
        del changed["clock_source_regulator_limit"]
    else:
        changed["all_order_state_and_full_quantum_solution"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_counts_and_rejection_boundary():
    d = verify.build_report()
    assert d["named_exact_check_count"] == 54
    assert d["checked_scalar_entries"] == 57
    assert len(d["proof_checks"]) == 28
    assert d["controls"]["rejected_inputs"] == 70
    assert d["controls"]["both_oscillatory_endpoints_and_feedback_retained"] is True
    assert len(d["source_sha256"]) == 13
