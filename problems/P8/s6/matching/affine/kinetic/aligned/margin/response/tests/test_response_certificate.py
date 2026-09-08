"""Read-only response report with full local source and mutation controls."""
import copy
import json

import pytest
from p8_margin_response import verify

KEYS = ("schema", "claim", "date", "status", "prior_sha256", "source_sha256",
        "formulation", "written_proofs", "exact_residuals", "named_exact_check_count",
        "checked_scalar_entries", "proof_checks", "literal_model_state_and_source",
        "physical_source_variation", "regular_canonical_reduction",
        "normalized_two_component_system", "continuous_inverse_bound",
        "reconstruction_and_time_derivative_bounds", "exact_positive_physical_chart_lift",
        "approximate_metric_turning_minimum", "exact_derivative_majorants",
        "exact_operator_and_reconstruction_constants", "rounded_physical_scale_example",
        "actual_vector_scale_example", "quantum_feedback_boundary",
        "literature_scope", "controls", "verdict", "not_established", "verification_boundary")


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
        del changed["continuous_inverse_bound"]
    else:
        changed["full_quantum_feedback_is_controlled"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_counts_and_rejection_boundary():
    d = verify.build_report()
    assert d["named_exact_check_count"] == 22
    assert d["checked_scalar_entries"] == 23
    assert len(d["proof_checks"]) == 36
    assert d["controls"]["rejected_inputs"] == 34
    assert len(d["source_sha256"]) == 12
