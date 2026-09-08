"""Read-only selected-state profile report and source-manifest controls."""
import copy
import json

import pytest
from p8_clock_tadpole import verify

KEYS = ("schema", "claim", "date", "status", "prior_sha256", "source_sha256",
        "formulation", "written_proofs", "exact_residuals", "named_exact_check_count",
        "checked_scalar_entries", "proof_checks", "literal_new_action",
        "fixed_profile_definition", "actual_mode_profile_integrands",
        "unchanged_initial_data_examples", "literal_clock_tube_and_global_profile",
        "actual_point_chart_quadratic_terms", "window_derivative_envelopes",
        "finite_profile_error_budget", "background_stationarity_and_completeness_scope",
        "global_smoothness_and_local_budget_boundary", "vacuum_preservation_boundary",
        "retained_vector_principal_data", "prior_approximate_response_cone_boundary",
        "exact_vector_block_examples", "remaining_variation_boundary", "controls",
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
        del changed["remaining_variation_boundary"]
    else:
        changed["all_quantum_variations_cancelled"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_counts_and_rejection_boundary():
    data = verify.build_report()
    assert data["named_exact_check_count"] == 30
    assert data["checked_scalar_entries"] == 30
    assert len(data["proof_checks"]) == 50
    assert data["controls"]["rejected_inputs"] == 110
    assert len(data["source_sha256"]) == 14
