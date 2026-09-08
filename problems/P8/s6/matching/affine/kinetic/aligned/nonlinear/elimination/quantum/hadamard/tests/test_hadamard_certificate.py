"""Read-only controlled Hadamard-state report and complete local source manifest."""
import copy
import json

import pytest
from p8_vector_hadamard import verify

KEYS = ("schema", "claim", "date", "status", "prior_sha256", "source_sha256",
        "formulation", "written_proofs", "exact_residuals", "named_exact_check_count",
        "checked_scalar_entries", "proof_checks", "literal_component", "all_order_recurrence",
        "locally_finite_state_definition", "sixth_order_and_full_correction_envelopes",
        "first_higher_cutoff_controls", "explicit_Cauchy_data_examples", "positive_normalized_preparation",
        "Proca_reference_construction", "imported_theorem_and_hypotheses", "all_order_Hadamard_comparison",
        "physical_vector_reconstruction", "quantitative_observable_transfer",
        "finite_state_change_regulator_limit", "physical_scale_example", "physical_units_and_scope",
        "controls", "verdict", "not_established", "verification_boundary")


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
        del changed["all_order_Hadamard_comparison"]
    else:
        changed["old_state_relabelled_or_full_quantum_solution_established"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_counts_and_rejection_boundary():
    d = verify.build_report()
    assert d["named_exact_check_count"] == 36
    assert d["checked_scalar_entries"] == 61
    assert len(d["proof_checks"]) == 26
    assert d["controls"]["rejected_inputs"] == 113
    assert d["controls"]["three_positive_constrained_Proca_polarizations_not_four_scalars"] is True
    assert len(d["source_sha256"]) == 15
