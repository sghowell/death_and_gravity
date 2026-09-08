"""Read-only new-action report and complete local source manifest."""
import copy
import json

import pytest
from p8_aligned_margin import verify

KEYS = ("schema", "claim", "date", "status", "prior_sha256", "source_sha256",
        "formulation", "written_proofs", "exact_residuals", "named_exact_check_count",
        "checked_scalar_entries", "proof_checks", "literal_new_action",
        "background_and_full_volume_variation", "literal_constraints_and_regular_Hamiltonian",
        "two_chart_principal_derivation", "physical_characteristics",
        "chart_coverage_and_continuous_margins", "nonlinear_constraint_count_scope",
        "global_polynomial_bound_data", "center_data", "finite_potential_jet_diagnostic",
        "finite_potential_interpretation", "physical_scale_example",
        "unchanged_vector_first_variation_example", "controls", "verdict",
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
        del changed["two_chart_principal_derivation"]
    else:
        changed["full_quantum_cone_or_P8_closure"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_counts_and_rejection_boundary():
    d = verify.build_report()
    assert d["named_exact_check_count"] == 33
    assert d["checked_scalar_entries"] == 62
    assert len(d["proof_checks"]) == 15
    assert d["controls"]["rejected_inputs"] == 51
    assert len(d["source_sha256"]) == 13
