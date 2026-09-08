"""Read-only dimensional continuation report and mutation controls."""
import copy
import json

import pytest
from p8_vector_quadratic_matching import verify

KEYS = ("schema", "claim", "date", "status", "prior_sha256", "source_sha256",
        "formulation", "written_proofs", "exact_residuals", "named_exact_check_count",
        "checked_scalar_entries", "proof_checks", "literal_continuation",
        "continued_auxiliary_curvatures", "quadratic_pole_dimensional_jets",
        "actual_compact_dimensional_jets", "normalized_lapse_operator_dimensional_jets",
        "bare_counterterm_finite_evanescent_operators", "measure_variation_controls",
        "Gauss_Bonnet_boundary_control", "counterterm_convention", "scope", "controls",
        "verdict", "not_established", "verification_boundary")


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_sources_and_keys_pinned():
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
def test_missing_or_extra_field_rejected(mutation):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    if mutation == "missing":
        del changed["measure_variation_controls"]
    else:
        changed["full_finite_response_complete"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_and_control_counts():
    report = verify.build_report()
    assert report["named_exact_check_count"] == report["checked_scalar_entries"] == 41
    assert len(report["proof_checks"]) == 11
    assert report["controls"]["rejected_inputs"] == 189
    assert len(report["source_sha256"]) == 16
