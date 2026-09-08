"""Read-only rebuild, exact local manifest and semantic mutation rejection."""
import copy
import json

import pytest
from p8_vector_curvature import verify


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_local_sources_are_pinned():
    assert verify.build_report()["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}


@pytest.mark.parametrize("key", ("claim", "date", "status", "source_sha256", "prior_sha256", "exact_residuals",
                                  "named_exact_check_count", "checked_scalar_entries", "proof_checks", "literal_component",
                                  "determinant_prescription", "local_heat_coefficients", "regulator_and_pole_normalization",
                                  "actual_physical_metric_invariants_and_coefficients", "continuous_binomial_envelopes",
                                  "sharp_scalar_curvature_and_units", "boundary_terms", "physical_scale_example", "pole_ratio_boundary",
                                  "independent_S4_spectral_control", "S4_mode_and_domain_boundary", "controls", "verdict",
                                  "not_established", "verification_boundary"))
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
        del changed["boundary_terms"]
    else:
        changed["full_quantum_bounce_established"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_counts_and_rejection_boundary():
    d = verify.build_report()
    assert d["named_exact_check_count"] == 44
    assert d["checked_scalar_entries"] == 81
    assert len(d["proof_checks"]) == 21
    assert d["controls"]["rejected_inputs"] == 22
    assert d["controls"]["local_pole_and_S4_error_not_full_FLRW_in_in_error"] is True
