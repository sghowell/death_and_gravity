"""Exact source manifest, frozen ancestry and semantic replay rejection."""
import copy
import json

import pytest
from p8_vector_subtraction import verify


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_local_sources_are_pinned():
    assert verify.build_report()["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}


@pytest.mark.parametrize("key", ("claim", "date", "status", "source_sha256", "prior_sha256", "exact_residuals",
                                  "named_exact_check_count", "checked_scalar_entries", "proof_checks", "literal_component",
                                  "physical_reference_weights", "derivative_order_prescription", "generic_reference_and_exact_tail",
                                  "actual_full_subtraction_coefficients", "continuous_reference_tail_bounds", "reference_tail_and_momentum_integral",
                                  "full_subtracted_exact_mode_integral", "physical_units_and_scheme_boundary", "physical_scale_example",
                                  "independent_ordinary_Proca_integrated_control", "incomplete_subtraction_control", "incomplete_subtraction_consequence",
                                  "controls", "verdict", "not_established", "verification_boundary"))
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
        del changed["physical_units_and_scheme_boundary"]
    else:
        changed["full_covariant_quantum_matching_established"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_counts_and_rejection_boundary():
    d = verify.build_report()
    assert d["named_exact_check_count"] == 11
    assert d["checked_scalar_entries"] == 11
    assert len(d["proof_checks"]) == 13
    assert d["controls"]["rejected_inputs"] == 50
    assert d["controls"]["full_subtracted_integral_not_automatic_covariant_matching_or_VGB"] is True
