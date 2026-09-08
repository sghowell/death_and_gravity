"""Read-only dimensional matching replay and exact source manifest."""
import copy
import json

import pytest
from p8_vector_dimensional import verify


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_local_sources_are_pinned():
    assert verify.build_report()["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}


@pytest.mark.parametrize("key", ("claim", "date", "status", "source_sha256", "prior_sha256", "exact_residuals",
                                  "named_exact_check_count", "checked_scalar_entries", "literal_component",
                                  "dimensional_canonical_normalization", "radial_Gamma_prescription",
                                  "actual_lapse_radial_poles_and_finite_parts", "physical_pressure_radial_poles_and_finite_parts",
                                  "radial_finite_part_boundary", "dimension_dependent_local_heat_coefficients_and_finite_action",
                                  "ordinary_counterterm_variation", "ordinary_covariantly_matched_finite_local_diagnostic",
                                  "ordinary_finite_matching_boundary", "actual_clock_mass_covariant_pole_reconstruction",
                                  "actual_clock_mass_pole_boundary", "finite_extension_boundary", "local_boundary_terms",
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
        del changed["finite_extension_boundary"]
    else:
        changed["actual_finite_clock_quantum_matching_established"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_counts_and_rejection_boundary():
    d = verify.build_report()
    assert d["named_exact_check_count"] == 71
    assert d["checked_scalar_entries"] == 71
    assert d["controls"]["rejected_inputs"] == 33
    assert d["controls"]["ordinary_finite_diagnostic_not_actual_clock_finite_matching"] is True
