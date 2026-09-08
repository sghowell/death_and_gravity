"""Read-only ancestry, exact source manifest and semantic report controls."""
import copy
import json

import pytest
from p8_vector_state import verify


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_local_sources_are_pinned():
    assert verify.build_report()["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}


@pytest.mark.parametrize("key", ("claim", "date", "status", "source_sha256", "prior_sha256", "exact_residuals",
                                  "named_exact_check_count", "checked_scalar_entries", "proof_checks", "literal_component",
                                  "actual_clock_mass_direction", "physical_canonical_energy_and_pressure", "local_pole_control",
                                  "local_pole_boundary", "domain_and_state_preparation", "exact_WKB_reference_coefficients",
                                  "continuous_WKB_envelopes", "reference_and_residual_bounds", "exact_evolution_comparison",
                                  "all_momentum_integral_and_units", "physical_scale_example", "controls", "verdict",
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
        del changed["domain_and_state_preparation"]
    else:
        changed["full_renormalized_energy_established"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_counts_and_rejection_boundary():
    d = verify.build_report()
    assert d["named_exact_check_count"] == 40
    assert d["checked_scalar_entries"] == 45
    assert len(d["proof_checks"]) == 17
    assert d["controls"]["rejected_inputs"] == 40
    assert d["controls"]["evolution_difference_not_full_renormalized_energy_or_VGB"] is True
