"""Read-only ancestry, exact source manifest and semantic mutation controls."""
import copy
import json

import pytest
from p8_vector_bubble import verify


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_local_sources_are_pinned():
    assert verify.build_report()["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}


@pytest.mark.parametrize("key", ("claim", "date", "status", "source_sha256", "prior_sha256", "exact_residuals",
                                  "named_exact_check_count", "checked_scalar_entries", "proof_checks",
                                  "literal_component", "regulator_and_normalization", "general_pole_polynomial",
                                  "actual_mass_direction_and_pole", "continuous_pole_bound", "radial_poles_and_finite_parts",
                                  "finite_logarithmic_subtraction", "finite_remainder_constants",
                                  "continuous_finite_domain_and_bound", "local_subtraction_boundary", "physical_scale_example",
                                  "source_functional_boundary", "controls", "verdict", "not_established", "verification_boundary"))
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
        del changed["local_subtraction_boundary"]
    else:
        changed["full_quantum_completion_established"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_control_count_and_limited_admissibility_boundary():
    d = verify.controls()
    assert d["rejected_inputs"] == 49
    assert d["subtracted_loop_not_unknown_higher_local_matching_operators"] is True
    assert d["no_stationary_to_curved_in_in_transfer_asserted"] is True


def test_scalar_and_matrix_check_counts_remain_distinct():
    d = verify.build_report()
    assert d["named_exact_check_count"] == 27
    assert d["checked_scalar_entries"] == 42
