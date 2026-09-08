"""Read-only finite local mass-response report and mutation controls."""
import copy
import json

import pytest
from p8_vector_mass_adiabatic import verify

KEYS = ["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","moving_mode_source_variations","actual_radial_Laurent_coefficients","actual_finite_local_operators","actual_compact_finite_local_actions","omitted_counterterm_adjoint_controls","continuous_operator_coefficient_envelopes","scale_example","normalization_and_prescription","remaining_integral","controls","verdict","not_established","verification_boundary"]


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_sources_and_report_keys_pinned():
    data = verify.build_report()
    assert data["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}
    assert set(KEYS) == set(data)


@pytest.mark.parametrize("key", KEYS)
def test_report_mutation_rejected(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_or_extra_report_field_rejected(mutation):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    if mutation == "missing":
        del changed["remaining_integral"]
    else:
        changed["nonlocal_response_bounded"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_and_control_counts():
    data = verify.build_report()
    assert data["named_exact_check_count"] == 43
    assert data["checked_scalar_entries"] == 43
    assert len(data["proof_checks"]) == 14
    assert data["controls"]["rejected_inputs"] == 92
    assert len(data["source_sha256"]) == 16
