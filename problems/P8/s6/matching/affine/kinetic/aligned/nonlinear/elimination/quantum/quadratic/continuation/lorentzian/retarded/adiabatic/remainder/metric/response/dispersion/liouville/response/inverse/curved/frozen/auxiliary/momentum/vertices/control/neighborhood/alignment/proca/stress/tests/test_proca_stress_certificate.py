"""Read-only new-profile certificate, pinned sources and complete field mutation."""
import copy
import json

import pytest
from p8_proca_stress import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","actual_new_energy_rows_and_continuous_envelopes","actual_differentiated_subtraction_and_integrable_reference_tails","actual_lower_reference_failure_control","actual_new_local_energy_derivative_envelopes","continuous_new_stress_and_energy_change_bounds","actual_fixed_continuum_integrand_definition","explicit_new_fixed_profile_and_old_profile_difference","global_mixed_profile_and_local_contact_bounds","controls","verdict","not_established","verification_boundary"]


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()),verify.build_report())


def test_all_source_files_and_report_keys_pinned():
    actual=verify.build_report()
    assert actual["source_sha256"]=={str(path.relative_to(verify.ROOT)):verify.sha(path) for path in verify.source_files()}
    assert set(KEYS)==set(actual)


@pytest.mark.parametrize("key",KEYS)
def test_report_field_mutation_is_rejected(key):
    actual=verify.build_report()
    changed=copy.deepcopy(actual)
    changed[key]="not the certified value"
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


@pytest.mark.parametrize("mutation",("missing","extra"))
def test_missing_or_extra_field_is_rejected(mutation):
    actual=verify.build_report()
    changed=copy.deepcopy(actual)
    if mutation=="missing":
        del changed["explicit_new_fixed_profile_and_old_profile_difference"]
    else:
        changed["full_UV_matching_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_identity_gate_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==170
    assert actual["checked_scalar_entries"]==186
    assert len(actual["proof_checks"])==46
    assert actual["controls"]["rejected_inputs"]==94
    assert len(actual["source_sha256"])==17
