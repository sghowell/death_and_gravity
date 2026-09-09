"""Read-only actual prepared Proca inverse report and complete field mutation."""
import copy
import json

import pytest
from p8_proca_prepared_inverse import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","actual_new_time_covariance_and_frozen_profile_currents","actual_literal_tree_zero_charge_reduction_and_all_derivative_coefficients","new_full_dimensional_current_pairs","new_curved_scalar_leading_log_and_fixed_finite_contact_matching","actual_prepared_block_inverse_and_original_force_reconstruction","exact_fixed_mass_coupling_normalizations","controls","verdict","not_established","verification_boundary"]


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
        del changed["actual_prepared_block_inverse_and_original_force_reconstruction"]
    else:
        changed["full_UV_matching_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_identity_gate_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==116
    assert actual["checked_scalar_entries"]==154
    assert len(actual["proof_checks"])==25
    assert actual["controls"]["rejected_inputs"]==39
    assert len(actual["source_sha256"])==17
