"""Read-only complete field mutation and source-bound certificate replay."""
import copy
import json

import pytest
from p8_offclock_scalar import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","actual_constraint_compatible_homogeneous_family","continuous_offclock_auxiliary_domain","complete_six_channel_phase_Hessian","actual_canonical_gamma_blocks","actual_Lagrangian_high_frequency_orders","arbitrary_finite_time_jet_determinant_check","Euler_first_frozen_diagnostic","old_clock_negative_controls","canonical_and_frame_conventions","quantifiers_and_boundary","controls","verdict","not_established","verification_boundary"]


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
        del changed["Euler_first_frozen_diagnostic"]
    else:
        changed["actual_bounce_instability_or_UV_exclusion_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_identity_gate_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==45
    assert actual["checked_scalar_entries"]==97
    assert len(actual["proof_checks"])==21
    assert actual["controls"]["rejected_inputs"]==36
    assert len(actual["source_sha256"])==16
