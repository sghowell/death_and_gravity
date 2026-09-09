"""Full ancestor replay, exact source manifest and every-field mutation."""
import copy
import json

import pytest
from p8_proca_finite_width import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","actual_all_real_frequency_chart","actual_initial_state_and_eighth_order_coefficients","whole_interval_jet_and_remainder_bounds","explicit_actual_finite_width_inequality","controls","verdict","not_established","verification_boundary"]

def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()),verify.build_report())

def test_every_source_and_report_field_is_pinned():
    actual=verify.build_report()
    assert actual["source_sha256"]=={str(p.relative_to(verify.ROOT)):verify.sha(p) for p in verify.source_files()}
    assert set(actual)==set(KEYS)

@pytest.mark.parametrize("key",KEYS)
def test_every_report_field_mutation_rejected(key):
    actual=verify.build_report()
    changed=copy.deepcopy(actual)
    changed[key]="not the certified value"
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)

@pytest.mark.parametrize("mutation",("missing","extra"))
def test_missing_or_extra_report_fields_rejected(mutation):
    actual=verify.build_report()
    changed=copy.deepcopy(actual)
    if mutation=="missing":
        del changed["whole_interval_jet_and_remainder_bounds"]
    else:
        changed["complete_covariant_physical_stress_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)

def test_frozen_identity_scalar_gate_control_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==30
    assert actual["checked_scalar_entries"]==261
    assert len(actual["proof_checks"])==17
    assert actual["controls"]["rejected_inputs"]==38
    assert len(actual["source_sha256"])==20
