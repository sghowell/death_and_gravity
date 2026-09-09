"""Independent complete replay and every report-field mutation."""
import copy
import json

import pytest
from p8_proca_relational_qei import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","actual_normalized_modes_and_nonzero_relational_clock_residue","actual_worldline_domain_and_interluminal_control","finite_reference_difference_QEI_and_short_sampling_coefficients","fixed_background_test_stress_and_actual_coupled_source","controls","verdict","not_established","verification_boundary"]


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()),verify.build_report())


def test_all_source_files_and_report_fields_pinned():
    actual=verify.build_report()
    assert actual["source_sha256"]=={str(p.relative_to(verify.ROOT)):verify.sha(p) for p in verify.source_files()}
    assert set(KEYS)==set(actual)


@pytest.mark.parametrize("key",KEYS)
def test_every_report_field_mutation_rejected(key):
    actual=verify.build_report()
    changed=copy.deepcopy(actual)
    changed[key]="not the certified value"
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


@pytest.mark.parametrize("mutation",("missing","extra"))
def test_missing_or_extra_report_field_rejected(mutation):
    actual=verify.build_report()
    changed=copy.deepcopy(actual)
    if mutation=="missing":
        del changed["fixed_background_test_stress_and_actual_coupled_source"]
    else:
        changed["optimal_cosmological_physical_stress_QEI_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_frozen_identity_scalar_gate_control_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==31
    assert actual["checked_scalar_entries"]==56
    assert len(actual["proof_checks"])==17
    assert actual["controls"]["rejected_inputs"]==66
    assert len(actual["source_sha256"])==18
