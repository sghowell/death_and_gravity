"""Read-only new-action certificate and complete report-field mutation tests."""
import copy
import json

import pytest
from p8_proca_retuned_margin import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","literal_new_covariant_action_and_fixed_phase_lapse_change","new_exact_rational_constraint_and_flow_signatures","fresh_whole_real_and_complex_coordinate_enclosures","complete_physical_proper_time_Hubble_acceleration","new_Picard_solution_conserved_charge_and_principal_health","preserved_global_clock_and_literal_deformation_budgets","central_positive_matter_continuum_and_weaker_margin_control","controls","verdict","not_established","verification_boundary"]


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
        del changed["central_positive_matter_continuum_and_weaker_margin_control"]
    else:
        changed["interacting_parent_matching_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_frozen_identity_scalar_gate_control_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==62
    assert actual["checked_scalar_entries"]==62
    assert len(actual["proof_checks"])==54
    assert actual["controls"]["rejected_inputs"]==69
    assert len(actual["source_sha256"])==18
