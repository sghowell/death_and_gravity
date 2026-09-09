"""Independent read-only global support replay and every report-field mutation."""
import copy
import json

import pytest
from p8_proca_global_support import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","actual_retuned_global_clock_background_and_coefficient_substitution","complete_regular_density_phase_and_original_local_source","complete_unitary_and_gamma_Euler_systems","full_chart_intertwining_and_polynomial_endpoint_weights","explicit_global_chart_cover_and_complex_pole_separation","exact_principal_modes_and_real_frequency_growth_cancellation","new_global_quadratic_scalar_state_and_causal_Kubo_response","controls","verdict","not_established","verification_boundary"]


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
        del changed["explicit_global_chart_cover_and_complex_pole_separation"]
    else:
        changed["self_consistent_quantum_global_bounce_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_frozen_identity_scalar_gate_control_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==55
    assert actual["checked_scalar_entries"]==349
    assert len(actual["proof_checks"])==20
    assert actual["controls"]["rejected_inputs"]==54
    assert len(actual["source_sha256"])==20
