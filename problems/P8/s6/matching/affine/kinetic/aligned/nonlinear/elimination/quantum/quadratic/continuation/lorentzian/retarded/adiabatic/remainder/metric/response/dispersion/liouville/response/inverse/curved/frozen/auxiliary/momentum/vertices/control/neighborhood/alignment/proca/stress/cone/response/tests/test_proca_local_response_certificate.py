"""Read-only local ordinary Proca report, pinned sources and complete field mutation."""
import copy
import json

import pytest
from p8_proca_local_response import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","new_finite_covariant_heat_quadratic_densities","actual_physical_coordinate_Euler_currents","actual_normalized_local_physical_stress_response","actual_nonlinear_clock_chart_derivatives_and_time_jets","actual_metric_pullback_and_fixed_local_profile_densities","actual_clock_chart_local_Euler_currents","new_rank_one_fourth_derivative_coefficient_and_null_direction","continuous_local_response_coefficient_envelopes","same_scale_new_local_response_C4_to_C0_bounds","controls","verdict","not_established","verification_boundary"]


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
        del changed["new_rank_one_fourth_derivative_coefficient_and_null_direction"]
    else:
        changed["full_UV_matching_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_identity_gate_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==246
    assert actual["checked_scalar_entries"]==250
    assert len(actual["proof_checks"])==15
    assert actual["controls"]["rejected_inputs"]==100
    assert len(actual["source_sha256"])==16

