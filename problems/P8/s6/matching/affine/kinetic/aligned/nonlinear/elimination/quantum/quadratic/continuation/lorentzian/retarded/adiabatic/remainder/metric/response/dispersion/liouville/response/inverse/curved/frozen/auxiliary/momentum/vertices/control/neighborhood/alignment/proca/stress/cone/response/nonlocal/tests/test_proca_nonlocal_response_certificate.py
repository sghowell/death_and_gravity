"""Read-only prepared ordinary Proca report, pinned sources and complete field mutation."""
import copy
import json

import pytest
from p8_proca_nonlocal_response import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","actual_new_mode_source_vertices","actual_new_normalized_physical_readout_vertices","actual_new_varied_adiabatic_readouts","new_reference_low_residuals_and_full_tenth_source_coefficients","new_continuous_varied_reference_envelopes","new_integrable_reference_readout_tails","new_transport_constants_with_original_nonzero_initial_mixing","new_complete_physical_vector_response_bounds","actual_already_fixed_profile_vertices_and_total_chart_contact","continuous_actual_clock_chart_C10_lift","new_complete_background_cancelled_response_bounds_in_both_charts","controls","verdict","not_established","verification_boundary"]


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
        del changed["new_complete_background_cancelled_response_bounds_in_both_charts"]
    else:
        changed["full_UV_matching_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_identity_gate_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==112
    assert actual["checked_scalar_entries"]==133
    assert len(actual["proof_checks"])==41
    assert actual["controls"]["rejected_inputs"]==146
    assert len(actual["source_sha256"])==19

