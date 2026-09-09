"""Independent report replay and mutation of every declared field."""
import copy
import json

import pytest
from p8_proca_global_hadamard import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","physical_volume_operator_and_Green_CCR_conventions","fresh_complete_global_canonical_Laurent_chart","actual_initial_slab_graph_and_center_jets","all_order_Riccati_inverse_and_four_order_native_tests","exact_positive_graph_Cauchy_covariance_and_infrared_completion","explicit_all_finite_strip_frequency_gaps","complete_order_zero_chart_transitions","global_two_cone_wavefront_condition","controls","verdict","not_established","verification_boundary"]


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
        del changed["global_two_cone_wavefront_condition"]
    else:
        changed["renormalized_physical_stress_and_original_P8_closed"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_frozen_identity_scalar_gate_control_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==65
    assert actual["checked_scalar_entries"]==599
    assert len(actual["proof_checks"])==19
    assert actual["controls"]["rejected_inputs"]==64
    assert len(actual["source_sha256"])==22
