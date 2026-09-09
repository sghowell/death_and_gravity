"""Read-only finite-frequency packet certificate and complete field mutation."""
import copy
import json

import pytest
from p8_proca_nearby_packets import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","exact_finite_frequency_canonical_system","action_normalized_modes_and_complete_transport","actual_complex_solution_domain_and_annuli","literal_finite_frequency_Laurent_coefficient_majorants","explicit_analytic_basis_and_inverse_majorants","actual_near_identity_normalform_and_uniform_constants","local_relational_observable_pole_and_reconstruction","exact_three_dimensional_classical_packets_and_normalized_errors","controls","verdict","not_established","verification_boundary"]


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
        del changed["exact_three_dimensional_classical_packets_and_normalized_errors"]
    else:
        changed["full_UV_matching_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_identity_gate_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==45
    assert actual["checked_scalar_entries"]==273
    assert len(actual["proof_checks"])==51
    assert actual["controls"]["rejected_inputs"]==137
    assert len(actual["source_sha256"])==19
