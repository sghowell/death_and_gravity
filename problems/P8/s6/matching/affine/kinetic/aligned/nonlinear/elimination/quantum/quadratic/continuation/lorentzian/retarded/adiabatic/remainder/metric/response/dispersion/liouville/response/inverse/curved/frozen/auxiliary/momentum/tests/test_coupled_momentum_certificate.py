"""Read-only source-bound replay and mutation of every certificate field."""
import copy
import json

import pytest
from p8_coupled_momentum import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","canonical_matter_boundary","physical_momentum_definition","exact_York_operator_inverse_and_norms","mixed_fixture_metadata","omission_negative_controls","full_quadratic_phase_channels","full_quadratic_phase_matrix","actual_linear_lapse_invariant_coefficients","compact_scalar_Hamiltonian","nonzero_transfer_and_quadratic_scope","controls","verdict","not_established","verification_boundary"]


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
        del changed["full_quadratic_phase_matrix"]
    else:
        changed["physical_interacting_cutoff_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_identity_gate_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==147
    assert actual["checked_scalar_entries"]==273
    assert len(actual["proof_checks"])==17
    assert actual["controls"]["rejected_inputs"]==65
    assert len(actual["source_sha256"])==15
