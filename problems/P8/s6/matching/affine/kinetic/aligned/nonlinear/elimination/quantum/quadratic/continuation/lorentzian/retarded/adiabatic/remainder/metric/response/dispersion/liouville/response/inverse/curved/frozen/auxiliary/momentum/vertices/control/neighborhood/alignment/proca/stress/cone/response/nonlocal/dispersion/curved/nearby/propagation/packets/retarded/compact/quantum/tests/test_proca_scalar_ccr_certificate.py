"""Read-only scalar CCR certificate and every-field mutation rejection."""
import copy
import json

import pytest
from p8_proca_scalar_ccr import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation",
    "written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries",
    "proof_checks","literal_density_Hamiltonian_symplectic_map_and_source_normalization",
    "actual_all_momentum_transfer_growth_and_Schwartz_test_space",
    "explicit_positive_Gaussian_Cauchy_covariance_and_field_state",
    "state_independent_Kubo_response_and_compact_commutator_lower",
    "controls","verdict","not_established","verification_boundary"]


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()),verify.build_report())


def test_complete_source_manifest_and_report_keys_pinned():
    actual=verify.build_report()
    assert actual["source_sha256"]=={str(path.relative_to(verify.ROOT)):verify.sha(path) for path in verify.source_files()}
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
        del changed["state_independent_Kubo_response_and_compact_commutator_lower"]
    else:
        changed["Hadamard_and_UV_matching_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_frozen_native_identity_gate_control_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==22
    assert actual["checked_scalar_entries"]==175
    assert len(actual["proof_checks"])==12
    assert actual["controls"]["rejected_inputs"]==17
    assert len(actual["source_sha256"])==17
