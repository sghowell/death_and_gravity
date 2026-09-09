"""Read-only new-action certificate, pinned sources and complete field mutation."""
import copy
import json

import pytest
from p8_constant_proca import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","literal_affine_mass_source_and_contact_update","actual_new_Hamiltonian","actual_new_lapse_jet_inventory","continuous_new_auxiliary_and_coefficient_bounds","same_scale_finite_hard_tree_change","constrained_physical_Proca_symbol","actual_nearby_vector_canonical_blocks","changed_quantum_readouts_and_operator_jets","new_constant_mass_local_quantum_coefficients_and_boundary","controls","verdict","not_established","verification_boundary"]


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
        del changed["same_scale_finite_hard_tree_change"]
    else:
        changed["full_UV_matching_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_identity_gate_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==85
    assert actual["checked_scalar_entries"]==1330
    assert len(actual["proof_checks"])==23
    assert actual["controls"]["rejected_inputs"]==34
    assert len(actual["source_sha256"])==18
