"""Frozen field-by-field report replay with source pins."""
import copy
import json

import pytest
from p8_vector_dispersion_inverse import verify

KEYS = ["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","first_sheet_block_and_endpoint_matrices","exact_cut_density_and_positive_Gram","causal_kernel_endpoint_constants_and_representation","prepared_primitive_examples","function_space_boundary","variable_h_policy","controls","verdict","not_established","verification_boundary"]


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
        del changed["exact_cut_density_and_positive_Gram"]
    else:
        changed["full_curved_quantum_health_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_identity_gate_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==40
    assert actual["checked_scalar_entries"]==100
    assert len(actual["proof_checks"])==15
    assert actual["controls"]["rejected_inputs"]==38
    assert len(actual["source_sha256"])==14
