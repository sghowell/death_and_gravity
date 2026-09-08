"""Pinned report replay and independent mutation of every report field."""
import copy
import json

import pytest
from p8_prepared_volterra import verify

KEYS = ["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","physical_current_UV_coefficients","full_dimensional_fourth_adiabatic_contacts","literal_classical_matter_reduction","local_four_primitive_kernels","singular_power_four_primitive_shapes","function_space_estimates_and_domains","coupling_and_preparation_policy","controls","verdict","not_established","verification_boundary"]


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()),verify.build_report())


def test_all_sources_and_report_keys_pinned():
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
        del changed["full_dimensional_fourth_adiabatic_contacts"]
    else:
        changed["quantum_stability_and_VGB_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_identity_gate_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==47
    assert actual["checked_scalar_entries"]==72
    assert len(actual["proof_checks"])==18
    assert actual["controls"]["rejected_inputs"]==40
    assert len(actual["source_sha256"])==15
