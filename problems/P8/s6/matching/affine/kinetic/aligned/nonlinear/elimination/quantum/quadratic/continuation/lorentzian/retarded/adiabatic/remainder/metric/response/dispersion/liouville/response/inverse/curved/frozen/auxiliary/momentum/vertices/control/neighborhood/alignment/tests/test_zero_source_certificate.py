"""Read-only new-action certificate, pinned sources and complete field mutation."""
import copy
import json

import pytest
from p8_zero_source import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","literal_new_action_and_Hamiltonian","all_connection_affine_lift","actual_new_five_lapse_jets_for_fourteen_inventory_rows","continuous_new_auxiliary_domain","universal_new_coefficient_bounds_and_positive_scaling","finite_same_scale_hard_tree_change","fixed_light_Gaussian_and_changed_source_functional","old_offclock_mechanism_comparison","canonical_comparison_convention","controls","verdict","not_established","verification_boundary"]


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
        del changed["finite_same_scale_hard_tree_change"]
    else:
        changed["full_UV_matching_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_identity_gate_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==94
    assert actual["checked_scalar_entries"]==406
    assert len(actual["proof_checks"])==31
    assert actual["controls"]["rejected_inputs"]==34
    assert len(actual["source_sha256"])==20
