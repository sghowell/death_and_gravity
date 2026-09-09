"""Read-only rank-one ordinary Proca report, pinned sources and complete field mutation."""
import copy
import json

import pytest
from p8_proca_rank_one_inverse import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","actual_new_physical_flat_pair_and_fixed_canonical_contacts","new_exact_rank_one_massive_scalar_block_and_radial_representation","full_massive_positive_axis_and_removable_zero_controls","new_positive_scalar_inverse_cut_density_and_exact_moments","exact_open_cut_controls","new_causal_scalar_L1_kernel_and_primitive","positive_fixed_mass_moment_and_primitive_bounds","controls","verdict","not_established","verification_boundary"]


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
        del changed["new_causal_scalar_L1_kernel_and_primitive"]
    else:
        changed["full_UV_matching_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_identity_gate_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==50
    assert actual["checked_scalar_entries"]==67
    assert len(actual["proof_checks"])==21
    assert actual["controls"]["rejected_inputs"]==78
    assert len(actual["source_sha256"])==16
