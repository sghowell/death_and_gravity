"""Read-only scalar-cone report, pinned sources and complete field mutation."""
import copy
import json

import pytest
from p8_proca_scalar_cone import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","actual_homogeneous_background_time_jets","full_spatial_two_leg_time_jet_reduction","actual_canonical_and_Euler_first_blocks","physical_principal_matrices_and_characteristic_factorization","continuous_central_positive_subluminal_domain","actual_finite_amplitude_luminality_boundary_and_superluminal_control","controls","verdict","not_established","verification_boundary"]


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
        del changed["physical_principal_matrices_and_characteristic_factorization"]
    else:
        changed["full_UV_matching_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_identity_gate_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==46
    assert actual["checked_scalar_entries"]==91
    assert len(actual["proof_checks"])==32
    assert actual["controls"]["rejected_inputs"]==68
    assert len(actual["source_sha256"])==17

