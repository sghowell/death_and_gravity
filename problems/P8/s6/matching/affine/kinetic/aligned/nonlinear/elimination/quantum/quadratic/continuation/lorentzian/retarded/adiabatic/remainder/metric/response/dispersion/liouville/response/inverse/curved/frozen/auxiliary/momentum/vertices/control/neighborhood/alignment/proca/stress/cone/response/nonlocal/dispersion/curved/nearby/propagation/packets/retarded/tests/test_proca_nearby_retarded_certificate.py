"""Read-only retarded-source certificate and complete field-mutation tests."""
import copy
import json

import pytest
from p8_proca_nearby_retarded import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","literal_local_relational_source_and_complete_canonical_maps","uniform_two_time_transfer_with_both_endpoint_changes","exact_clock_and_matter_front_amplitudes","complete_frequency_Fourier_and_L2_remainder_bounds","actual_matter_spacelike_compact_probe_neighborhoods","controls","verdict","not_established","verification_boundary"]


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
        del changed["complete_frequency_Fourier_and_L2_remainder_bounds"]
    else:
        changed["full_UV_matching_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_identity_gate_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==30
    assert actual["checked_scalar_entries"]==78
    assert len(actual["proof_checks"])==30
    assert actual["controls"]["rejected_inputs"]==85
    assert len(actual["source_sha256"])==18
