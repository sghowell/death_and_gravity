"""Read-only sharper response certificate and complete field-mutation tests."""
import copy
import json

import pytest
from p8_proca_sharp_response import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation",
    "written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries",
    "proof_checks","actual_complex_whole_box_logarithmic_jets",
    "small_exact_action_basis_and_complete_source_output_factor",
    "complete_literal_finite_frequency_Laurent_majorants",
    "sharper_all_frequency_transfer_and_two_endpoint_normal_form",
    "new_rational_compact_probes_complete_remainder_and_spatial_band",
    "controls","verdict","not_established","verification_boundary"]


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()),verify.build_report())


def test_all_source_files_and_report_fields_pinned():
    actual=verify.build_report()
    assert actual["source_sha256"]=={str(p.relative_to(verify.ROOT)):verify.sha(p) for p in verify.source_files()}
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
        del changed["new_rational_compact_probes_complete_remainder_and_spatial_band"]
    else:
        changed["interacting_parent_matching_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_frozen_identity_scalar_gate_control_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==25
    assert actual["checked_scalar_entries"]==70
    assert len(actual["proof_checks"])==76
    assert actual["controls"]["rejected_inputs"]==66
    assert len(actual["source_sha256"])==17
