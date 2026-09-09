"""Read-only spectator screening certificate and complete field-mutation tests."""
import copy
import json

import pytest
from p8_proca_spectator_repair import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","actual_fixed_physical_matter_frame_light_canonical_forms","finite_Hermitian_principal_extension_and_explicit_embedding","necessary_two_orientation_light_repair_and_negative_controls","actual_ordinary_Proca_constraint_decoupling_and_finite_species","controls","verdict","not_established","verification_boundary"]


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
        del changed["necessary_two_orientation_light_repair_and_negative_controls"]
    else:
        changed["interacting_parent_matching_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_frozen_identity_scalar_gate_control_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==31
    assert actual["checked_scalar_entries"]==40
    assert len(actual["proof_checks"])==16
    assert actual["controls"]["rejected_inputs"]==64
    assert len(actual["source_sha256"])==15
