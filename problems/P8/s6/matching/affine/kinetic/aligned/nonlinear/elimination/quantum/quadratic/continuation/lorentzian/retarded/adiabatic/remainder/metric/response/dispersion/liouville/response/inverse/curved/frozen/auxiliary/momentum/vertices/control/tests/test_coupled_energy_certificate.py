"""Read-only complete field mutation and source-bound certificate replay."""
import copy
import json

import pytest
from p8_coupled_energy import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","actual_two_chart_scalar_matrices","actual_background_jet_enclosures","continuous_chart_margins","canonical_initial_shift_negative_control","normalized_scalar_column_domain","unchanged_selected_vector_columns","both_tensor_columns","finite_hard_tree_estimate","exact_phase_homothety","fixed_fiber_observable","controls","verdict","not_established","verification_boundary"]


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
        del changed["finite_hard_tree_estimate"]
    else:
        changed["Wilsonian_interacting_cutoff_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_identity_gate_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==32
    assert actual["checked_scalar_entries"]==71
    assert len(actual["proof_checks"])==51
    assert actual["controls"]["rejected_inputs"]==76
    assert len(actual["source_sha256"])==18
