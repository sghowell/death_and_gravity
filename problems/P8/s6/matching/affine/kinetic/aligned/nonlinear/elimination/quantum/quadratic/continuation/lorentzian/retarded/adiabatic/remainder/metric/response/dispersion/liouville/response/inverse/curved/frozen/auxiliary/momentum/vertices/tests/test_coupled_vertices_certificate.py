"""Read-only complete field mutation and source-bound certificate replay."""
import copy
import json

import pytest
from p8_coupled_vertices import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","actual_lapse_monomial_derivatives","original_boundary_and_Q_lapse_jets","actual_time_coefficient_bounds","physical_vertex_fixture_kernels","independent_electric_quartic","reality_and_permutation_fixtures","raw_physical_majorant_example","hard_wave_raw_majorant_example","vertex_domain_and_normalization","actual_stationary_formulas","controls","verdict","not_established","verification_boundary"]


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
        del changed["raw_physical_majorant_example"]
    else:
        changed["normalized_interacting_cutoff_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_identity_gate_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==154
    assert actual["checked_scalar_entries"]==455
    assert len(actual["proof_checks"])==17
    assert actual["controls"]["rejected_inputs"]==98
    assert len(actual["source_sha256"])==16
