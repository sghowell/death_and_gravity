"""Full field replay and source-hash coverage."""
import copy
import json

import pytest
from p8_bounce_frequency import verify

KEYS=["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","physical_background_bounce_jets","constant_source_chart","frozen_tree_Euler_coefficients","frozen_local_Gaussian_Euler_coefficients","fixed_tadpole_Euler_contact","actual_pair_vertices_and_Gram","continuous_frequency_bracket_bounds","diagnostic_definition_and_domain","controls","verdict","not_established","verification_boundary"]


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
        del changed["continuous_frequency_bracket_bounds"]
    else:
        changed["actual_bounce_instability_proved"]=True
    with pytest.raises(ValueError,match="differs"):
        verify.validate_report(changed,actual)


def test_identity_gate_and_source_counts():
    actual=verify.build_report()
    assert actual["named_exact_check_count"]==24
    assert actual["checked_scalar_entries"]==64
    assert len(actual["proof_checks"])==26
    assert actual["controls"]["rejected_inputs"]==28
    assert len(actual["source_sha256"])==12
