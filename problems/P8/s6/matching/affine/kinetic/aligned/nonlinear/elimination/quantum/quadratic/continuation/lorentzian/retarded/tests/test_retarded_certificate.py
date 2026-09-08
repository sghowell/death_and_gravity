"""Read-only retarded mass-response report and mutation checks."""
import copy
import json

import pytest
from p8_vector_retarded import verify

KEYS = ["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","unchanged_clock_profiles","constraint_and_contact","canonical_insertions","general_momentum_connected_integrands","mode_normalization","retarded_formula","canonical_covariance_identity","flat_frequency_kernel","flat_pole_and_omitted_contact_controls","prepared_source_and_causal_boundary","controls","verdict","not_established","verification_boundary"]


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_sources_and_report_keys_pinned():
    data = verify.build_report()
    assert data["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}
    assert set(KEYS) == set(data)


@pytest.mark.parametrize("key", KEYS)
def test_report_mutation_rejected(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_or_extra_report_field_rejected(mutation):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    if mutation == "missing":
        del changed["retarded_formula"]
    else:
        changed["finite_off_clock_cone_proved"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_and_control_counts():
    data = verify.build_report()
    assert data["named_exact_check_count"] == 65
    assert data["checked_scalar_entries"] == 147
    assert len(data["proof_checks"]) == 13
    assert data["controls"]["rejected_inputs"] == 140
    assert len(data["source_sha256"]) == 16
