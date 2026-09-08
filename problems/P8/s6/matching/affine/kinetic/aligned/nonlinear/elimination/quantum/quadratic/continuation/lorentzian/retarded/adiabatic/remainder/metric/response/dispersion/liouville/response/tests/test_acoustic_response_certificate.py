"""Read-only report replay and complete report-field mutation coverage."""
import copy
import json

import pytest
from p8_vector_acoustic_response import verify

KEYS = ["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","generic_covariance_transport","prepared_sector_covariance_sources","physical_readout_response_and_contact_rows","generic_physical_readout_variation","full_dimensional_pump_squares","retarded_initial_policy","finite_limit_policy","controls","verdict","not_established","verification_boundary"]


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_source_files_and_report_keys_pinned():
    data = verify.build_report()
    assert data["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}
    assert set(KEYS) == set(data)


@pytest.mark.parametrize("key", KEYS)
def test_report_field_mutation_is_rejected(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_or_extra_field_is_rejected(mutation):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    if mutation == "missing":
        del changed["physical_readout_response_and_contact_rows"]
    else:
        changed["coupled_quantum_inverse_proved"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_gate_and_source_counts():
    data = verify.build_report()
    assert data["named_exact_check_count"] == 46
    assert data["checked_scalar_entries"] == 82
    assert len(data["proof_checks"]) == 13
    assert data["controls"]["rejected_inputs"] == 40
    assert len(data["source_sha256"]) == 15
