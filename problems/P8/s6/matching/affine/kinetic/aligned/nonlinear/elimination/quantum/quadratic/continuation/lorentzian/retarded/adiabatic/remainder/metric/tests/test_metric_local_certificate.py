"""Read-only full metric local report and field mutation coverage."""
import copy
import json

import pytest
from p8_vector_metric_local import verify

KEYS = ["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","full_metric_vertices","finite_Euler_current_coefficient_rows","physical_stress_coefficient_majorants","finite_local_scale_example","omitted_counterterm_failure_controls","second_mass_vertex_control","fourth_derivative_local_symbol_control","dimensional_prescription","source_and_output_domain","controls","verdict","not_established","verification_boundary"]


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
        del changed["source_and_output_domain"]
    else:
        changed["full_nonlocal_metric_response_proved"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_and_control_counts():
    data = verify.build_report()
    assert data["named_exact_check_count"] == 132
    assert data["checked_scalar_entries"] == 135
    assert len(data["proof_checks"]) == 22
    assert data["controls"]["rejected_inputs"] == 126
    assert len(data["source_sha256"]) == 19
