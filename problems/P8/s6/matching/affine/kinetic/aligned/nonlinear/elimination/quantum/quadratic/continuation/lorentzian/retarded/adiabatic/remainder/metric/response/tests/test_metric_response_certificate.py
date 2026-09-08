"""Read-only finite prepared metric report and full field mutation coverage."""
import copy
import json

import pytest
from p8_vector_metric_response import verify

KEYS = ["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","eighth_order_two_source_coefficient_bounds","varied_reference_residual_and_normalization_bounds","physical_readout_vertices","varied_physical_readout_tail_bounds","exact_mixing_and_reference_constants","finite_vector_metric_response_scale_example","existing_tadpole_physical_vertices","complete_background_cancelled_response_scale_example","initial_mixing_and_source_norm_controls","reference_mode_policy","dimensional_limit","prepared_source_domain","controls","verdict","not_established","verification_boundary"]


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
        del changed["prepared_source_domain"]
    else:
        changed["no_loss_full_quantum_feedback_proved"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_and_control_counts():
    data = verify.build_report()
    assert data["named_exact_check_count"] == 59
    assert data["checked_scalar_entries"] == 69
    assert len(data["proof_checks"]) == 38
    assert data["controls"]["rejected_inputs"] == 204
    assert len(data["source_sha256"]) == 17
