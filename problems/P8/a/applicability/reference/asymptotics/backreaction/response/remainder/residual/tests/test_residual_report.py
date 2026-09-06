import copy
import json

import pytest
from p8a_residual import verify


@pytest.fixture(scope="module")
def replay():
    return verify.build_report()


def test_residual_new_certificate_replays_without_writing(replay):
    before = verify.REPORT.read_bytes()
    verify.validate_report(json.loads(before), replay)
    assert verify.REPORT.read_bytes() == before
    assert replay["prior_A7_sha256"] == verify.PRIOR_SHA


def test_residual_parent_covariant_audit_and_independent_assembly_are_pinned(replay):
    assert "tests/test_residual_covariant_audit.py" in replay["source_sha256"]
    assert "src/p8a_residual/independent.py" in replay["source_sha256"]


@pytest.mark.parametrize("field,value", [("status", "EXACT_SEE_SOLUTION"),
    ("prior_A7_sha256", "changed"), ("source_sha256", {}), ("verification_boundary", [])])
def test_residual_hash_and_scope_tampering_rejected(replay, field, value):
    changed = copy.deepcopy(replay)
    changed[field] = value
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, replay)


@pytest.mark.parametrize("field", ["generic_gamma_or_other_lambda_covered_by_numeric_constants",
    "absolute_stress_scheme_independent", "positive_EED_for_arbitrary_Hadamard_states",
    "extra_curvature_sources_automatically_bounded", "exact_or_nearby_SEE_solution_or_stability",
    "self_consistent_past_initial_data_certified", "small_residual_on_entire_off_shell_preparation_interval",
    "second_order_corrected_metric_family_certified", "perturbed_QSEI_or_all_sampler_bound",
    "all_geodesic_focusing_or_P8a_completion"])
def test_residual_unearned_promotion_rejected(replay, field):
    changed = copy.deepcopy(replay)
    changed["scope"][field] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, replay)


def test_residual_literal_frozen_denominator_tampering_rejected(replay):
    changed = copy.deepcopy(replay)
    changed["independently_replayed_finite_bounds"]["literal_frozen_denominator"] = 1922**2
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, replay)


def test_residual_actual_quantum_bound_cannot_be_changed(replay):
    changed = copy.deepcopy(replay)
    changed["independently_replayed_finite_bounds"]["first_constants"]["density"] = "1"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, replay)


def test_residual_prior_certificate_hash_tampering_rejected(tmp_path, monkeypatch):
    changed = tmp_path/"changed.json"
    changed.write_text("{}")
    monkeypatch.setattr(verify.prior, "REPORT", changed)
    with pytest.raises(ValueError, match="Pinned A.7"):
        verify.prior_checks()


def test_residual_nonzero_identity_and_lost_negative_control_rejected():
    with pytest.raises(ValueError, match="Nonzero"):
        verify.verified_residuals({"wrong": 1})
    with pytest.raises(ValueError, match="vanished"):
        verify.verified_exclusion("wrong", 0)
