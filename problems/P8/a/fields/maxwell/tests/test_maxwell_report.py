"""Read-only certificate, parent pin and false-promotion rejection tests."""

import copy
import hashlib
import json
from functools import cache

import pytest
from p8a_maxwell import verify


@cache
def replayed():
    return verify.build_report()


def test_certificate_exact_replay_and_no_certificate_mutation():
    before = hashlib.sha256(verify.REPORT.read_bytes()).hexdigest()
    verify.validate_report(json.loads(verify.REPORT.read_text()), replayed())
    assert hashlib.sha256(verify.REPORT.read_bytes()).hexdigest() == before


def test_old_geometric_and_clock_reports_are_exactly_pinned():
    assert verify.sha(verify.geometric_prior.REPORT) == verify.PINS["A1"]
    assert verify.sha(verify.clock_prior.REPORT) == verify.PINS["A12"]
    assert replayed()["prior_sha256"] == verify.PINS


def test_every_new_source_and_reserved_covariant_audit_is_hashed():
    source_hashes = replayed()["source_sha256"]
    assert "tests/test_maxwell_covariant_audit.py" in source_hashes
    assert "src/p8a_maxwell/independent.py" in source_hashes
    for name, digest in source_hashes.items():
        assert verify.sha(verify.ROOT/name) == digest


@pytest.mark.parametrize("key", ["global_target_state_extension_assumed",
                                 "finite_beta_inherited_from_scalar_gamma",
                                 "optimal_coefficient_claimed"])
def test_field_scope_false_promotions_are_rejected(key):
    wrong = copy.deepcopy(replayed())
    wrong["theorem"][key] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(wrong, replayed())


@pytest.mark.parametrize("key", ["new_SEE_solution_asserted",
                                 "initial_pointwise_Ricci_requirement_removed",
                                 "all_A1_extended_normal_windows_verified",
                                 "naive_small_quantum_constant_example_focuses"])
def test_conditional_focusing_is_not_promoted_to_a_cosmological_witness(key):
    wrong = copy.deepcopy(replayed())
    wrong["conditional_geometric_dictionary"][key] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(wrong, replayed())


def test_reference_anomaly_or_extra_source_cannot_disappear_from_report():
    wrong = copy.deepcopy(replayed())
    wrong["derived_constants"]["radiation_control"]["reference_EED_pi_squared_t_fourth_over_hbar"] = "0"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(wrong, replayed())
    wrong = copy.deepcopy(replayed())
    wrong["derived_constants"]["source_controls"]["Lambda=1,other_lower=-2"]["coarsened_source"] = "1"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(wrong, replayed())


def test_fraction_disagreement_fails_the_actual_bridge(monkeypatch):
    expected = verify.independent.replay()
    expected["radiation_control"]["absolute_zero_times_t_fourth"] = "0"
    monkeypatch.setattr(verify.independent, "replay", lambda: expected)
    with pytest.raises(ValueError, match="independent Fraction"):
        verify.checked_constants()


def test_changed_parent_pin_fails_before_replay(monkeypatch):
    monkeypatch.setitem(verify.PINS, "A12", "0"*64)
    verify.prior_checks.cache_clear()
    with pytest.raises(ValueError, match="pinned"):
        verify.prior_checks()
    verify.prior_checks.cache_clear()
