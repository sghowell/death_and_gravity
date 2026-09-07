"""Read-only replay and rejection of missing physical/geometric hypotheses."""

import copy
import hashlib
import json
from functools import cache

import pytest
from p8a_maxwell_focusing import verify


@cache
def replayed():
    return verify.build_report()


def test_stored_report_exact_replay_is_read_only():
    before = hashlib.sha256(verify.REPORT.read_bytes()).hexdigest()
    verify.validate_report(json.loads(verify.REPORT.read_text()), replayed())
    assert hashlib.sha256(verify.REPORT.read_bytes()).hexdigest() == before


def test_actual_Maxwell_parent_and_every_new_source_are_pinned():
    assert verify.sha(verify.prior.REPORT) == verify.PARENT_SHA
    hashes = replayed()["source_sha256"]
    assert "tests/test_maxwell_focusing_covariant_audit.py" in hashes
    for name, value in hashes.items():
        assert verify.sha(verify.ROOT/name) == value


@pytest.mark.parametrize("key,value", [
    ("future_caps_derived_from_initial_patch_or_QSEI", True),
    ("quantitative_contraction_history_assumed", False),
    ("initial_pointwise_SEC_or_Ricci_sign_assumed", True),
    ("compact_Cauchy_theorem_imported_on_R3", True),
    ("larger_extension_or_curvature_inextendibility_proved", True),
])
def test_hypothesis_or_scope_mutations_fail(key, value):
    wrong = copy.deepcopy(replayed())
    wrong["theorem"][key] = value
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(wrong, replayed())


def test_complete_geometric_control_cannot_be_promoted_to_a_quantum_solution():
    wrong = copy.deepcopy(replayed())
    wrong["countercontrol"]["actual_allowed_Maxwell_SEE_state_witness"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(wrong, replayed())


def test_an_actual_cosmology_or_EFT_claim_is_not_certified():
    wrong = copy.deepcopy(replayed())
    wrong["named_calibration"]["observational_or_fundamental_EFT_validity_verified"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(wrong, replayed())


def test_independent_fraction_bridge_rejects_tampered_cost(monkeypatch):
    wrong = verify.independent.replay(json.loads(verify.prior.REPORT.read_text()))
    wrong["constants"]["past_cost"] = "0"
    monkeypatch.setattr(verify.independent, "replay", lambda parent: wrong)
    with pytest.raises(ValueError, match="independent Fraction"):
        verify.checked_constants()


def test_changed_parent_hash_is_rejected_before_replay(monkeypatch):
    monkeypatch.setattr(verify, "PARENT_SHA", "0"*64)
    verify.prior_checks.cache_clear()
    with pytest.raises(ValueError, match="pinned"):
        verify.prior_checks()
    verify.prior_checks.cache_clear()
