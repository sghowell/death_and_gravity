"""Read-only certificate, exact identities and non-overclaim controls."""

import copy
import json
from functools import cache

import pytest
import sympy as sp
from p8a_maxwell_cosmology import calibration, controls, dictionary, verify


@cache
def candidate():
    return verify.build_report()


def test_frozen_certificate_replays_without_writing_any_file():
    before = verify.REPORT.read_bytes()
    report = json.loads(before)
    verify.validate_report(report, candidate())
    assert verify.REPORT.read_bytes() == before


def test_exact_residual_count_and_scope_flags():
    report = candidate()
    assert sum(len(group) for group in report["exact_residuals"].values()) == 35
    assert report["theorem"]["future_caps_derived_from_history_or_QSEI"] is False
    assert report["scope_assessment"]["root_P8_completion_declared"] is False
    assert report["scope_assessment"]["actual_SEE_witness_required_for_conditional_theorem"] is False
    assert report["geometric_nonvacuity_controls"]["future_timelike_and_null_complete"]
    assert not report["geometric_nonvacuity_controls"]["actual_allowed_Maxwell_SEE_solution"]


def test_each_new_source_hash_and_reserved_audit_are_included():
    hashes = candidate()["source_sha256"]
    assert "tests/test_maxwell_cosmology_covariant_audit.py" in hashes
    assert "notes/closure.md" in hashes
    for name, expected in hashes.items():
        assert verify.sha(verify.ROOT/name) == expected


@pytest.mark.parametrize("section,key,value", [
    ("theorem", "strict_margin_lower", "1/4"),
    ("theorem", "future_caps_derived_from_history_or_QSEI", True),
    ("observer_and_physical_dictionary", "actual_H0_times_tau_unanchored", ["2", "2"]),
    ("geometric_nonvacuity_controls", "actual_allowed_Maxwell_SEE_solution", True),
    ("scope_assessment", "root_P8_completion_declared", True),
])
def test_numerical_or_claim_promotions_fail_full_replay(section, key, value):
    damaged = copy.deepcopy(candidate())
    damaged[section][key] = value
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(damaged, candidate())


@pytest.mark.parametrize("bad", [True, 0.1, sp.oo, sp.nan, sp.Symbol("unknown")])
def test_all_exact_numerical_gate_inputs_reject_ambiguous_values(bad):
    for function in (lambda: calibration.cost_at_beta(bad),
                     lambda: calibration.theorem_gate(bad, 0),
                     lambda: dictionary.proper_scales(bad),
                     lambda: controls.complete_geometry(bad)):
        with pytest.raises((TypeError, ValueError)):
            function()


def test_theorem_parent_pin_is_required_before_replay(monkeypatch):
    verify.prior_checks.cache_clear()
    try:
        monkeypatch.setattr(verify, "PARENT_SHA", "0"*64)
        with pytest.raises(ValueError, match="pinned"):
            verify.prior_checks()
    finally:
        verify.prior_checks.cache_clear()
