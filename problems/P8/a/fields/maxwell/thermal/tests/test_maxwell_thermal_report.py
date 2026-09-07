"""Pinned read-only replay, identity inventory and claim-omission controls."""

import copy
import json
from functools import cache

import pytest
from p8a_maxwell_thermal import verify


@cache
def candidate():
    return verify.build_report()


def test_frozen_thermal_certificate_replays_without_writing():
    before = verify.REPORT.read_bytes()
    verify.validate_report(json.loads(before), candidate())
    assert verify.REPORT.read_bytes() == before


def test_exact_identity_inventory_and_transitive_source_pin():
    report = candidate()
    assert sum(map(len, report["exact_residuals"].values())) == 27
    assert report["prior_sha256"] == {"A18": verify.PARENT_SHA, "A16_transitive": verify.PHOTON_SHA}
    assert report["actual_history_embedding"]["actual_SEE_history_and_state_realization"]
    assert not report["named_source_and_solution"]["order_reduced_surrogate_or_only_formal_solution"]
    assert report["relationship_to_QSEI_objective"]["actual_history_realization_is_optional_strengthening"]
    assert not report["relationship_to_QSEI_objective"]["reopens_scoped_P8a_completion"]


def test_all_source_hashes_include_reserved_covariant_audit_unchanged():
    hashes = candidate()["source_sha256"]
    assert "tests/test_maxwell_thermal_covariant_audit.py" in hashes
    assert "notes/proof.md" in hashes and "notes/sources.md" in hashes
    for name, expected in hashes.items():
        assert verify.sha(verify.ROOT/name) == expected


@pytest.mark.parametrize("section,key,value", [
    ("actual_state", "b_T_is_photon_beta_M_or_scalar_gamma", True),
    ("actual_state", "thermal_amplitude", "hbar*pi^2/(30*b_T^4)"),
    ("actual_state", "physical_temperature", "k_B*T_physical=1/b_T"),
    ("named_source_and_solution", "vacuum_H4_anomaly_omitted", True),
    ("named_source_and_solution", "independent_pressure_equation_checked", False),
    ("actual_history_embedding", "all_radiation_to_matter_p_values_realized", True),
    ("endpoint", "endpoint_admitted_as_smooth_metric_or_state_initial_data", True),
    ("endpoint", "A18_future_samplers_or_caps_evaluated_through_endpoint", True),
    ("relationship_to_QSEI_objective", "new_QEI_only_endpoint_proof", True),
])
def test_physical_omissions_and_scope_promotions_fail_full_replay(section, key, value):
    damaged = copy.deepcopy(candidate())
    damaged[section][key] = value
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(damaged, candidate())


@pytest.mark.parametrize("name", ["PARENT_SHA", "PHOTON_SHA"])
def test_old_certificate_pin_cannot_be_bypassed(monkeypatch, name):
    verify.prior_checks.cache_clear()
    try:
        monkeypatch.setattr(verify, name, "0"*64)
        with pytest.raises(ValueError, match="pinned"):
            verify.prior_checks()
    finally:
        verify.prior_checks.cache_clear()
