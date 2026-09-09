"""Read-only parent replay, exact manifest and every report-key mutation."""

import copy
import json
from functools import cache

import pytest
from p8a_nonminimal import verify


@cache
def candidate():
    return verify.build_report()


def test_complete_frozen_report_rebuilds_without_writing():
    before = verify.REPORT.read_bytes()
    verify.validate_report(json.loads(before), candidate())
    assert verify.REPORT.read_bytes() == before


def test_complete_source_and_scientific_inventory():
    report = candidate()
    assert set(report) == set(verify.REPORT_KEYS)
    assert len(report["source_sha256"]) == 20
    assert report["named_exact_check_count"] == 57
    assert report["checked_scalar_entries"] == 72
    assert len(report["proof_checks"]) == 23
    assert report["controls"]["rejected_inputs"] == 55
    for name, expected in report["source_sha256"].items():
        assert verify.sha(verify.ROOT / name) == expected


@pytest.mark.parametrize("key", verify.REPORT_KEYS)
@pytest.mark.parametrize("mode", ["omit", "replace"])
def test_every_report_field_mutation_fails(key, mode):
    damaged = copy.deepcopy(candidate())
    if mode == "omit":
        del damaged[key]
    else:
        damaged[key] = {"unapproved_replacement": True}
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(damaged, candidate())


def test_unapproved_extra_report_claim_fails():
    damaged = copy.deepcopy(candidate())
    damaged["original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(damaged, candidate())


def test_actual_source_and_scope_promotions_are_absent():
    report = candidate()
    assert not report["actual_thermal_history"][
        "future_state_cap_verified_on_every_shorter_segment"
    ]
    assert not report["coherent_obstruction"][
        "actual_self_consistent_SEE_counterexample"
    ]
    assert not report["conformal_source"]["reference_trace_anomaly_omitted"]
    assert not report["cosmological_theorem"][
        "future_caps_derived_from_past_or_effective_Newton_positivity"
    ]
    assert (
        report["prior_sha256"]["A19_fully_rebuilt_with_A18_A17_A16_ancestry"]
        == verify.PARENT_SHA
    )


def test_changed_parent_pin_is_rejected(monkeypatch):
    verify.prior_checks.cache_clear()
    try:
        monkeypatch.setattr(verify, "PARENT_SHA", "0" * 64)
        with pytest.raises(ValueError, match="pinned"):
            verify.prior_checks()
    finally:
        verify.prior_checks.cache_clear()
