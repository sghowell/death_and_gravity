"""Full read-only report replay and every field mutation control."""

import copy
import json
from functools import cache

import pytest
from p8a_null_focusing import verify


@cache
def candidate():
    return verify.build_report()


def test_complete_report_and_frozen_source_bytes_replay_without_writes():
    before = verify.REPORT.read_bytes()
    verify.validate_report(json.loads(before), candidate())
    assert verify.REPORT.read_bytes() == before
    hashes = candidate()["source_sha256"]
    assert len(hashes) == 19
    for path, expected in hashes.items():
        assert verify.sha(verify.ROOT / path) == expected


def test_inventory_and_precise_null_theorem_scope():
    report = candidate()
    assert set(report) == set(verify.REPORT_KEYS)
    assert report["named_exact_check_count"] == report["checked_scalar_entries"] == 35
    assert len(report["proof_checks"]) == 28
    assert report["controls"]["rejected_inputs"] == 24
    assert not report["prior_sha256"]["null_line_QEI_or_state_homogeneity_assumed"]
    assert not report["complete_geometry_control"]["allowed_small_source_SEE_witness"]
    assert not report["null_theorem"][
        "future_caps_inferred_from_actual_past_or_effective_Newton_constant"
    ]


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


def test_added_root_closure_claim_fails():
    damaged = copy.deepcopy(candidate())
    damaged["P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(damaged, candidate())


def test_parent_pin_cannot_be_bypassed(monkeypatch):
    verify.prior_checks.cache_clear()
    try:
        monkeypatch.setattr(verify, "PARENT_SHA", "0" * 64)
        with pytest.raises(ValueError, match="pinned"):
            verify.prior_checks()
    finally:
        verify.prior_checks.cache_clear()
