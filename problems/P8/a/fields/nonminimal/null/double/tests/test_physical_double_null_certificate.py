"""Read-only replay, frozen source bytes and complete mutation controls."""

import copy
import json
from functools import cache

import pytest
from p8a_double_null import verify


@cache
def candidate():
    return verify.build_report()


def test_complete_report_and_source_bytes_replay_without_writes():
    before = verify.REPORT.read_bytes()
    verify.validate_report(json.loads(before), candidate())
    assert verify.REPORT.read_bytes() == before
    hashes = candidate()["source_sha256"]
    assert len(hashes) == 17
    for path, expected in hashes.items():
        assert verify.sha(verify.ROOT / path) == expected


def test_science_inventory_and_precise_scope():
    report = candidate()
    assert set(report) == set(verify.REPORT_KEYS)
    assert report["named_exact_check_count"] == 39
    assert report["checked_scalar_entries"] == 51
    assert len(report["proof_checks"]) == 18
    assert report["controls"]["rejected_inputs"] == 39
    assert not report["prior_sha256"]["single_null_cap_bound_inferred_from_parent"]
    assert not report["normalized_product_bound"][
        "optimal_sampler_or_optimal_physical_QEI_claimed"
    ]
    assert not report["curved_conformal_transport"][
        "single_ray_focusing_or_SEE_incompleteness_deduced"
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
