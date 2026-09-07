"""Read-only candidate/report integrity and exact proof-to-code bridges."""

import copy
import json

import pytest
import sympy as sp
from p8_variable_forced import verify


@pytest.fixture(scope="module")
def report():
    return verify.build_report()


def test_all_residuals_and_margins(report):
    assert all(value == "0" for group in report["exact_residuals"].values() for value in group.values())
    assert all(sp.Rational(value) > 0 for group in report["strict_continuous_margins"].values() for value in group.values())


def test_frozen_preparation_and_original_contract_replayed(report):
    assert report["prior_sha256"]["S6_23_prepared_and_replayed_S6_21_S6_20_lineage"] == verify.PREPARED_SHA
    assert report["prior_sha256"]["adopted_original_S6_contract"] == verify.CONTRACT_SHA


def test_every_current_source_is_hashed(report):
    assert report["source_sha256"] == {str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}
    assert "tests/test_forced_independent_audit.py" in report["source_sha256"]


def test_independent_continuous_engine_is_required(report):
    evidence = report["independent_Fraction_replay"]
    assert len(evidence["exact_arithmetic_regression_rows"]) == 3
    assert evidence["continuous_polynomial_replay"]["positive_polynomial_coefficients"] > 100


def test_state_source_domain_and_nonlocal_scope(report):
    assert "zero at -r" in report["unchanged_physical_parent"]["initial_data"]
    assert "0<=K" in report["unchanged_physical_parent"]["domain"]
    assert "retarded" in report["causal_reduction"]["inverses"]
    assert any("No local-in-time" in value for value in report["not_established"])
    assert any("No delta-convergence" in value for value in report["not_established"])
    assert any("No identification" in value for value in report["not_established"])
    assert report["checked_controls"]["rejected_exact_domain_calls"] == 26


@pytest.mark.parametrize("key", ["status", "source_sha256", "causal_reduction", "not_established"])
def test_changed_scope_hash_or_bound_rejected(report, key):
    modified = copy.deepcopy(report)
    modified[key] = "OMITTED"
    with pytest.raises(ValueError):
        verify.validate_report(modified, report)


def test_extra_unreplayed_key_rejected(report):
    modified = copy.deepcopy(report)
    modified["low_energy_UV_completion"] = True
    with pytest.raises(ValueError):
        verify.validate_report(modified, report)


def test_inexact_report_data_rejected():
    for value in (0.1, sp.Float(".1"), sp.oo, sp.nan):
        with pytest.raises(TypeError):
            verify.serialize(value)


def test_pinned_report_read_only(report):
    before = verify.REPORT.read_bytes()
    verify.validate_report(json.loads(before), report)
    assert verify.REPORT.read_bytes() == before
