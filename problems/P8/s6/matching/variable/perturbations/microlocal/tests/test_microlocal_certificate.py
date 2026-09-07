import copy
import json
from pathlib import Path

import pytest
import sympy as sp
from p8_variable_microlocal import verify


@pytest.fixture(scope="module")
def candidate():
    return verify.build_report()


def test_actual_recursive_prior_pin_and_exact_counts(candidate):
    assert verify.sha(verify.prior.REPORT) == verify.PRIOR_SHA
    assert candidate["continuous_Bernstein_coefficient_count"] == 611
    assert candidate["independent_replay"]["total_exact_scalar_comparisons"] == 1373
    assert candidate["exact_time_normalization_and_energy"]["order_k_energy_skew_identity_count"] == 36


def test_complete_local_manifest(candidate):
    paths = sorted(verify.ROOT.glob("src/p8_variable_microlocal/*.py"))
    paths += sorted(verify.ROOT.glob("tests/*.py"))+sorted(verify.ROOT.glob("*.md"))+sorted(verify.ROOT.glob("notes/*.md"))
    expected = {str(path.relative_to(verify.ROOT)): verify.sha(path) for path in paths}
    assert candidate["source_sha256"] == expected
    assert "tests/test_microlocal_independent_audit.py" in expected
    assert all(".." not in Path(path).parts for path in expected)


def test_report_replays_when_frozen(candidate):
    if not verify.REPORT.is_file():
        pytest.skip("No report is written before the root's complete proof review")
    verify.validate_report(json.loads(verify.REPORT.read_text()), candidate)


@pytest.mark.parametrize("section", [
    "status", "center_physical_g_principal", "generic_rational_symbol",
    "continuous_positive_principal_forms", "constructive_uniform_Cauchy_inverse",
    "prior_context_sha256", "source_sha256",
])
def test_corrupted_report_rejected(candidate, section):
    corrupted = copy.deepcopy(candidate)
    corrupted[section] = {"corruption": "exact replay must reject this"}
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(corrupted, candidate)


def test_extra_claim_rejected(candidate):
    corrupted = copy.deepcopy(candidate)
    corrupted["full_UV_health"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(corrupted, candidate)


def test_exact_matrix_serialization_and_inexact_rejection():
    assert verify.serialize(sp.Matrix([[sp.Rational(1, 3)]])) == [["1/3"]]
    for value in (0.1, sp.Float("0.1"), sp.oo, sp.Matrix([[sp.Float("0.1")]])):
        with pytest.raises(TypeError):
            verify.serialize(value)


def test_candidate_is_read_only(candidate, monkeypatch):
    def forbidden(*_args, **_kwargs):
        raise AssertionError("Read-only verification attempted a file write")
    monkeypatch.setattr(Path, "write_text", forbidden)
    monkeypatch.setattr(Path, "write_bytes", forbidden)
    verify.build_report.cache_clear()
    verify.validate_report(candidate, verify.build_report())
