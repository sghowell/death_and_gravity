import copy
import json
from pathlib import Path

import pytest
import sympy as sp
from p8_variable_vacuum import verify


@pytest.fixture(scope="module")
def candidate():
    return verify.build_report()


def test_actual_recursive_prior_and_independent_counts(candidate):
    assert verify.sha(verify.prior.REPORT) == verify.PRIOR_SHA
    assert verify.sha(verify.P8/"s6"/"FORMULATION.md") == verify.CONTRACT_SHA
    assert candidate["exact_symbolic_residual_count"] == 23
    assert candidate["independent_replay"]["total_exact_scalar_comparisons"] == 171
    assert candidate["checked_controls"]["invalid_calls_rejected"] == 17


def test_all_local_sources_and_reserved_audit_are_inventoried(candidate):
    paths = sorted(verify.ROOT.glob("src/p8_variable_vacuum/*.py"))+sorted(verify.ROOT.glob("tests/*.py"))
    paths += sorted(verify.ROOT.glob("*.md"))+sorted(verify.ROOT.glob("notes/*.md"))
    actual = {str(path.relative_to(verify.ROOT)): verify.sha(path) for path in paths}
    assert actual == candidate["source_sha256"]
    assert "tests/test_local_vacuum_independent_audit.py" in actual
    assert all(".." not in Path(name).parts for name in actual)


def test_frozen_report_replays_when_present(candidate):
    if not verify.REPORT.is_file():
        pytest.skip("No report before the complete root science/source audit")
    verify.validate_report(json.loads(verify.REPORT.read_text()), candidate)


@pytest.mark.parametrize("section", [
    "status", "actual_action_and_domain", "full_stationary_equations", "continuous_proof",
    "prior_context_sha256", "source_sha256", "adopted_V_B_interpretation", "not_established",
])
def test_corrupted_claim_or_data_rejected(candidate, section):
    changed = copy.deepcopy(candidate)
    changed[section] = "A changed hypothesis/result must not replay"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, candidate)


def test_added_global_exclusion_is_rejected(candidate):
    changed = copy.deepcopy(candidate)
    changed["no_off_interval_vacuum_or_UV_completion"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, candidate)


def test_exact_serialization_rejects_floats_and_nonfinite():
    assert verify.serialize(sp.Rational(2, 3)) == "2/3"
    for value in (0.1, sp.Float("0.1"), sp.oo, -sp.oo, sp.zoo, sp.nan):
        with pytest.raises(TypeError):
            verify.serialize(value)


def test_candidate_rebuild_is_read_only(candidate, monkeypatch):
    def forbidden(*_args, **_kwargs):
        raise AssertionError("Read-only replay attempted to write a file")
    monkeypatch.setattr(Path, "write_text", forbidden)
    monkeypatch.setattr(Path, "write_bytes", forbidden)
    verify.build_report.cache_clear()
    verify.validate_report(candidate, verify.build_report())
