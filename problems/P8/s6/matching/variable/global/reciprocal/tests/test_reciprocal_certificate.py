import copy
import json
from pathlib import Path

import pytest
from p8_reciprocal_geometry import verify


@pytest.fixture(scope="module")
def candidate():
    return verify.build_report()


def test_actual_ancestor_and_continuous_proof_inputs(candidate):
    assert verify.sha(verify.prior.REPORT) == verify.PRIOR_SHA
    assert candidate["exact_residual_count"] == 8
    assert candidate["checked_controls"]["invalid_exact_bound_inputs_rejected"] == 16
    assert set(candidate["exact_source_and_differential_residuals"].values()) == {"0"}


def test_complete_source_manifest_includes_separately_authored_audit(candidate):
    sources = sorted(verify.ROOT.glob("src/p8_reciprocal_geometry/*.py"))
    sources += sorted(verify.ROOT.glob("tests/*.py"))
    sources += sorted(verify.ROOT.glob("*.md"))+sorted(verify.ROOT.glob("notes/*.md"))
    expected = {str(path.relative_to(verify.ROOT)): verify.sha(path) for path in sources}
    assert candidate["source_sha256"] == expected
    assert "tests/test_reciprocal_independent_audit.py" in expected
    assert all(".." not in Path(path).parts for path in expected)


def test_frozen_report_matches_fresh_candidate(candidate):
    if not verify.REPORT.is_file():
        pytest.skip("The report is absent until complete scientific review")
    verify.validate_report(json.loads(verify.REPORT.read_text()), candidate)


@pytest.mark.parametrize("section", [
    "status", "hypotheses", "continuous_proof_chain", "prior_context_sha256",
    "source_sha256", "not_established", "checked_controls",
])
def test_changed_hypothesis_or_conclusion_is_rejected(candidate, section):
    altered = copy.deepcopy(candidate)
    altered[section] = "A broader unsupported claim is not this certificate"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(altered, candidate)


def test_extra_physical_g_or_uv_exclusion_is_rejected(candidate):
    altered = copy.deepcopy(candidate)
    altered["all_physical_g_bounces_excluded"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(altered, candidate)


def test_replay_is_read_only(candidate, monkeypatch):
    def forbidden(*_args, **_kwargs):
        raise AssertionError("Read-only certificate replay attempted a write")
    monkeypatch.setattr(Path, "write_text", forbidden)
    monkeypatch.setattr(Path, "write_bytes", forbidden)
    verify.build_report.cache_clear()
    verify.validate_report(candidate, verify.build_report())
