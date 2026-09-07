import copy
import json

import pytest
import sympy as sp
from p8_variable_constraints import verify


@pytest.fixture(scope="module")
def candidate():
    return verify.build_report()


def test_candidate_counts_and_scoped_status(candidate):
    assert candidate["claim"] == "P8-S6.22.CONSTRAINTS"
    assert "SCALAR_HEALTH_AND_ORIGINAL_P8_OPEN" in candidate["status"]
    assert candidate["exact_residual_count"] == 440
    assert candidate["primary_independent_output_bridges"]["total_exact_output_comparisons"] == 1264
    assert candidate["checked_omissions_and_invalid_inputs"]["invalid_or_inexact_calls_rejected_after_cache_warmup"] == 18


def test_manifest_pins_only_this_child_and_actual_ancestor(candidate):
    expected = sorted(verify.ROOT.glob("src/p8_variable_constraints/*.py"))
    expected += sorted(verify.ROOT.glob("tests/*.py"))
    expected += sorted(verify.ROOT.glob("*.md"))+sorted(verify.ROOT.glob("notes/*.md"))
    assert candidate["source_sha256"] == {str(path.relative_to(verify.ROOT)): verify.sha(path) for path in expected}
    assert all(not path.startswith(("../", "response/", "global/")) for path in candidate["source_sha256"])
    assert verify.sha(verify.prior.REPORT) == verify.PRIOR_SHA


@pytest.mark.parametrize("section,key,replacement", [
    ("regular_scalar_reduction", "strict_rank", "scalar health proved"),
    ("complete_vector_result", "physical_principal_speed", "all modes subluminal"),
    ("nonuniform_scalar_guard", "inference_rejected", "center polynomial is a cone"),
])
def test_corrupted_scope_or_formula_is_rejected(candidate, section, key, replacement):
    altered = copy.deepcopy(candidate)
    altered[section][key] = replacement
    with pytest.raises(ValueError):
        verify.validate_report(altered, candidate)


def test_build_is_read_only_after_cache_clear(candidate):
    paths = [verify.ROOT/path for path in candidate["source_sha256"]]
    paths += [verify.prior.REPORT]
    before = {path: verify.sha(path) for path in paths}
    verify.build_report.cache_clear()
    verify.validate_report(candidate, verify.build_report())
    assert {path: verify.sha(path) for path in paths} == before


def test_report_serializer_rejects_inexact_values():
    for value in (0.1, sp.Float("0.1"), sp.oo, sp.nan):
        with pytest.raises(TypeError):
            verify.serialize(value)


def test_frozen_report_replays_exactly(candidate):
    expected = json.loads(verify.REPORT.read_text())
    verify.validate_report(expected, candidate)
