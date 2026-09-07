"""Read-only report replay, actual-input bridges and promotion controls."""
import copy
import json

import pytest
import sympy as sp
from p8_variable_reduction import bridges, dictionary, stationary, verify


def test_frozen_report_rebuilds_with_ordinary_ancestry():
    expected = json.loads(verify.REPORT.read_text())
    verify.validate_report(expected, verify.build_report())


def test_every_source_and_scientific_audit_is_pinned():
    report = verify.build_report()
    sources = {str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}
    assert report["source_sha256"] == sources
    assert len(sources) == 18
    for name in ("test_reduction_independent_audit.py", "test_reduction_center_independent_audit.py"):
        assert "tests/"+name in sources


def test_exact_counts_and_continuous_ancestry_bridges():
    report = verify.build_report()
    assert [len(group) for group in report["exact_residuals"].values()] == [56, 27, 29, 26]
    assert report["exact_residual_count"] == 138
    assert set(bridges.checks().values()) == {0}
    assert bridges.target_input()["source_hash_count"] == 12
    assert all(value > 0 for value in bridges.domain_margins().values())


@pytest.mark.parametrize("key", ("status", "prior_sha256", "source_sha256", "exact_residuals",
                                  "stationary_action_calibration", "actual_background_center_jets",
                                  "continuous_fourth_jet_margin_polynomial", "not_established"))
def test_changed_evidence_is_not_admitted(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "changed"
    with pytest.raises(ValueError):
        verify.validate_report(changed, actual)


def test_scope_and_actual_component_derivative_controls():
    report = verify.build_report()
    assert "ORIGINAL_S6_P8_OPEN" in report["status"]
    assert report["controls"]["rejected_inputs"] == 26
    center = dictionary.center()
    assert not center["controlled_reduction_claimed"]
    assert not center["matter_background_retuned"]
    jets = stationary.center_jets()
    assert "c=1" not in jets["domain"]
    assert jets["metric00_remainder_u0_u2_u4"][0] == (stationary.C-2)**2
    assert jets["remainder00_varphi4_limit"] > 18


@pytest.mark.parametrize("bad", (0.1, sp.Float("0.1"), sp.oo, sp.nan, sp.I))
def test_report_rejects_rounded_or_nonreal_data(bad):
    with pytest.raises(TypeError):
        verify.serialize(bad)


def test_symbolic_algebra_is_not_a_numeric_domain_admission():
    x = sp.Symbol("x", real=True)
    assert verify.serialize(x**2) == "x**2"
    with pytest.raises(TypeError):
        dictionary.center(x=x)
