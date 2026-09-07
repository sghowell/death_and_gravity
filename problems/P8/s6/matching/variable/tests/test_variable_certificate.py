import copy
import json
from fractions import Fraction

import pytest
import sympy as sp
from p8_variable_beta import verify


def test_frozen_report_replays_without_writing_and_hashes_the_exact_child_inventory():
    before = verify.REPORT.read_bytes()
    expected = json.loads(before)
    actual = verify.build_report()
    verify.validate_report(expected, actual)
    assert verify.REPORT.read_bytes() == before
    inventory = sorted(verify.ROOT.glob("src/p8_variable_beta/*.py"))+sorted(verify.ROOT.glob("tests/*.py"))
    inventory += sorted(verify.ROOT.glob("*.md"))+sorted(verify.ROOT.glob("notes/*.md"))
    assert set(expected["source_sha256"]) == {str(path.relative_to(verify.ROOT)) for path in inventory}
    assert len(inventory) == 16
    assert all(verify.sha(verify.ROOT/name) == digest for name, digest in expected["source_sha256"].items())
    assert "tests/test_variable_independent_audit.py" in expected["source_sha256"]


def test_actual_prior_replay_preserves_constant_tree_scope_and_adopted_contract():
    report = verify.build_report()
    pins = report["prior_context_sha256"]
    assert pins["S6_17_constant_beta_tree_scope_and_replay"] == verify.TREE_SHA
    assert pins["adopted_S6_matching_contract"] == verify.CONTRACT_SHA
    assert len(pins) == 4
    assert "S6.19" not in json.dumps(pins)
    assert sum(len(group) for group in report["exact_residuals"].values()) == 64
    assert report["status"].endswith("ORIGINAL_P8_OPEN")
    assert "not the original DHOST" in report["actual_local_solution"]["original_target_boundary"]
    assert "No equal retarded data" in report["shrinking_inner_window_theorem"]["forcing_control"]


def test_independent_actual_and_differentiated_primary_bridges_are_not_preset_verdicts():
    bridges = verify.independent_bridges()
    profiles, jets = bridges["actual_background_profiles"], bridges["canonical_second_jets"]
    assert len(profiles) == 7 and len(jets) == 5
    assert sum(len(case["primary_minus_independent"]) for case in profiles) == 98
    assert sum(len(case["primary_minus_independent"]) for case in jets) == 50
    assert all(value == "0" for group in bridges.values() for case in group
               for value in case["primary_minus_independent"].values())


def test_corrupted_identity_or_claim_cannot_pass_exact_read_only_comparison():
    actual = verify.build_report()
    wrong = copy.deepcopy(actual)
    wrong["exact_residuals"]["canonical_action_boundary_and_adjoint"]["light_equation_keeps_adjoint_connection"] = "1"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(wrong, actual)
    wrong = copy.deepcopy(actual)
    wrong["shrinking_inner_window_theorem"]["window_guard"] = "UNIFORM_FIXED_PHYSICAL_WINDOW"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(wrong, actual)


def test_exact_report_round_trip_and_omission_controls_preserve_domains():
    report = verify.build_report()
    assert json.loads(json.dumps(report)) == report
    assert verify.serialize({"q": Fraction(1, 3)}) == {"q": "1/3"}
    for invalid in (0.1, sp.Float("0.1"), sp.oo, sp.nan):
        with pytest.raises(TypeError):
            verify.serialize(invalid)
    controls = report["checked_controls"]
    assert controls["invalid_exact_domain_calls_rejected"] == 20
    assert controls["limiting_operator_omission_control"]["omitted_derivative_residual"] == "1/5"
    assert controls["limiting_operator_omission_control"]["particular_to_algebraic_ratio"] == "5/6"
    assert controls["source_aware_background_controls"]["omitted_clock_exchange_at_c1_u1over10"] != "0"
