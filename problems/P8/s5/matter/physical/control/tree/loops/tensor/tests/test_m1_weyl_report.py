import copy
import json

import pytest
from p8_m1_weyl import verify


def test_certificate_replay_and_pinned_prior():
    actual = verify.build_report()
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual)
    assert actual["prior_S5_9_CD_sha256"] == verify.PRIOR_SHA
    assert "tests/test_m1_weyl_action_audit.py" in actual["source_sha256"]
    assert "QUANTUM_CAUSALITY_OPEN" in actual["status"]


def test_claim_scope_mutation_is_rejected():
    original = {"exact_fourth_order_branch": False, "finite_coefficient_input": True}
    mutant = copy.deepcopy(original)
    mutant["exact_fourth_order_branch"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(original, mutant)


def test_all_symbolic_controls_are_replayed():
    assert all(value == 0 for group in verify.residuals().values() for value in group.values())
    assert verify.control_checks()["naive_y_speed_used_as_physical_at_bounce_error"] == 64
