import copy
import json

import pytest
from p8_s5 import verify


def test_previous_linear_source_pins_unchanged():
    assert len(verify.check_prior_sources()) == 64


def test_saved_report_replays_and_retains_open_control_gate():
    actual = verify.build_report()
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual)
    assert "CONTROL_OPEN" in actual["status"]
    assert "physical 2-to-2 amplitudes" in actual["not_established"]


def test_mutated_report_fails_closed():
    actual = verify.build_report()
    wrong = copy.deepcopy(actual)
    wrong["lapse_hessian"]["at_bounce"] = "0"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(wrong, actual)


def test_nonzero_residual_fails_closed():
    with pytest.raises(ValueError, match="Nonzero"):
        verify.require_zero({"negative_control": 1})
