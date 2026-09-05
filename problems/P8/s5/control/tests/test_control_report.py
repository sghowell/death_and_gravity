import copy
import json

import pytest
from p8_control import verify


def test_control_report_replays():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_changed_threshold_or_scope_rejected():
    actual = verify.build_report()
    bad = copy.deepcopy(actual)
    bad["all_time_bounds"]["q_threshold"] = "1"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, actual)
    bad = copy.deepcopy(actual)
    bad["not_established"] = []
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, actual)
