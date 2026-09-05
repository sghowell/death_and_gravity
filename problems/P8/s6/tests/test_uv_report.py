import copy
import json

import pytest
from p8_uv import verify


def test_uv_report_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_changed_vacuum_coefficient_and_verdict_rejected():
    actual = verify.build_report()
    bad = copy.deepcopy(actual)
    bad["vacuum_amplitude"]["b2"] = "0"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, actual)
    bad = copy.deepcopy(actual)
    bad["status"] = "UV_COMPLETE"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, actual)
