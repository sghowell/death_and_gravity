import copy
import json

import pytest
from p8_m1_physical import verify


def test_pinned_report_exact_replay():
    report = verify.build_report()
    verify.validate_report(json.loads(verify.REPORT.read_text()), report)
    assert report["claim"] == "P8-S5.6.CD"
    assert len(report["phase_kernel_examples"]) == 26
    assert len(report["unnormalized_velocity_examples"]) == 10
    assert all(not row["normalized"] for row in report["unnormalized_velocity_examples"])


def test_changed_report_and_wrong_scope_rejected():
    report = verify.build_report()
    for key, wrong in (("claim", "P8_COMPLETE"), ("status", "M1_CONTROL_PROVED")):
        tampered = copy.deepcopy(report)
        tampered[key] = wrong
        with pytest.raises(ValueError, match="differs"):
            verify.validate_report(tampered, report)
