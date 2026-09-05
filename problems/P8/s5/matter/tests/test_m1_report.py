import copy
import json

import pytest
from p8_m1 import verify


def test_report_matches_read_only_replay():
    expected = json.loads(verify.REPORT.read_text())
    verify.validate_report(expected, verify.build_report())


def test_prior_sources_and_certificates_are_pinned():
    assert verify.prior_checks() == verify.PRIOR


def test_tampered_result_is_rejected():
    report = verify.build_report()
    tampered = copy.deepcopy(report)
    tampered["negative_control"]["at_bounce"] = "0"
    with pytest.raises(ValueError):
        verify.validate_report(tampered, report)


def test_no_physical_or_loop_completion_is_claimed():
    report = verify.build_report()
    assert "PHYSICAL_INTERACTION_CONTROL_OPEN" in report["status"]
    assert report["local_units"]["coefficient_count"] == 85
    assert "calculated radiative corrections or technical naturalness" in report["not_established"]
