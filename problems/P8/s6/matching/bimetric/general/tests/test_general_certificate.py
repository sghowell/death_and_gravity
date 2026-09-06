import copy
import json

import pytest
from p8_bimetric_general import verify


def test_general_certificate_and_beta1_input_replay():
    actual = verify.build_report()
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual)
    assert actual["prior_S6_3_beta1_sha256"] == verify.PRIOR_SHA
    assert actual["claim"] == "P8-S6.4.HR"


@pytest.mark.parametrize("key", ["root_equation", "nonzero_P_weighted_equation", "robust_CD_mismatch"])
def test_changed_branch_or_window_certificate_rejected(key):
    actual = verify.build_report()
    altered = copy.deepcopy(actual)
    altered[key] = "root branch or matter lapse omitted"
    with pytest.raises(ValueError, match="differs from exact replay"):
        verify.validate_report(altered, actual)


def test_exact_positive_and_control_checks():
    assert verify.positive_checks()
    assert verify.control_checks()
    assert verify.residuals()
