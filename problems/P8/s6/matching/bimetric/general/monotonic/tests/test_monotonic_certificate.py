import copy
import json

import pytest
from p8_bimetric_monotonic import verify


def test_certificate_and_pinned_general_parent_replay():
    actual = verify.build_report()
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual)
    assert actual["prior_S6_4_HR_sha256"] == verify.PRIOR_SHA
    assert actual["claim"] == "P8-S6.5.HR"


@pytest.mark.parametrize("key", ["identically_zero_exception", "boundary_root_argument", "sharp_CD_endpoint_mismatch"])
def test_wrong_case_or_endpoint_claim_rejected(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "discard zero polynomial or root boundary"
    with pytest.raises(ValueError, match="differs from exact replay"):
        verify.validate_report(changed, actual)


def test_exact_positive_domain_and_residual_checks():
    assert verify.positive_checks()
    assert verify.residuals()
