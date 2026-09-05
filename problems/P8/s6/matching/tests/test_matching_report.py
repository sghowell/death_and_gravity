import copy
import json

import pytest
from p8_matching import verify


def test_matching_report_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_changed_bounds_and_overbroad_verdict_rejected():
    actual = verify.build_report()
    bad = copy.deepcopy(actual)
    bad["nonminimal_parent_error_bounds"]["common_C2_error_necessary_lower_bound"] = "0"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, actual)
    bad = copy.deepcopy(actual)
    bad["exact_quartic_Horndeski_map"]["result"] = "ORIGINAL_D_BOUNCE_IS_SINGULAR"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, actual)
    bad = copy.deepcopy(actual)
    bad["status"] = "ALL_D_UV_COMPLETIONS_EXCLUDED"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, actual)
