import copy
import json

import pytest
from p8_composite_response import verify


def test_frozen_source_hashed_report_replays_exactly():
    actual = verify.build_report()
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual)
    assert actual["prior_S6_10_composite_asymptotics_sha256"] == verify.PRIOR_SHA
    assert "tests/test_composite_response_covariant_audit.py" in actual["source_sha256"]


def test_changed_parameter_or_source_contract_does_not_replay():
    actual = verify.build_report()
    modified = copy.deepcopy(actual)
    modified["domain"]["parameter"]["inclusive_upper"] = "1/10000"
    with pytest.raises(ValueError):
        verify.validate_report(modified, actual)
    modified = copy.deepcopy(actual)
    modified["physical_probe"]["literal_u_action_source"] = "M² epsilon sigma p/4"
    with pytest.raises(ValueError):
        verify.validate_report(modified, actual)
