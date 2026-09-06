import copy
import json

import pytest
from p8_m1_weyl_scalar import verify


def test_checked_certificate_replays_full_pinned_lineage():
    expected = json.loads(verify.REPORT.read_text())
    verify.validate_report(expected, verify.build_report())


def test_claim_and_scope_changes_fail_closed():
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed["operational_contract"]["data_distinctions"] = "all initial data are identical"
    with pytest.raises(ValueError):
        verify.validate_report(changed, actual)
