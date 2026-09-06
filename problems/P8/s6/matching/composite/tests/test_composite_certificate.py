import copy
import json

import pytest
from p8_composite import verify


def test_certificate_and_pinned_old_scope_replay():
    actual = verify.build_report()
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual)
    assert actual["prior_S6_5_separate_matter_sha256"] == verify.PRIOR_SHA
    assert actual["claim"] == "P8-S6.6.COMPOSITE"


@pytest.mark.parametrize("key", ["prior_scope", "undivided_Bianchi", "local_CD", "not_established"])
def test_changed_scope_or_branch_claim_rejected(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "full healthy CD in the old matter frame"
    with pytest.raises(ValueError, match="differs from exact replay"):
        verify.validate_report(changed, actual)
