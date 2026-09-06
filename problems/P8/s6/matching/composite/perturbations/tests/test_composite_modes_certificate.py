import copy
import json

import pytest
from p8_composite_modes import verify


def test_certificate_and_immutable_background_lineage_replay():
    actual = verify.build_report()
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual)
    assert actual["prior_S6_6_composite_background_sha256"] == verify.PRIOR_SHA
    assert actual["claim"] == "P8-S6.7.COMPOSITE"


@pytest.mark.parametrize("key", ["domain", "finite_shift", "vector", "analytic_bounce_screen",
                               "source_normalization_audit", "not_established"])
def test_changed_domain_normalization_or_exclusion_claim_rejected(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "all composite models and ultraviolet completions excluded"
    with pytest.raises(ValueError, match="differs from exact replay"):
        verify.validate_report(changed, actual)
