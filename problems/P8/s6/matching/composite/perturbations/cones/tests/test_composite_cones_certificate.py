import copy
import json

import pytest
from p8_composite_cones import verify


def test_certificate_and_full_immutable_lineage_replay():
    actual = verify.build_report()
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual)
    assert actual["prior_S6_7_composite_modes_sha256"] == verify.PRIOR_SHA
    assert actual["claim"] == "P8-S6.8.COMPOSITE"


@pytest.mark.parametrize("key", ["domain", "implication_chain", "tensor_principal_cones",
                               "actual_zero_vector_outside_contract_control", "not_established"])
def test_inflated_scope_or_changed_contract_fails_replay(key):
    actual = verify.build_report()
    altered = copy.deepcopy(actual)
    altered[key] = "all ultraviolet completions excluded"
    with pytest.raises(ValueError, match="differs from exact replay"):
        verify.validate_report(altered, actual)
