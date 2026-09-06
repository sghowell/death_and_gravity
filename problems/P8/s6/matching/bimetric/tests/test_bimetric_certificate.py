import copy
import json

import pytest
from p8_bimetric import verify


def test_certificate_replays_both_pinned_inputs_and_new_checks():
    actual = verify.build_report()
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual)
    assert actual["prior_certificate_sha256"] == {"S6.2": verify.PARENT_SHA, "S5.11.CD": verify.SCALAR_SHA}
    assert actual["status"].endswith("GENERAL_S6_MATCHING_OPEN")


def test_mutated_bounce_conclusion_rejected():
    actual = verify.build_report()
    altered = copy.deepcopy(actual)
    altered["conclusion"] = "CD bounce matched"
    with pytest.raises(ValueError, match="differs from exact replay"):
        verify.validate_report(altered, actual)


def test_mutated_source_trace_or_pole_rejected():
    actual = verify.build_report()
    altered = copy.deepcopy(actual)
    altered["flat_conserved_source_response"] = "massive exchange uses the massless trace projector"
    with pytest.raises(ValueError, match="differs from exact replay"):
        verify.validate_report(altered, actual)


def test_exact_positive_expressions_and_controls():
    assert verify.positive_checks()
    assert verify.control_checks()
    assert verify.residuals()
