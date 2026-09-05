import copy
import json

import pytest
from p8_physical import verify


def test_full_exact_report_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_altered_contact_coefficient_is_rejected():
    actual = verify.build_report()
    bad = copy.deepcopy(actual)
    bad["canonically_normalized_velocity_examples"][2]["rational_kernel"] = "0"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, actual)


def test_altered_source_or_scope_is_rejected():
    actual = verify.build_report()
    for key, value in (("source_sha256", {}), ("not_established", [])):
        bad = copy.deepcopy(actual)
        bad[key] = value
        with pytest.raises(ValueError, match="differs"):
            verify.validate_report(bad, actual)
