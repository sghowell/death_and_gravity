import copy
import json

import pytest
from p8_a26_vacuum import verify


def test_frozen_report_replays_exactly_without_writes():
    expected = json.loads(verify.REPORT.read_text())
    before = verify.sha(verify.REPORT)
    verify.validate_report(expected, verify.build_report())
    assert verify.sha(verify.REPORT) == before


def test_source_manifest_matches_every_new_direct_source():
    report = verify.build_report()
    paths = list(verify.ROOT.glob("src/p8_a26_vacuum/*.py"))+list(verify.ROOT.glob("tests/*.py"))
    paths += list(verify.ROOT.glob("*.md"))+list(verify.ROOT.glob("notes/*.md"))
    assert set(report["source_sha256"]) == {str(path.relative_to(verify.ROOT)) for path in paths}
    for name, value in report["source_sha256"].items():
        assert verify.sha(verify.ROOT/name) == value
    assert "tests/test_a26_independent_audit.py" in report["source_sha256"]


def test_parent_contract_is_pinned_and_replayed():
    pins = verify.prior_checks()
    assert pins["S6_1_smooth_extension_context"] == verify.PRIOR_SHA
    assert pins["adopted_S6_contract"] == verify.CONTRACT_SHA


@pytest.mark.parametrize("key", ["status", "theorem", "source_regularity", "not_established", "source_sha256"])
def test_report_corruption_is_rejected(key):
    actual = verify.build_report()
    broken = copy.deepcopy(actual)
    broken[key] = "corrupted"
    with pytest.raises(ValueError):
        verify.validate_report(broken, actual)


def test_domain_boundaries_are_not_promoted_to_uv_or_instability():
    actual = verify.build_report()
    assert actual["checked_omission_controls"]["invalid_inputs_rejected"] == 8
    assert not actual["checked_omission_controls"]["regularity_is_physical_instability"]
    assert actual["checked_omission_controls"]["r_zero_generic_verdict"] is False
    assert "MODIFIED_MATCHING_AND_P8_B_OPEN" in actual["status"]
    assert len(actual["not_established"]) == 7


def test_separate_fraction_primary_bridges_are_replayed():
    actual = verify.bridges()
    assert actual["Laurent_fixtures_checked"] == 10
    assert actual["direct_quartic_minus_Fraction_log_IBP"] == "0"
