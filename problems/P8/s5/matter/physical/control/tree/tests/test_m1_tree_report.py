import copy
import json

import pytest
from p8_m1_tree import verify


def test_read_only_exact_tree_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_tree_scope_and_missing_mode_channel_mutations_rejected():
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed["not_established"] = []
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)
    changed = copy.deepcopy(actual)
    changed["operational_contract"]["quartic_order"] = "one H4 only"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_changed_free_control_prior_rejected(monkeypatch):
    original = verify.sha
    monkeypatch.setattr(verify, "sha", lambda path: "changed" if path == verify.prior.REPORT else original(path))
    with pytest.raises(ValueError, match="certificate changed"):
        verify.prior_checks()
