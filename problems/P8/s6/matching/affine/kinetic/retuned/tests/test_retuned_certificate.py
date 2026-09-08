"""Strict read-only ancestry, manifest and tamper replay."""
import copy
import json

import pytest
from p8_affine_retuned import verify


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_local_sources_are_pinned():
    assert verify.build_report()["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}


@pytest.mark.parametrize("key", ("claim", "status", "source_sha256", "prior_sha256", "exact_residuals",
                                  "proof_checks", "literal_new_operator", "matching", "primary_constraint",
                                  "positive_principal_range", "characteristic", "units", "controls",
                                  "verdict", "not_established", "verification_boundary"))
def test_report_mutation_is_rejected(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_control_count_and_original_admissibility_boundary():
    assert verify.controls()["rejected_inputs"] == 34
    assert verify.controls()["matter_cone_test_not_replaced_by_positive_energy"] is True
