"""Full source-hashed read-only replay and independent tamper rejection."""
import copy
import json

import pytest
from p8_affine_traces import verify


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_source_manifest_covers_every_local_proof_and_test():
    actual = verify.build_report()
    assert actual["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}


@pytest.mark.parametrize("key", ("claim", "status", "source_sha256", "prior_sha256",
                                  "exact_residuals", "proof_checks", "controls", "family",
                                  "geometry", "rank_one", "rank_two", "interfaces_and_units",
                                  "verdict", "not_established", "verification_boundary"))
def test_changed_report_field_is_rejected(key):
    actual = verify.build_report()
    mutated = copy.deepcopy(actual)
    mutated[key] = "not the frozen value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(mutated, actual)


def test_exact_input_and_exceptional_rank_controls():
    assert verify.controls()["rejected_inputs"] == 32
