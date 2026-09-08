"""Read-only ancestry/source replay and strict report tamper rejection."""
import copy
import json

import pytest
from p8_affine_ricci import verify


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_every_local_source_is_hashed():
    assert verify.build_report()["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}


@pytest.mark.parametrize("key", ("claim", "status", "source_sha256", "prior_sha256",
                                  "exact_residuals", "proof_checks", "controls", "operator",
                                  "curved_reduction", "coupled_scalar_calibration", "units",
                                  "verdict", "not_established", "verification_boundary"))
def test_mutated_report_is_rejected(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_rejected_inputs_and_rank_controls():
    assert verify.controls()["rejected_inputs"] == 28
