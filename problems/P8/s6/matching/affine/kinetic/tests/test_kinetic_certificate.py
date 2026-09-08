"""Read-only reproducibility, source manifest and tamper rejection."""
import copy
import json

import pytest
from p8_affine_kinetic import verify


def test_readonly_certificate_replay():
    expected = json.loads(verify.REPORT.read_text())
    verify.validate_report(expected, verify.build_report())


def test_report_includes_every_local_source():
    actual = verify.build_report()
    assert actual["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}


@pytest.mark.parametrize("key", ("claim", "status", "source_sha256", "prior_sha256",
                                  "exact_residuals", "proof_checks", "controls",
                                  "homothetic_verdict", "trace_form_verdict",
                                  "quotient_vector", "coupled_scalar", "nonunit_and_parent_interfaces",
                                  "not_established", "verification_boundary"))
def test_report_mutation_is_rejected(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "not the frozen value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_rejected_input_controls():
    assert verify.controls()["rejected_inputs"] == 24
