"""Read-only ancestry, exact source manifest and report-mutation controls."""
import copy
import json

import pytest
from p8_affine_nonlinear import verify


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_local_sources_are_pinned():
    assert verify.build_report()["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}


@pytest.mark.parametrize("key", ("claim", "status", "source_sha256", "prior_sha256", "exact_residuals",
                                  "proof_checks", "literal_action", "complete_ADM_chart",
                                  "joint_trace_Legendre_map", "Maxwell_and_spatial_constraints",
                                  "actual_background_auxiliary_Jacobian", "clock_bounds", "local_constraint_theorem",
                                  "complete_physical_count", "controls", "verdict", "not_established", "verification_boundary"))
def test_report_mutation_is_rejected(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_control_count_and_admissibility_boundary():
    assert verify.controls()["rejected_inputs"] == 7
    assert verify.controls()["local_count_distinguished_from_nonlinear_health_and_X_zero_vacuum"] is True
