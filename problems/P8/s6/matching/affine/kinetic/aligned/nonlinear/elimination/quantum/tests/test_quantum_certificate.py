"""Read-only ancestry, strict local manifest and report mutation controls."""
import copy
import json

import pytest
from p8_aligned_quantum import verify


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_local_sources_are_pinned():
    assert verify.build_report()["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}


@pytest.mark.parametrize("key", ("claim", "status", "source_sha256", "prior_sha256", "exact_residuals",
                                  "proof_checks", "literal_action", "complete_constant_coefficient_kernel",
                                  "complete_determinant", "regulator_and_scheme", "pole_and_finite_weights",
                                  "potential_normalization", "actual_clock_coefficient_jets_at_mu_m0",
                                  "continuous_domain_and_bound", "physical_scale_example", "reference_scale_boundary",
                                  "source_and_sign_boundary", "controls", "verdict", "not_established", "verification_boundary"))
def test_report_mutation_is_rejected(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_control_count_and_admissibility_boundary():
    assert verify.controls()["rejected_inputs"] == 26
    assert verify.controls()["finite_MSbar_potential_not_a_UV_or_full_quantum_verdict"] is True
