"""Read-only ancestry, exact source manifest and report mutation rejection."""
import copy
import json

import pytest
from p8_aligned_elimination import verify


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_local_sources_are_pinned():
    assert verify.build_report()["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}


@pytest.mark.parametrize("key", ("claim", "status", "source_sha256", "prior_sha256", "exact_residuals",
                                  "proof_checks", "literal_action", "ordered_elimination", "state_and_boundary",
                                  "full_source_domain", "actual_ODE_remainders", "continuous_constants",
                                  "source_and_local_density_bounds", "preparation", "Fourier_boundary",
                                  "normalized_example", "physical_example", "quantum_boundary", "controls",
                                  "verdict", "not_established", "verification_boundary"))
def test_report_mutation_is_rejected(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_control_count_and_admissibility_boundary():
    assert verify.controls()["rejected_inputs"] == 42
    assert verify.controls()["source_bound_not_an_inverse_or_nonlinear_solution_bound"] is True
