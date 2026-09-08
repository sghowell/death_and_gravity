"""Read-only acoustic reduction report with complete field mutation coverage."""
import copy
import json

import pytest
from p8_vector_liouville import verify

KEYS = ["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","positive_physical_mode_reduction","longitudinal_scalar_and_momentum_potentials","original_clock_coefficients","prepared_source_chart","prepared_potential_variation","continuous_prepared_variation_majorants","fixed_scale_examples","controls","verdict","not_established","verification_boundary"]


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_source_and_report_keys_pinned():
    data = verify.build_report()
    assert data["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}
    assert set(KEYS) == set(data)


@pytest.mark.parametrize("key", KEYS)
def test_each_report_field_mutation_is_rejected(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_or_extra_field_is_rejected(mutation):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    if mutation == "missing":
        del changed["prepared_potential_variation"]
    else:
        changed["full_renormalized_stress_no_loss_norm"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_gate_and_source_counts():
    data = verify.build_report()
    assert data["named_exact_check_count"] == 30
    assert data["checked_scalar_entries"] == 36
    assert len(data["proof_checks"]) == 12
    assert data["controls"]["rejected_inputs"] == 93
    assert len(data["source_sha256"]) == 14
