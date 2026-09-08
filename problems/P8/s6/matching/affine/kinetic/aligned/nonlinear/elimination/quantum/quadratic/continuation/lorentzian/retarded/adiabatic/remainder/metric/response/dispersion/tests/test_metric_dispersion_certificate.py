"""Read-only report replay, source pins and every report-field mutation."""
import copy
import json

import pytest
from p8_vector_metric_dispersion import verify

KEYS = ["schema","claim","date","status","prior_sha256","source_sha256","formulation","written_proofs","exact_residuals","named_exact_check_count","checked_scalar_entries","proof_checks","physical_flat_pair_vertices","physical_flat_contacts","spectral_matrix","radial_moments","thrice_subtracted_normalized_remainder","finite_large_frequency_constant","frozen_pole_and_finite_fourth_symbols","isolated_positive_real_axis_bounds","regular_chart","one_prepared_tree_response_bounds","controls","verdict","not_established","verification_boundary"]


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_sources_and_keys_pinned():
    data = verify.build_report()
    assert data["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}
    assert set(KEYS) == set(data)


@pytest.mark.parametrize("key", KEYS)
def test_report_mutation_rejected(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_or_extra_field_rejected(mutation):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    if mutation == "missing":
        del changed["isolated_positive_real_axis_bounds"]
    else:
        changed["full_causal_matrix_inverse_proved"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_and_source_counts():
    data = verify.build_report()
    assert data["named_exact_check_count"] == 56
    assert data["checked_scalar_entries"] == 75
    assert len(data["proof_checks"]) == 14
    assert data["controls"]["rejected_inputs"] == 52
    assert len(data["source_sha256"]) == 14
