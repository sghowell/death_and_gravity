"""Read-only report replay and mutation/manifest rejection."""
import copy
import json

import pytest
from p8_vector_quadratic import verify

KEYS = ("schema", "claim", "date", "status", "prior_sha256", "source_sha256",
        "formulation", "written_proofs", "exact_residuals", "named_exact_check_count",
        "checked_scalar_entries", "proof_checks", "conventions", "accepted_pole_bilinears",
        "independent_scalar_conformal_check", "auxiliary_metric_definition",
        "independent_coordinate_geometry", "local_auxiliary_metric_quadratic_density",
        "withheld_literal_reference", "actual_mass_profiles", "actual_clock_lapse_local_coefficients",
        "matching_boundary", "controls", "verdict", "not_established", "verification_boundary")


def test_readonly_certificate_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_all_local_sources_and_report_keys_are_pinned():
    report = verify.build_report()
    assert report["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}
    assert set(KEYS) == set(report)


@pytest.mark.parametrize("key", KEYS)
def test_report_mutation_is_rejected(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_or_extra_report_fields_are_rejected(mutation):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    if mutation == "missing":
        del changed["withheld_literal_reference"]
    else:
        changed["finite_retarded_response_completed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_identity_and_rejection_counts():
    data = verify.build_report()
    assert data["named_exact_check_count"] == 63
    assert data["checked_scalar_entries"] == 153
    assert len(data["proof_checks"]) == 16
    assert data["controls"]["rejected_inputs"] == 170
    assert len(data["source_sha256"]) == 17
