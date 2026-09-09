"""Independent full ancestry, exact source manifest and every-field mutations."""

import copy
import json

import pytest
from p8_proca_seven_modes import verify

KEYS = [
    "schema",
    "claim",
    "date",
    "status",
    "prior_sha256",
    "source_sha256",
    "formulation",
    "written_proofs",
    "exact_residuals",
    "named_exact_check_count",
    "checked_scalar_entries",
    "proof_checks",
    "actual_fourteen_phase_hamiltonian",
    "actual_cartesian_polarization_projectors",
    "positive_reference_covariances",
    "smooth_reference_preparation",
    "actual_global_product_state",
    "canonical_physical_stress_bridge",
    "controls",
    "verdict",
    "not_established",
    "verification_boundary",
]


def test_readonly_complete_report_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_every_own_source_and_report_field_is_pinned():
    actual = verify.build_report()
    assert actual["source_sha256"] == {
        str(p.relative_to(verify.ROOT)): verify.sha(p) for p in verify.source_files()
    }
    assert set(actual) == set(KEYS)


@pytest.mark.parametrize("key", KEYS)
def test_each_report_field_mutation_rejected(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_or_extra_report_field_rejected(mutation):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    if mutation == "missing":
        del changed["actual_global_product_state"]
    else:
        changed["renormalized_conserved_physical_source_and_UV_completion"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, actual)


def test_report_exact_counts_and_scope():
    actual = verify.build_report()
    assert actual["named_exact_check_count"] == 69
    assert actual["checked_scalar_entries"] == 1360
    assert len(actual["proof_checks"]) == 22
    assert actual["controls"]["rejected_inputs"] == 16
    assert len(actual["source_sha256"]) == 21
    assert actual["claim"] == "P8-S6.103.ACTUAL_GLOBAL_SEVEN_MODE_STATE"
