"""Complete native ancestry, source pinning and all-field mutation controls."""

import copy
import json

import pytest
from p8_polynomial_vacuum import verify

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
    "polynomial_model",
    "exact_finite_regulator_integration",
    "power_counting_counterterm_basis",
    "positive_one_loop_potential",
    "exact_calibrations",
    "controls",
    "verdict",
    "not_established",
    "verification_boundary",
]


def test_complete_read_only_report_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_own_source_manifest_and_complete_fields():
    d = verify.build_report()
    assert d["source_sha256"] == {
        str(p.relative_to(verify.ROOT)): verify.sha(p) for p in verify.source_files()
    }
    assert set(d) == set(KEYS)


@pytest.mark.parametrize("key", KEYS)
def test_every_field_mutation_rejected(key):
    d = verify.build_report()
    changed = copy.deepcopy(d)
    changed[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, d)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_or_extra_field_rejected(mutation):
    d = verify.build_report()
    changed = copy.deepcopy(d)
    if mutation == "missing":
        del changed["positive_one_loop_potential"]
    else:
        changed["full_V_G_B_and_original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, d)


def test_exact_counts_and_original_scope():
    d = verify.build_report()
    assert d["claim"] == "P8-S6.110.POLYNOMIAL_VACUUM_QUANTUM_PREPARATION"
    assert d["named_exact_check_count"] == 37
    assert d["checked_scalar_entries"] == 51
    assert len(d["proof_checks"]) == 39
    assert d["controls"]["rejected_inputs"] == 87
    assert len(d["source_sha256"]) == 18
