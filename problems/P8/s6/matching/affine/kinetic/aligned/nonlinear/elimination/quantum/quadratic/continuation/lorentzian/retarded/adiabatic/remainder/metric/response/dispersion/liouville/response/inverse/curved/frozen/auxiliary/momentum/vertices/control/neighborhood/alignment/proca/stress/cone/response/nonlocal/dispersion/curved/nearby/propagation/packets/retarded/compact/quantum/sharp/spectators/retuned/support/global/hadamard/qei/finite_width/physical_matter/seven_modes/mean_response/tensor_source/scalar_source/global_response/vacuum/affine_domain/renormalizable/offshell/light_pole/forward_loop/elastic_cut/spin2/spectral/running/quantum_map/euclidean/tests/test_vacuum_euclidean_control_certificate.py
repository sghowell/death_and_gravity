"""Complete read-only ancestry and all report-field mutation tests."""

import copy
import json

import pytest
from p8_vacuum_euclidean_control import verify

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
    "actual_global_Euclidean_kernel",
    "parent_light_and_stress_normalization",
    "restricted_known_kernel_insertion_tails",
    "actual_composite_coordinate_domain",
    "actual_rational_global_and_pointwise_bounds",
    "controls",
    "verdict",
    "not_established",
    "verification_boundary",
]


def test_complete_read_only_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_own_source_manifest_and_all_report_fields():
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
def test_missing_or_extra_report_field_rejected(mutation):
    d = verify.build_report()
    changed = copy.deepcopy(d)
    if mutation == "missing":
        del changed["actual_composite_coordinate_domain"]
    else:
        changed["full_V_G_B_and_original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, d)


def test_exact_counts_and_original_scope():
    d = verify.build_report()
    assert d["claim"] == "P8-S6.119.POLYNOMIAL_VACUUM_ONE_LOOP_EUCLIDEAN_CONTROL"
    assert d["named_exact_check_count"] == d["checked_scalar_entries"] == 49
    assert len(d["proof_checks"]) == 34
    assert d["controls"]["rejected_inputs"] == 35
    assert len(d["source_sha256"]) == 17
