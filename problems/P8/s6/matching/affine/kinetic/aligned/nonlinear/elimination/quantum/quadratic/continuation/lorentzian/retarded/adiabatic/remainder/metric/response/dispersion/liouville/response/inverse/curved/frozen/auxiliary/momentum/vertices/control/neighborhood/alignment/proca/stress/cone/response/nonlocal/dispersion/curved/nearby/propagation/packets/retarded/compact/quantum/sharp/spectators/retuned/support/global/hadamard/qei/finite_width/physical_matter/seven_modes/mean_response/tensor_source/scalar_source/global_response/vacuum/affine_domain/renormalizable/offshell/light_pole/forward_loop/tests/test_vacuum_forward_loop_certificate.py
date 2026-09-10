"""Complete ancestry, source manifest and every report-field mutation."""

import copy
import json

import pytest
from p8_vacuum_forward_loop import verify

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
    "actual_full_vertex_normalization",
    "holomorphic_forward_domain",
    "uniform_angular_remainder",
    "subtracted_radial_derivative_bound",
    "fixed_local_counterterm_scheme",
    "exact_error_calibrations",
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
        del changed["exact_error_calibrations"]
    else:
        changed["full_V_G_B_and_original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, d)


def test_exact_counts_and_original_scope():
    d = verify.build_report()
    assert d["claim"] == "P8-S6.113.POLYNOMIAL_VACUUM_ONE_LOOP_FORWARD_COEFFICIENT"
    assert d["named_exact_check_count"] == 54
    assert d["checked_scalar_entries"] == 69
    assert len(d["proof_checks"]) == 44
    assert d["controls"]["rejected_inputs"] == 16
    assert len(d["source_sha256"]) == 19
