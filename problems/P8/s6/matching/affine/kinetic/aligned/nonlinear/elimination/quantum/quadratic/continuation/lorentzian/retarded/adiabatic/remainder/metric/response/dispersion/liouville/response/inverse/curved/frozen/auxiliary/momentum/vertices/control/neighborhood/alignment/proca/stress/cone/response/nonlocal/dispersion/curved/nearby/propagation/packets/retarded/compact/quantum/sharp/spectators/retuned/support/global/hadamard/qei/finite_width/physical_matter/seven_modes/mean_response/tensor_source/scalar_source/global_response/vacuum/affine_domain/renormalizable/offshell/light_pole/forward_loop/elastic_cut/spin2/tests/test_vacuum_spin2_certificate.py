"""Complete ancestry, own source manifest and all report-field mutations."""

import copy
import json

import pytest
from p8_vacuum_spin2 import verify

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
    "literal_canonical_stress_projection",
    "both_actual_triangles_Ward_identity_slope_and_transfer_domain",
    "literal_graviton_contraction_and_fixed_negative_t_limit",
    "exact_selected_and_transfer_calibrations",
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
        del changed["literal_graviton_contraction_and_fixed_negative_t_limit"]
    else:
        changed["full_V_G_B_and_original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, d)


def test_exact_counts_and_original_scope():
    d = verify.build_report()
    assert d["claim"] == "P8-S6.115.POLYNOMIAL_VACUUM_SPIN_TWO_FORM_FACTOR"
    assert d["named_exact_check_count"] == d["checked_scalar_entries"] == 34
    assert len(d["proof_checks"]) == 28
    assert d["controls"]["rejected_inputs"] == 30
    assert len(d["source_sha256"]) == 16
