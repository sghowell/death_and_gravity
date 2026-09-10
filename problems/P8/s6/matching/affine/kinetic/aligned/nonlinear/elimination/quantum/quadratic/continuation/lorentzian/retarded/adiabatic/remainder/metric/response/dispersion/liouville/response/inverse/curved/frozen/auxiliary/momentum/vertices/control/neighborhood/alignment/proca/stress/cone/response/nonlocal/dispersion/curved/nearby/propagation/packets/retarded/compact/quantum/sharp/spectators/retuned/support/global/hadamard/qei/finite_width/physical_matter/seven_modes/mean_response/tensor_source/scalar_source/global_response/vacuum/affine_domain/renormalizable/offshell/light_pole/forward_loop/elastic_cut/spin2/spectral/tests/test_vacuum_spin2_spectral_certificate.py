"""Complete ancestry and all source/report-field mutation checks."""

import copy
import json

import pytest
from p8_vacuum_spin2_spectral import verify

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
    "actual_positive_vertex_spectral_representation",
    "both_actual_unequal_mass_cut_densities",
    "independent_spin_two_phase_space_and_positive_series",
    "strict_low_window_and_full_species_moments",
    "exact_stable_spectral_calibrations",
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
        del changed["strict_low_window_and_full_species_moments"]
    else:
        changed["full_V_G_B_and_original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, d)


def test_exact_counts_and_original_scope():
    d = verify.build_report()
    assert d["claim"] == "P8-S6.116.POLYNOMIAL_VACUUM_SPIN_TWO_SPECTRAL_CUTS"
    assert d["named_exact_check_count"] == d["checked_scalar_entries"] == 52
    assert len(d["proof_checks"]) == 30
    assert d["controls"]["rejected_inputs"] == 47
    assert len(d["source_sha256"]) == 17
