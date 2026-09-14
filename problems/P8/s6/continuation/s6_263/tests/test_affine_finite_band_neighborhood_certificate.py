"""Frozen complete classical neighborhood and finite-band linear energy certificate contract."""

import copy
import json

import pytest
from p8_vacuum_affine_finite_band_neighborhood import verify

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
    "whole_complete_sources_classical_neighborhood_and_physical_turn",
    "whole_full_reference_energy_and_finite_band_variational_bounds",
    "observable_and_scope",
    "primitive_and_matching_frontier",
    "controls",
    "verdict",
    "not_established",
    "verification_boundary",
]


def test_complete_read_only_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_exact_source_and_field_manifest():
    data = verify.build_report()
    assert set(data) == set(KEYS)
    assert data["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path)
        for path in verify.source_files()
    }


@pytest.mark.parametrize("key", KEYS)
def test_every_report_field_mutation_rejected(key):
    data = verify.build_report()
    bad = copy.deepcopy(data)
    bad[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_and_extra_report_fields_rejected(mutation):
    data = verify.build_report()
    bad = copy.deepcopy(data)
    if mutation == "missing":
        del bad["whole_complete_sources_classical_neighborhood_and_physical_turn"]
    else:
        bad["original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_exact_counts_and_unchanged_quantum_and_cutoff_boundary():
    data = verify.build_report()
    assert len(data["source_sha256"]) == 18
    assert (data["named_exact_check_count"], data["checked_scalar_entries"]) == (
        44,
        6613,
    )
    assert len(data["proof_checks"]) == 98
    assert len(data["controls"]) == 9 and data["controls"]["rejected_inputs"] == 288
    assert (
        "NOT_INTERACTING_QUANTUM_MEAN_REGULATOR_CUTOFF_NONLINEAR_INHOMOGENEOUS_B_ORIGINAL_V_G_B_OR_P8"
        in data["status"]
    )
    assert (
        data["prior_sha256"][
            "distinct_classical_family_not_fixed_interacting_quantum_mean"
        ]
        is True
    )
