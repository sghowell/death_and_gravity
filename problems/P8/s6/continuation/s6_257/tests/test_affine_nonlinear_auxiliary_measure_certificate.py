"""Frozen complete nonlinear auxiliary and measure certificate contract."""

import copy
import json

import pytest
from p8_vacuum_affine_nonlinear_auxiliary_measure import verify

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
    "whole_current_nonlinear_canonical_Hamiltonian_and_local_auxiliary_branch",
    "whole_cotangent_boundary_second_class_measure_and_source_centering",
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
        del bad[
            "whole_current_nonlinear_canonical_Hamiltonian_and_local_auxiliary_branch"
        ]
    else:
        bad["original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_exact_counts_and_unchanged_quantum_and_cutoff_boundary():
    data = verify.build_report()
    assert len(data["source_sha256"]) == 17
    assert (data["named_exact_check_count"], data["checked_scalar_entries"]) == (
        53,
        451,
    )
    assert len(data["proof_checks"]) == 37
    assert len(data["controls"]) == 9 and data["controls"]["rejected_inputs"] == 294
    assert (
        "NOT_FULL_GAUGE_AFFINE_COVARIANT_QUANTUM_MEASURE_FIXED_MEAN_ORIGINAL_V_G_B_OR_P8"
        in data["status"]
    )
    assert (
        data["prior_sha256"][
            "classical_comparison_not_replacement_fixed_quantum_preparation"
        ]
        is True
    )
