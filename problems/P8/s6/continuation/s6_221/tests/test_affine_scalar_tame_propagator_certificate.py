"""Immutable native report, exact manifests and every-field rejection."""

import copy
import json

import pytest
from p8_vacuum_affine_scalar_tame_propagator import verify

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
    "actual_coefficients_complete_charts_and_uniform_majorants",
    "energy_and_all_transfer_tame_phase_comparison",
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
        str(p.relative_to(verify.ROOT)): verify.sha(p) for p in verify.source_files()
    }


@pytest.mark.parametrize("key", KEYS)
def test_every_report_field_mutation_rejected(key):
    data = verify.build_report()
    bad = copy.deepcopy(data)
    bad[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_and_extra_fields_rejected(mutation):
    data = verify.build_report()
    bad = copy.deepcopy(data)
    if mutation == "missing":
        del bad["energy_and_all_transfer_tame_phase_comparison"]
    else:
        bad["same_space_quantum_inverse"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_exact_counts_and_uncompleted_quantum_boundary():
    data = verify.build_report()
    assert len(data["source_sha256"]) == 18
    assert (
        data["named_exact_check_count"] == 45 and data["checked_scalar_entries"] == 126
    )
    assert len(data["proof_checks"]) == 51
    assert len(data["controls"]) == 9 and data["controls"]["rejected_inputs"] == 176
    assert "NOT_SAME_SPACE_QUANTUM_INVERSE" in data["status"]
    assert (
        data["prior_sha256"]["S6_182_fixed_QG1_and_original_five_stress_jets_unchanged"]
        is True
    )
