"""Read-only report replay, explicit manifests and every-field rejection."""

import copy
import json

import pytest
from p8_vacuum_affine_homogeneous_trace_anchor import verify

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
    "full_trace_Hamiltonian_domination_and_original_comparison",
    "general_dimension_volume_current_and_actual_homogeneous_anchor",
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
    d = verify.build_report()
    assert set(d) == set(KEYS)
    assert d["source_sha256"] == {
        str(p.relative_to(verify.ROOT)): verify.sha(p) for p in verify.source_files()
    }


@pytest.mark.parametrize("key", KEYS)
def test_every_report_field_mutation_rejected(key):
    d = verify.build_report()
    bad = copy.deepcopy(d)
    bad[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, d)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_or_extra_field_rejected(mutation):
    d = verify.build_report()
    bad = copy.deepcopy(d)
    if mutation == "missing":
        del bad["general_dimension_volume_current_and_actual_homogeneous_anchor"]
    else:
        bad["full_interacting_UV_completion"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, d)


def test_exact_counts_and_full_parent_scope():
    d = verify.build_report()
    assert len(d["source_sha256"]) == 19
    assert d["named_exact_check_count"] == 56 and d["checked_scalar_entries"] == 110
    assert len(d["proof_checks"]) == 46
    assert len(d["controls"]) == 9 and d["controls"]["rejected_inputs"] == 167
    assert "NOT_NONZERO_TRANSFER_SCALAR_REMAINDERS_REDUCED_INVERSE" in d["status"]
    assert (
        d["prior_sha256"][
            "same_actual_reference_Hessian_not_full_nonlinear_source_free_parent"
        ]
        is True
    )
