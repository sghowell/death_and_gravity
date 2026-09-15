"""Immutable complete-source and twenty-field corrected hybrid report contract."""

import copy
import json

import pytest
from p8_vacuum_affine_canonical_boundary_corrected_hybrid import verify

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
    "whole_canonical_boundary_erratum_source_and_corrected_geometry",
    "whole_corrected_quantum_hybrid_and_regulator_comparisons",
    "observable_and_scope",
    "primitive_and_matching_frontier",
    "controls",
    "verdict",
    "not_established",
    "verification_boundary",
]


def test_complete_read_only_native_report_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_exact_complete_source_and_twenty_field_manifest():
    data = verify.build_report()
    assert list(data) == KEYS
    assert data["source_sha256"] == {
        str(path.relative_to(verify.ROOT)): verify.sha(path)
        for path in verify.source_files()
    }


@pytest.mark.parametrize("key", KEYS)
def test_every_whole_corrected_report_field_mutation_is_rejected(key):
    data = verify.build_report()
    bad = copy.deepcopy(data)
    bad[key] = "not the certified corrected value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_and_extra_corrected_report_fields_are_rejected(mutation):
    data = verify.build_report()
    bad = copy.deepcopy(data)
    if mutation == "missing":
        del bad["whole_canonical_boundary_erratum_source_and_corrected_geometry"]
    else:
        bad["original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_exact_counts_and_explicit_archived_physical_qualification():
    data = verify.build_report()
    assert len(data["source_sha256"]) == 19
    assert (data["named_exact_check_count"], data["checked_scalar_entries"]) == (
        111,
        1493,
    )
    assert len(data["proof_checks"]) == 94
    assert len(data["controls"]) == 13 and data["controls"]["rejected_inputs"] == 623
    assert "NOT_ARCHIVED_SAME_PHYSICAL_STATE_OR_ORIGINAL_V_G_B_P8" in data["status"]
    assert (
        data["prior_sha256"][
            "archived_numeric_replay_does_not_reendorse_refuted_physical_dictionary"
        ]
        is True
    )
    assert (
        data["prior_sha256"][
            "new_corrected_solutions_not_identical_to_archived_solutions"
        ]
        is True
    )
    assert (
        data["prior_sha256"][
            "S261_refutation_and_S265_integrability_boundary_unchanged"
        ]
        is True
    )
    assert len(data["observable_and_scope"]["historical_qualification"]) == 6
