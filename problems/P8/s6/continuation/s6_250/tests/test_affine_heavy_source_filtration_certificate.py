"""Immutable full conditional source report and every-field mutation controls."""

import copy
import json

import pytest
from p8_vacuum_affine_heavy_source_filtration import verify

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
    "whole_conditional_source_and_influence",
    "formal_filtration_and_counterterm_boundary",
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
        del bad["formal_filtration_and_counterterm_boundary"]
    else:
        bad["full_interacting_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, data)


def test_exact_counts_and_conditional_formal_boundary():
    data = verify.build_report()
    assert len(data["source_sha256"]) == 17
    assert (data["named_exact_check_count"], data["checked_scalar_entries"]) == (
        152,
        201,
    )
    assert len(data["proof_checks"]) == 26
    assert len(data["controls"]) == 9 and data["controls"]["rejected_inputs"] == 254
    assert (
        "EXACT_CONDITIONAL_GAUSSIAN_SOURCE_INFLUENCE_AND_FORMAL_ANCESTRY_PRESERVING_LOOP_FILTRATION"
        in data["status"]
    )
    assert (
        "NOT_COMPLETE_CURVED_COUNTERFUNCTIONAL_INTERACTING_STATE_NONZERO_RENORMALIZED_COEFFICIENT_OMITTED_LOOP_BOUND_NONLINEAR_UV_REGGE_ORIGINAL_V_G_B_OR_P8"
        in data["status"]
    )
    assert (
        data["prior_sha256"][
            "missing_curved_counterfunctional_not_inferred_from_flat_jet_prescription"
        ]
        is True
    )
