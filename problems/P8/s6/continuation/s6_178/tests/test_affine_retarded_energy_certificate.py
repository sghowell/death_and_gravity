"""Immutable exact report, explicit manifest and every-field mutation rejection."""

import copy
import json

import pytest
from p8_vacuum_affine_retarded_energy import verify

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
    "physical_retarded_energy_and_full_source_bounds",
    "conditional_mean_energy_and_causal_light_force",
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
        del bad["conditional_mean_energy_and_causal_light_force"]
    else:
        bad["quantum_parent_matching_complete"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(bad, d)


def test_exact_counts_and_scope():
    d = verify.build_report()
    assert len(d["source_sha256"]) == 18
    assert d["named_exact_check_count"] == 51
    assert d["checked_scalar_entries"] == 51
    assert len(d["proof_checks"]) == 41
    assert len(d["controls"]) == 9
    assert d["controls"]["rejected_inputs"] == 100
    assert "NOT_METRIC_NOISE_FULL_PARENT_CUTOFF" in d["status"]
