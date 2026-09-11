"""Read-only report, source manifest and all-field mutation checks."""

import copy
import json

import pytest
from p8_vacuum_full_heavy_source import verify

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
    "complete_stationary_source_ownership_and_matching",
    "mass_differentiated_sunset_and_source_bounds",
    "actual_complete_second_H_source_enclosure",
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
    changed = copy.deepcopy(d)
    changed[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, d)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_or_extra_field_rejected(mutation):
    d = verify.build_report()
    changed = copy.deepcopy(d)
    if mutation == "missing":
        del changed["complete_stationary_source_ownership_and_matching"]
    else:
        changed["original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, d)


def test_exact_counts_and_scope():
    d = verify.build_report()
    assert len(d["source_sha256"]) == 18
    assert d["named_exact_check_count"] == 44
    assert d["checked_scalar_entries"] == 44
    assert len(d["proof_checks"]) == 23
    assert d["controls"]["rejected_inputs"] == 138
    assert "NOT_OTHER_COORDINATE_DICTIONARIES_GLOBAL_POTENTIAL" in d["status"]
