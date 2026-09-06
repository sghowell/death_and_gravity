import copy
import json

import pytest
from p8a_extended_qsei import verify
from p8a_extension import verify as actual
from p8a_see_qsei import verify as qsei_prior


@pytest.fixture(scope="module")
def replayed():
    return verify.build_report()


def test_stored_certificate_replays(replayed):
    verify.validate_report(json.loads(verify.REPORT.read_text()), replayed)


def test_both_distinct_analytic_ancestors_are_pinned(replayed):
    assert replayed["prior_sha256"] == {"A14": verify.sha(actual.REPORT), "A12": verify.sha(qsei_prior.REPORT)}
    assert len(replayed["exact_residuals"]) == 6
    assert all(value == "0" for group in replayed["exact_residuals"].values() for value in group.values())


def test_complete_child_source_hashes_are_current(replayed):
    paths = replayed["source_sha256"]
    assert "tests/test_extended_qsei_covariant_audit.py" in paths
    assert "notes/proof.md" in paths and "FORMULATION.md" in paths
    for relative, expected in paths.items():
        assert verify.sha(verify.ROOT/relative) == expected


@pytest.mark.parametrize("key", ["status", "theorem", "derived_constants"])
def test_report_mutation_is_rejected(replayed, key):
    changed = copy.deepcopy(replayed)
    changed[key] = "unproved stronger statement"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, replayed)


def test_target_and_completion_boundary(replayed):
    theorem = replayed["theorem"]
    assert "every Hadamard" in theorem["target_states"]
    assert theorem["old_short_interval_coefficient_transferred_without_recalibration"] is False
    assert "no positive" in theorem["reference_EED_bound"]
    assert replayed["focusing_boundary"]["sufficient_index_trigger_available"] is False
    assert "P8_OPEN" in replayed["status"]
