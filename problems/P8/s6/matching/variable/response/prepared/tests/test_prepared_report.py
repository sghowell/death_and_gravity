import copy
import hashlib
import json

import pytest
from p8_variable_prepared import verify


@pytest.fixture(scope="module")
def candidate():
    return verify.build_report()


def test_candidate_replays_all_priors_and_exact_bridges(candidate):
    assert candidate["claim"] == "P8-S6.23.PREPARED"
    assert candidate["prior_sha256"] == {"S6_21_full_physical_response_and_its_replayed_lineage": verify.RESPONSE_SHA}
    assert candidate["independent_Fraction_replay"]["all_calibrations_equal"] is True
    assert sum(len(group) for group in candidate["exact_residuals"].values()) == 33
    assert all(value == "0" for group in candidate["exact_residuals"].values() for value in group.values())


def test_every_source_hash_is_actual_and_reserves_root_audit(candidate):
    hashes = candidate["source_sha256"]
    assert "tests/test_prepared_independent_audit.py" in hashes
    assert len(hashes) == len(verify.source_files())
    assert all(hashlib.sha256((verify.ROOT/path).read_bytes()).hexdigest() == value
               for path, value in hashes.items())
    assert all(not path.startswith("../") for path in hashes)


def test_build_report_is_read_only(candidate):
    before = {path: path.read_bytes() for path in verify.source_files()}
    candidate_again = verify.build_report()
    verify.validate_report(candidate, candidate_again)
    assert before == {path: path.read_bytes() for path in verify.source_files()}


@pytest.mark.parametrize("field", ["claim", "prior_sha256", "source_sha256", "exact_residuals", "not_established"])
def test_report_omission_or_mutation_is_rejected(candidate, field):
    wrong = copy.deepcopy(candidate)
    wrong.pop(field)
    with pytest.raises(ValueError):
        verify.validate_report(wrong, candidate)


def test_no_unscoped_health_or_completion_claim(candidate):
    assert "ORIGINAL_S6_P8_OPEN" in candidate["status"]
    assert "1<=K<=4" in candidate["unchanged_physical_parent"]["separate_transfer_domain"]
    assert any("off-shell" in item for item in candidate["not_established"])
    assert any("fixed" in item and "retuning" in item for item in candidate["not_established"])
    assert candidate["checked_controls"]["rejected_exact_domain_calls"] >= 25


def test_pinned_report_matches_read_only_replay(candidate):
    verify.validate_report(json.loads(verify.REPORT.read_text()), candidate)
