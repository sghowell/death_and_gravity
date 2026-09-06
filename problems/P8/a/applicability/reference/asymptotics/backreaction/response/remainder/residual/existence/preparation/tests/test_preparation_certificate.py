import copy
import json
from fractions import Fraction

import pytest
from p8a_preparation import independent, verify


def test_checked_in_certificate_replays_exactly():
    expected = json.loads(verify.REPORT.read_text())
    verify.validate_report(expected, verify.build_report())


def test_pinned_prior_inputs_are_unchanged():
    assert verify.prior_checks() == verify.PINS


def test_report_does_not_promote_local_theorem_to_full_p8():
    data = verify.build_report()
    assert "WHOLE_TARGET_AND_P8_OPEN" in data["status"]
    assert data["theorem"]["source_globally_compact_in_time"] is False
    assert "L/2" in data["theorem"]["unforced_SEE_domain"]
    assert "ordinary radiation" in data["theorem"]["fixed_classical_source"]
    assert "same slab" in data["theorem"]["smoothness"]
    assert any("A9 QSEI" in item for item in data["not_established"])


def test_report_tampering_rejected():
    expected = verify.build_report()
    altered = copy.deepcopy(expected)
    altered["theorem"]["source_globally_compact_in_time"] = True
    with pytest.raises(ValueError):
        verify.validate_report(expected, altered)


def test_main_and_independent_majorants_agree():
    result = verify.checked_constants()
    replay = result["independent_Fraction_replay"]
    assert all(Fraction(v) > 0 for v in replay["strict_margins"].values())
    assert Fraction(replay["constants"]["highest_jet_contraction"]) < Fraction(3, 25)


def test_independent_engine_integrates_polynomials_not_samples():
    assert independent.integral({1: 3, 3: 5}, Fraction(2, 3)) == Fraction(74, 81)


def test_source_hashes_cover_independent_and_written_proofs():
    sources = verify.build_report()["source_sha256"]
    required = {"src/p8a_preparation/independent.py", "notes/regularity.md",
                "notes/contraction.md", "notes/construction.md",
                "tests/test_preparation_covariant_audit.py"}
    assert required <= sources.keys()
