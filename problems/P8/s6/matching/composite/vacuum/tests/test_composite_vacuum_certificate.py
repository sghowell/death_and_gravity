import copy
import json

import pytest
from p8_composite_vacuum import verify


def test_frozen_report_equals_exact_read_only_replay():
    actual = verify.build_report()
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual)
    assert actual["claim"] == "P8-S6.12.COMPOSITE"
    assert len(actual["exact_residuals"]) == 5
    assert sum(map(len, actual["exact_residuals"].values())) == 33
    assert "tests/test_composite_vacuum_covariant_audit.py" in actual["source_sha256"]


def test_certificate_rejects_a_false_heavy_physical_source_residue():
    actual = verify.build_report()
    bad = copy.deepcopy(actual)
    bad["conditional_proportional_vacuum"]["physical_source_response"] = "nonzero massive physical source residue"
    with pytest.raises(ValueError):
        verify.validate_report(bad, actual)


def test_certificate_keeps_all_three_original_source_pins():
    pins = verify.prior_checks()
    assert pins == {"S6_11_composite_response": verify.PRIOR_SHA,
                    "original_32_row_classification": verify.CLASSIFICATION_SHA,
                    "original_CD_matter_witness": verify.CD_M1_SHA}


def test_certificate_controls_keep_resonance_and_the_approximate_budget_domain():
    data = verify.controls()["checked_values"]
    assert data["nonzero_solution_with_zero_Dirichlet_data"] == "1"
    assert data["oversized_chi_error_rejected"] == "True"
    assert data["exact_actual_null_residual_budget"] == "801/100"
    assert data["same_potential_gap_for_epsilon_at_most_1_over_10000"] == "3749/10000"
