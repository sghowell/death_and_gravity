import copy
import json

import pytest
from p8_m1_loops import verify


def test_report_replays_scoped_certificate_and_prior_lineage():
    actual = verify.build_report()
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual)
    assert actual["prior_S5_8_CD_sha256"] == verify.PRIOR_SHA
    assert actual["claim"] == "P8-S5.9.CD"
    assert "FULL_QUANTUM_CONTROL_OPEN" in actual["status"]
    assert "tests/test_m1_loops_lapse_variation.py" in actual["source_sha256"]
    assert "tests/test_m1_loops_weyl_variation.py" in actual["source_sha256"]


def test_report_mutation_cannot_promote_isolated_log_to_total_loop_control():
    original = {"status": "isolated", "finite_matching_coefficients_bounded": False}
    mutated = copy.deepcopy(original)
    mutated["finite_matching_coefficients_bounded"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(original, mutated)


def test_every_algebraic_and_omission_control_is_replayed():
    assert all(value == 0 for group in verify.residuals().values() for value in group.values())
    controls = verify.control_checks()
    assert controls["heat"]["conformal_R2_coefficient"] == "0"
    assert controls["geometry"]["bulk_A2_variation_bounce_not_zero"] == "4"
