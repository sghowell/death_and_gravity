import copy
import json
from fractions import Fraction as Q

import pytest
from p8a_existence import independent, log_inverse, mode_lipschitz, verify


def test_pinned_response_inverse_certificate_replay():
    report = verify.build_report()
    verify.validate_report(json.loads(verify.REPORT.read_text()), report)
    assert report["claim"] == "P8-A.10"
    assert report["prior_sha256"]["A9"] == verify.PRIOR_SHA


@pytest.mark.parametrize("key", ["mode_theorem", "inverse_theorem", "conditional_fixed_point_gate", "not_established"])
def test_report_cannot_promote_unproved_SEE_or_change_the_operator(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "P8 solved without compatible data or interval enclosure"
    with pytest.raises(ValueError, match="differs from exact replay"):
        verify.validate_report(changed, actual)


def test_independent_exact_division_and_calibration_margins():
    assert independent.exact_polynomial_checks()["pi_integral_polynomial_part"] == "22/7"
    assert independent.inverse_calibration()["rounded_gap"] > 0
    assert independent.inverse_calibration()["abstract_distance"] == Q(3, 20)
    assert independent.inverse_calibration()["norm_upper"] == log_inverse.dyadic_norm_bound(1, 512)["norm_upper"]


def test_independent_two_mode_allocations_use_both_frequencies():
    for n in range(3, 13):
        separate = independent.dyson(n, Q(2), Q(3, 5))
        primary = mode_lipschitz.dyson_coefficients(n, 2, "3/5")
        assert {key: str(value) for key, value in separate.items()} == {
            key: str(primary[key]) for key in separate}


def test_literal_first_prepared_derivative_jet_is_checked_against_pin():
    original = json.loads(verify.remainder.REPORT.read_text())
    changed = copy.deepcopy(original)
    changed["derived_finite_geometry"]["U_derivative_bounds_through_5_per_delta"][1] = "1/2"
    assert independent.mode_calibration(changed) != independent.mode_calibration(original)
    constants = verify.checked_constants()
    assert constants["independent_Fraction_replay"]["mode_calibration"]["prepared_U1"] == "61013499/8192"
