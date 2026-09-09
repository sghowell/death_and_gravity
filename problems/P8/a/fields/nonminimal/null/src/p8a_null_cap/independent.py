"""Independent Fraction proof constants, with no production or SymPy imports."""

from fractions import Fraction as F


def replay():
    e, w, err = F(1, 1000), F(1, 10000), F(1, 100)
    second = 12 * e * e + 5 * e
    value, first = 3 + second / 2, 1 + second
    errors = (
        second / 2 + w * w * value / 2,
        second + w * w * first / 2 + w * w * value,
        second + 2 * w * w * first + w * w * value,
    )
    negative = (1 - err) ** 2 - 3 * err * err
    correction = (3 + err) * err + 3 * err * err
    wick = -2 * negative / 3
    null = (-negative + correction) / 3
    conformal = (-4 * negative + 2 * correction) / 9
    bounds = {
        "epsilon": e,
        "bump_width": w,
        "discrete_second_derivative_bound": second,
        "mollified_value_error": errors[0],
        "mollified_first_error": errors[1],
        "mollified_second_upper": errors[2],
        "common_jet_error_bound": err,
        "strict_error_margins": [err - v for v in errors],
        "Wick_square_upper_per_mode_amplitude_squared": wick,
        "null_stress_upper_all_xi_zero_to_quarter": null,
        "null_stress_upper_at_conformal_xi": conformal,
        "negative_Wick_margin_below_minus_three_fifths": -F(3, 5) - wick,
        "negative_uniform_null_margin_below_minus_three_tenths": -F(3, 10) - null,
        "negative_conformal_null_margin_below_minus_two_fifths": -F(2, 5) - conformal,
        "positive_low_frequency_support_margin": e - w,
        "disjoint_frequency_support_margin": e - 2 * w,
    }
    Iq, I0, I2max = (
        24 * e + 3 / e,
        20 + 2 / e**2,
        32 * e**2 + 5 + (20 + 2 / e**2) * w * w,
    )
    resources = {
        "energy_lower_coefficient_of_hbar_R_squared": I0 / (6 * Iq),
        "energy_upper_coefficient_of_hbar_R_squared": 2 * I0 / (3 * Iq),
        "energy_upper_constant_coefficient_of_hbar": I2max / (6 * Iq),
        "positive_ANEC_support_sum_gap": 2 * (e - w),
        "each_fixed_R_has_finite_energy_and_smooth_Hadamard_difference": True,
        "uniform_total_energy_or_transverse_momentum_cutoff_for_the_family": False,
        "global_two_sided_Wick_square_cap_for_the_family": False,
    }
    assert (
        max(errors) < err
        and wick < -F(3, 5)
        and null < -F(3, 10)
        and conformal < -F(2, 5)
    )
    return {"profile_bounds": bounds, "resource_bounds": resources}
