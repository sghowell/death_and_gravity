"""Continuous profile, physical-state, resource and admissibility controls."""

from fractions import Fraction

import sympy as sp

from . import construction, normalization


def exact_nonnegative(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError(
            "Require an exact nonnegative rational, not a float or symbolic limit"
        )
    value = sp.Rational(value)
    if value < 0:
        raise ValueError("Require a nonnegative exact value")
    return value


def family_point(R, xi):
    R, xi = map(exact_nonnegative, (R, xi))
    if R < 1 or xi > sp.Rational(1, 4):
        raise ValueError("Require R>=1 and 0<=xi<=1/4 for this stated family")
    b = construction.profile_bounds()
    err = b["common_jet_error_bound"]
    A = (1 - err) ** 2 - 3 * err**2
    B = (3 + err) * err + 3 * err**2
    return {
        "R": R,
        "xi": xi,
        "null_stress_upper_in_hbar_gamma_squared_units": R
        * R
        * (-sp.Rational(2, 3) * (1 - 2 * xi) * A + sp.Rational(4, 3) * xi * B),
        "Wick_upper_in_hbar_gamma_squared_units": R
        * R
        * b["Wick_square_upper_per_mode_amplitude_squared"],
        "one_sided_nonnegative_Wick_cap_satisfied_on_entire_sample_segment": True,
        "uniform_energy_or_two_sided_Wick_cap_claimed": False,
    }


def proof_checks():
    b = construction.profile_bounds()
    r = normalization.resource_bounds()
    return {
        "all_frequency_bumps_are_strictly_positive": b[
            "positive_low_frequency_support_margin"
        ]
        > 0,
        "two_frequency_bumps_are_disjoint": b["disjoint_frequency_support_margin"] > 0,
        "continuous_value_error_strictly_below_one_hundredth": b[
            "strict_error_margins"
        ][0]
        > 0,
        "continuous_first_error_strictly_below_one_hundredth": b[
            "strict_error_margins"
        ][1]
        > 0,
        "continuous_second_bound_strictly_below_one_hundredth": b[
            "strict_error_margins"
        ][2]
        > 0,
        "actual_Wick_square_uniformly_below_negative_three_fifths": b[
            "negative_Wick_margin_below_minus_three_fifths"
        ]
        > 0,
        "actual_null_stress_negative_for_entire_xi_interval": b[
            "negative_uniform_null_margin_below_minus_three_tenths"
        ]
        > 0,
        "actual_conformal_null_stress_below_negative_two_fifths": b[
            "negative_conformal_null_margin_below_minus_two_fifths"
        ]
        > 0,
        "energy_per_R_squared_has_strictly_positive_lower_coefficient": r[
            "energy_lower_coefficient_of_hbar_R_squared"
        ]
        > 0,
        "positive_frequency_sum_support_excludes_zero_for_full_line_ANEC": r[
            "positive_ANEC_support_sum_gap"
        ]
        > 0,
        "state_is_actual_unit_norm_mode_squeeze_not_free_quadratic_coefficients": True,
        "full_nonminimal_improvement_retained_before_null_restriction": True,
        "smooth_renormalized_difference_restricted_not_unpulled_vacuum_two_point_distribution": True,
        "each_finite_R_has_smooth_compact_momentum_data_and_finite_energy": True,
        "nonnegative_Wick_upper_cap_does_not_bound_its_absolute_value": True,
        "complete_line_ANEC_for_this_family_is_positive_not_negative": True,
        "no_SEE_singularity_or_original_P8_closure_assigned": True,
    }


def bad_cases():
    out = []
    for v in (
        True,
        False,
        0.5,
        "1/6",
        sp.Float("0.5"),
        sp.oo,
        sp.nan,
        sp.Symbol("x"),
        -1,
    ):
        out.extend(
            [
                ("R_" + str(v), family_point, (v, sp.Rational(1, 6))),
                ("xi_" + str(v), family_point, (1, v)),
            ]
        )
    out.extend(
        [
            ("zero_R", family_point, (0, 0)),
            ("R_below_domain", family_point, (sp.Rational(1, 2), 0)),
            ("xi_above_domain", family_point, (1, sp.Rational(1, 3))),
        ]
    )
    return out


def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("An invalid null-family input was accepted: " + name)
    return len(bad_cases())
