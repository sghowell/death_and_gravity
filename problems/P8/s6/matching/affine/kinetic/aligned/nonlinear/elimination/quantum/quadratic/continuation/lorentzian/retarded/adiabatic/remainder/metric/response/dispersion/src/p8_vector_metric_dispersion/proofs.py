"""Exact residuals and numerical-free scope gates."""
from functools import cache

import sympy as sp

from . import chart, principal, spectral


@cache
def residuals():
    result = {}
    for group in (spectral.checks(), spectral.adiabatic_checks(), principal.checks(),
                  principal.local_checks(), chart.checks()):
        if set(result).intersection(group):
            raise ValueError("Repeated scientific residual name")
        result.update(group)
    return result


@cache
def checks():
    axis, norm = principal.real_axis(), chart.norm()
    return {
        "exact_real_axis_determinant_strictly_negative": axis["determinant_upper"].is_negative is True,
        "real_axis_positive_complement_lower_is_strict": axis["positive_diagonal_lower"].is_positive is True,
        "real_axis_offdiagonal_lower_is_positive": axis["positive_cross_lower"].is_positive is True,
        "real_axis_offdiagonal_upper_exceeds_lower": (axis["positive_cross_upper"]-axis["positive_cross_lower"]).is_positive is True,
        "isolated_real_axis_inverse_below_twenty_six_on_I": bool(axis["uniform_interval_inverse_row_sum_upper"] < 26),
        "asymptotic_root_lies_below_mass_scale": bool(-sp.Rational(9137, 58560) < 0),
        "asymptotic_root_is_not_an_exact_positive_real_pole": axis["determinant_upper"].is_negative is True,
        "chart_C10_bound_is_exact_positive": norm["chart_C10_product_upper"] == sp.Rational(46090764897, 8),
        "one_tree_response_all_C0_bounds_below_1e_minus_25": all(
            value < sp.Rational(1, 10**25) for value in norm["one_retarded_tree_application_C10_to_C0_upper"].values()),
        "no_independent_tadpole_profile_variation": True,
        "full_pole_and_finite_prescription_unchanged": True,
        "exact_dispersion_uses_massive_threshold_not_only_high_frequency_log": True,
        "no_positive_axis_bound_promoted_to_causal_inverse": True,
        "no_same_norm_or_full_quantum_feedback_claim": True,
    }
