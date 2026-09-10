"""Strict low-cut slope moment and non-negligible heavy-threshold contribution."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model
from p8_vacuum_spin2 import triangles


@cache
def data():
    x, y, u, K = sp.symbols(
        "parameter positive_y positive_u positive_cut_gap", positive=True
    )
    M, g = sp.symbols(
        "positive_heavy_mass_squared positive_cubic_squared", positive=True
    )
    d1 = x * M + (1 - x) ** 2
    b = x * (1 - x) / d1
    T = sp.Symbol("positive_transfer", positive=True)
    map_u = M * y * (1 + y)
    jacobian_density = y * y / ((1 + y) ** 3 * M * (1 + 2 * y))
    low_kernel_upper = u * u / M**3
    ratio = sp.factor(jacobian_density / (low_kernel_upper.subs(u, map_u)))
    beta_integral = sp.integrate(u * u / sp.sqrt(K - u), (u, 0, K))
    low_lower = g / (11059200 * M**3)
    low_upper = g / (17280 * M**3)
    slope_lower = g / (46080 * M * M)
    p = model.data()["actual_parameters"]
    Mp = p["heavy_mass_squared"]
    gp = p["cubic_coupling_squared"]
    actual = {M: Mp, g: gp}
    S = triangles.data()["actual_spin_two_slope_rational_upper"]
    splitL = sp.Integral(g * (1 - x) * b * b / (96 * sp.pi**2), (x, 0, 1))
    splitH = sp.Integral(g * x * b * b / (96 * sp.pi**2), (x, 0, 1))
    return {
        "M": M,
        "g": g,
        "T": T,
        "light_cut_density_parameter_reduction": jacobian_density,
        "positive_cut_beta_integral": beta_integral,
        "low_window_slope_moment_lower": low_lower,
        "low_window_slope_moment_upper": low_upper,
        "total_slope_lower": slope_lower,
        "low_window_fraction_upper": sp.Rational(8, 3) / M,
        "light_full_slope_moment": splitL,
        "heavy_full_slope_moment": splitH,
        "both_species_fraction_strict_lower": sp.Rational(1, 20),
        "actual_low_window_slope_moment_lower": low_lower.subs(actual),
        "actual_low_window_slope_moment_upper": low_upper.subs(actual),
        "actual_total_slope_lower": slope_lower.subs(actual),
        "actual_total_slope_upper": S,
        "actual_low_window_fraction_upper": (sp.Rational(8, 3) / M).subs(actual),
        "actual_heavy_mass_squared": Mp,
        "scope": "The low window is transfer invariant [4,6], not the full cut. The heavy cut starts at 4M. Positive form-factor spectral weight is not a full-amplitude positivity bound.",
        "checks": {
            "low_cut_parameter_map_jacobian": sp.expand(
                sp.diff(map_u, y) - M * (1 + 2 * y)
            ),
            "low_cut_density_ratio_exact": ratio - 1 / ((1 + y) ** 5 * (1 + 2 * y)),
            "strict_lower_ratio_denominator_endpoint": (
                (1 + y) ** 5 * (1 + 2 * y)
            ).subs(y, sp.Rational(1, 2))
            - sp.Rational(243, 16),
            "strict_lower_ratio_comparison_gap": sp.Rational(16, 243)
            - sp.Rational(1, 16)
            - sp.Rational(13, 3888),
            "finite_beta_parameter_integral": sp.simplify(
                beta_integral - 16 * K ** sp.Rational(5, 2) / 15
            ),
            "full_low_window_inverse_square_weight": sp.integrate(1 / T**2, (T, 4, 6))
            - sp.Rational(1, 12),
            "positive_subwindow_inverse_square_weight": sp.integrate(
                1 / T**2, (T, 5, 6)
            )
            - sp.Rational(1, 30),
            "low_moment_lower_phase_pi_and_window_factors": low_lower
            - g / (92160 * M**3) / 4 / 30,
            "low_moment_upper_phase_pi_and_window_factors": low_upper
            - g / (480 * M**3) / 3 / 12,
            "parameter_denominator_upper_mass_gap": sp.expand(
                M - d1 - (1 - x) * (M - 1 + x)
            ),
            "total_slope_lower_parameter_integral": sp.integrate(
                x * x * (1 - x) ** 2, (x, 0, 1)
            )
            - sp.Rational(1, 30),
            "low_window_fraction_bound_ratio": sp.factor(
                low_upper / slope_lower - sp.Rational(8, 3) / M
            ),
            "light_slope_lower_parameter_integral": sp.integrate(
                x * x * (1 - x) ** 3, (x, 0, 1)
            )
            - sp.Rational(1, 60),
            "heavy_slope_lower_parameter_integral": sp.integrate(
                x**3 * (1 - x) ** 2, (x, 0, 1)
            )
            - sp.Rational(1, 60),
            "species_lower_to_total_upper_fraction": sp.Rational(288, 5760)
            - sp.Rational(1, 20),
            "two_species_parameter_slope_sum": sp.factor(
                g * (1 - x) * b * b / (96 * sp.pi**2)
                + g * x * b * b / (96 * sp.pi**2)
                - g * b * b / (96 * sp.pi**2)
            ),
        },
        "bounds": {
            "actual_heavy_threshold_above_low_window": 4 * Mp > 6,
            "actual_mass_domain_for_all_positive_gaps": Mp > 1,
            "strict_positive_low_window_lower": low_lower.subs(actual) > 0,
            "low_window_lower_below_upper": low_lower.subs(actual)
            < low_upper.subs(actual),
            "actual_total_slope_lower_below_parent_upper": 0
            < slope_lower.subs(actual)
            < S,
            "actual_low_window_fraction_below_one_e_minus_196": (
                sp.Rational(8, 3) / M
            ).subs(actual)
            < sp.Rational(1, 10**196),
            "actual_low_window_slope_upper_below_one_e_minus_603": low_upper.subs(
                actual
            )
            < sp.Rational(1, 10**603),
        },
    }
