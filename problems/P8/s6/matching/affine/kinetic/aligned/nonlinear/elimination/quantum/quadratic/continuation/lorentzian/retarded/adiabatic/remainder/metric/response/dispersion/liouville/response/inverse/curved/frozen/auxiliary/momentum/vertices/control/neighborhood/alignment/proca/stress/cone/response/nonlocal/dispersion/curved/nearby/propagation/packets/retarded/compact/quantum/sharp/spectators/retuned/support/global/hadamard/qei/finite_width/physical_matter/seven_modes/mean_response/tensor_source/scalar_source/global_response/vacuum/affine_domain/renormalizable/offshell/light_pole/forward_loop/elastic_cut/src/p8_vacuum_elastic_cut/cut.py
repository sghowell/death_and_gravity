"""Actual finite low-energy elastic cut and controlled one-loop subtraction."""

from functools import cache

import sympy as sp
from p8_vacuum_forward_loop import calibration as loop

from . import amplitude, angular


@cache
def data():
    d = amplitude.data()
    a = angular.data()
    s = d["s"]
    M = d["M"]
    g = d["g"]
    quartic = d["lambda4"]
    k = (s - 4) / 2
    B = M + k
    C = -quartic + g / (M - s)
    average = a["exact_even_squared_vertex_average"].subs(
        {a["B"]: B, a["k"]: k, a["g"]: g, a["C"]: C}, simultaneous=True
    )
    beta = sp.sqrt(1 - 4 / s)
    rho = beta * average / (32 * sp.pi)
    lam = d["actual_fixed_lambda"]
    lbound = sp.Rational(2, 5) * 42**2 / (32 * 4)
    ubound = sp.Rational(3, 5) * 73**2 / (32 * 3)
    wfull = sp.integrate(1 / (s - 2) ** 3, (s, 4, 6))
    wsub = sp.integrate(1 / (s - 2) ** 3, (s, 5, 6))
    v = sp.Symbol("crossing_center_variable", real=True)
    dispersive = 2 * v * v / (sp.pi * (s - 2) * ((s - 2) ** 2 - v * v))
    prev = loop.point()
    cut_upper = 3 * lam * lam
    remaining = (
        prev["actual_tree_b2"] - prev["total_one_loop_b2_error_upper"] - cut_upper
    )
    return {
        "physical_s": s,
        "exact_elastic_one_loop_cut_density_above_threshold": rho,
        "cut_density_at_threshold": sp.Integer(0),
        "selected_low_energy_cut": sp.Integral(
            2 * rho / (sp.pi * (s - 2) ** 3), (s, 4, 6)
        ),
        "forward_crossing_dispersion_kernel": dispersive,
        "actual_strict_cut_lower": lam * lam / 20,
        "actual_cut_upper": cut_upper,
        "actual_tree_plus_one_loop_minus_cut_lower": remaining,
        "low_energy_cut_window": [sp.Integer(4), sp.Integer(6)],
        "strict_density_subwindow": [sp.Integer(5), sp.Integer(6)],
        "cut_density_enclosures": "0<=rho<34 lambda^2 on 4<=s<=6; rho>5 lambda^2 on 5<=s<=6. At s=4 it vanishes by phase space. The full exact nonlocal tree amplitude, not its contact truncation, is squared.",
        "scope": "Actual complete one-loop elastic discontinuity on the finite window, and a strictly positive one-loop cut-subtracted coefficient. No all-orders error, global dispersive equality, high-energy contour, finite gravity or common bounce matching is inferred.",
        "checks": {
            "strict_density_lower_rational_constant": lbound - sp.Rational(441, 80),
            "uniform_density_upper_rational_constant": ubound - sp.Rational(5329, 160),
            "full_cut_window_kernel_integral": wfull - sp.Rational(3, 32),
            "positive_subwindow_kernel_integral": wsub - sp.Rational(7, 288),
            "coefficient_not_second_derivative_dispersion_normalization": sp.factor(
                sp.diff(dispersive, v, 2).subs(v, 0) / 2 - 2 / (sp.pi * (s - 2) ** 3)
            ),
            "coarse_positive_low_cut_lower_constant": sp.Rational(1, 2) * 5 * wsub
            - sp.Rational(35, 576),
            "coarse_positive_low_cut_upper_constant": sp.Rational(2, 3) * 34 * wfull
            - sp.Rational(17, 8),
        },
        "bounds": {
            "phase_space_lower_on_s_five_to_six": sp.Rational(1, 5)
            > sp.Rational(2, 5) ** 2,
            "phase_space_upper_on_s_four_to_six": sp.Rational(1, 3)
            < sp.Rational(3, 5) ** 2,
            "strict_density_lower_above_five_lambda_squared": lbound > 5,
            "uniform_density_upper_below_thirty_four_lambda_squared": ubound < 34,
            "strict_partial_cut_lower_above_one_twentieth_lambda_squared": sp.Rational(
                35, 576
            )
            > sp.Rational(1, 20),
            "partial_cut_upper_below_three_lambda_squared": sp.Rational(17, 8) < 3,
            "actual_tree_plus_complete_one_loop_minus_computed_cut_positive": remaining
            > 0,
            "same_one_millionth_relative_margin_after_cut_subtraction": prev[
                "total_one_loop_b2_error_upper"
            ]
            + cut_upper
            < prev["actual_tree_b2"] / 10**6,
        },
    }
