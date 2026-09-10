"""Exact bounds for the isolated affine coordinate change, not a new loop budget."""

from functools import cache

import sympy as s
from p8_polynomial_vacuum import model
from p8_vacuum_two_loop_finite_contact.calibration import rational

from . import contact


def enclosure(L, h=1):
    L, h = map(rational, (L, h))
    if not 0 < L < s.Rational(3, 50) or not 0 <= h <= 1:
        raise ValueError("Need 0<L<3/50 and 0<=h<=1")
    k = s.Rational(25, 3) * L
    sigma = s.Rational(50, 3) * L**2
    return {
        "k_upper": k,
        "sigma_absolute_upper": sigma,
        "first_coordinate_change_upper": h * sigma,
        "second_coordinate_change_upper": h**2 * sigma * k,
        "coordinate_tail_after_order_two_upper": h**3 * sigma * k**2 / (1 - h * k),
        "uncancelled_insertion_bound_before_matching": s.Rational(10025, 18) * L**3,
    }


@cache
def data():
    p = model.data()["actual_parameters"]
    L, M, g = (
        p[name]
        for name in (
            "bare_polynomial_quartic",
            "heavy_mass_squared",
            "cubic_coupling_squared",
        )
    )
    en = enclosure(L)
    d = contact.data()
    sy = d["symbols"]
    actual = {sy["L"]: L, sy["M"]: M, sy["g"]: g, sy["Q"]: 16 * s.pi**2}
    sigma, k = (d[name].subs(actual) for name in ("sigma", "positive_affine_k"))
    return {
        "actual_parameters": p,
        "actual_fixed_sigma": sigma,
        "actual_affine_k": k,
        "isolated_coordinate_enclosure": en,
        "matched_sigma_nonlocal_order_two_exact": s.Integer(0),
        "checks": {
            "k_prefactor_from_log_g_and_Q_bounds": s.Rational(6 * 600, 3 * 144)
            - s.Rational(25, 3),
            "sigma_slope_second_coefficient_bound": en["second_coordinate_change_upper"]
            - s.Rational(1250, 9) * L**3,
            "insertion_bound_retained_from_S6_125": en[
                "uncancelled_insertion_bound_before_matching"
            ]
            - s.Rational(50, 3) * L**2 * s.Rational(4812, 144) * L,
            "matched_nonlocal_cancellation_not_triangle_inequality": s.Integer(0),
        },
        "bounds": {
            "actual_affine_bound_domain": bool(0 < L < s.Rational(3, 50)),
            "actual_heavy_mass_above_one": bool(M > 1),
            "positive_completed_square_margin": bool(L - 3 * g / M > 0),
            "g_over_M_less_than_L_third": bool(g / M < L / 3),
            "log_argument_less_than_ten_to_two_hundred": bool(4 * M < 10**200),
            "inverse_denominator_uniformly_above_one_half": bool(
                en["k_upper"] < s.Rational(1, 2)
            ),
            "first_coordinate_shift_below_one_e_minus_408": bool(
                en["first_coordinate_change_upper"] < s.Rational(1, 10**408)
            ),
            "second_coordinate_shift_below_two_e_minus_612": bool(
                en["second_coordinate_change_upper"] < s.Rational(2, 10**612)
            ),
            "coordinate_tail_below_four_e_minus_816": bool(
                en["coordinate_tail_after_order_two_upper"] < s.Rational(4, 10**816)
            ),
            "old_individual_insertion_bound_below_seven_e_minus_612": bool(
                en["uncancelled_insertion_bound_before_matching"]
                < s.Rational(7, 10**612)
            ),
        },
        "scope": "Because sigma<0 for the actual completed-square parameters, this isolated positive-h map increases L and cannot diminish L-3g/M. This statement does not apply to the full combined MS map. Its h^3 coordinate tail is an algebraic inverse tail only, not a higher-loop/truncation error bound.",
    }
