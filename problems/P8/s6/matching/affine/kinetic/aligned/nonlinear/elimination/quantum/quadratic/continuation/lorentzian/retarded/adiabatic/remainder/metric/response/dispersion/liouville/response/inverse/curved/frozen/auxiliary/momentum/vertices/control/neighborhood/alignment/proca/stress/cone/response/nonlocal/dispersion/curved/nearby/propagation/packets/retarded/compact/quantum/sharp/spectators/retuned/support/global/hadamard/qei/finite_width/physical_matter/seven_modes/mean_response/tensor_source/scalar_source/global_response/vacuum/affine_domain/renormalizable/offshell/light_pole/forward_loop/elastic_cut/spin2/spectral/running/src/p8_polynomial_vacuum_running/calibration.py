"""Rational enclosure of the actual one-loop flow on the declared reference window."""

from fractions import Fraction

import sympy as sp
from p8_polynomial_vacuum import model


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational logarithmic reference time")
    return sp.Rational(value)


def point(logarithmic_time=9):
    ell = rational(logarithmic_time)
    if not 0 <= ell <= 9:
        raise ValueError("Require logarithmic reference time in [0,9]")
    p = model.data()["actual_parameters"]
    L = p["bare_polynomial_quartic"]
    g = p["cubic_coupling_squared"]
    M = p["heavy_mass_squared"]
    d = p["positive_completed_square_quartic_margin"]
    den = 1 - 3 * L * ell
    return {
        "logarithmic_reference_time": ell,
        "flow_root_cubed": den,
        "exact_running_polynomial_quartic": L / den,
        "cubic_squared_lower": g,
        "cubic_squared_upper": g / den**2,
        "heavy_mass_squared_lower": M + g * ell,
        "heavy_mass_squared_upper": M + 2 * g * ell,
        "completed_square_margin_lower": d,
        "completed_square_margin_upper": d + 12 * L * L * ell,
        "margin_relative_increase_upper": 12 * L * L * ell / d,
        "scope": "Enclosures of the exact displayed one-loop beta ODE. They do not enclose omitted higher-loop flow or change the fixed physical light pole mass.",
    }


def data():
    p = model.data()["actual_parameters"]
    L = p["bare_polynomial_quartic"]
    g = p["cubic_coupling_squared"]
    M = p["heavy_mass_squared"]
    d = p["positive_completed_square_quartic_margin"]
    end = point()
    start = point(0)
    delta, q = sp.symbols(
        "positive_denominator_decrement positive_flow_root", positive=True
    )
    exp_three_lower = sum(sp.Rational(3) ** j / sp.factorial(j) for j in range(4))
    return {
        "actual_initial_polynomial_quartic": L,
        "actual_initial_cubic_squared": g,
        "actual_initial_heavy_mass_squared": M,
        "actual_initial_positive_margin": d,
        "selected_Planck_mass_in_light_mass_units": sp.Integer(10) ** 400,
        "whole_reference_window": "1 <= nu <= 10^400 gives ell=ln(nu)/(16pi^2)<25/3<9; every stated enclosure covers ell in [0,9]",
        "whole_window_quartic_relative_increase_upper": 54 * L,
        "whole_window_cubic_squared_relative_increase_upper": 108 * L,
        "whole_window_heavy_mass_squared_increase_upper": 18 * g,
        "whole_window_heavy_mass_squared_relative_increase_upper": 18 * g / M,
        "whole_window_margin_relative_increase_upper": 108 * L * L / d,
        "actual_end_enclosure": end,
        "checks": {
            "actual_initial_margin_same_polynomial_scheme": L - 3 * g / M - d,
            "initial_quartic_anchor": start["exact_running_polynomial_quartic"] - L,
            "initial_heavy_mass_anchor": start["heavy_mass_squared_lower"] - M,
            "whole_window_mass_increment_bound": end["heavy_mass_squared_upper"]
            - M
            - 18 * g,
            "whole_window_margin_increment_bound": end["completed_square_margin_upper"]
            - d
            - 108 * L * L,
            "root_above_its_cube_positive_gap": sp.expand(
                q - q**3 - q * (1 - q) * (1 + q)
            ),
            "quartic_relative_denominator_bound_gap": sp.factor(
                2 * delta - delta / (1 - delta) - delta * (1 - 2 * delta) / (1 - delta)
            ),
            "cubic_relative_denominator_bound_gap": sp.factor(
                4 * delta
                - (1 / (1 - delta) ** 2 - 1)
                - delta * (2 - 7 * delta + 4 * delta**2) / (1 - delta) ** 2
            ),
            "reference_logarithm_window_rational_bound": sp.Rational(1200, 144)
            - sp.Rational(25, 3),
        },
        "bounds": {
            "actual_positive_initial_completed_square_margin": d > 0,
            "denominator_decrement_below_one_quarter": 0 < 27 * L < sp.Rational(1, 4),
            "running_quartic_below_one_e_minus_204": 2 * L < sp.Rational(1, 10**204),
            "running_cubic_squared_below_one_e_minus_seven": 2 * g
            < sp.Rational(1, 10**7),
            "actual_quartic_relative_running_below_one_e_minus_202": 54 * L
            < sp.Rational(1, 10**202),
            "actual_cubic_relative_running_below_one_e_minus_202": 108 * L
            < sp.Rational(1, 10**202),
            "actual_heavy_mass_squared_increment_below_three_e_minus_seven": 18 * g
            < sp.Rational(3, 10**7),
            "actual_heavy_mass_squared_relative_increment_below_one_e_minus_203": 18
            * g
            / M
            < sp.Rational(1, 10**203),
            "actual_positive_margin_relative_increase_below_one_e_minus_five": 108
            * L
            * L
            / d
            < sp.Rational(1, 10**5),
            "selected_heavy_mass_below_selected_Planck_mass_squared": M + 18 * g
            < sp.Integer(10) ** 800,
            "exponential_positive_series_proves_log_ten_below_three": exp_three_lower
            > 10,
            "reference_time_window_strictly_within_nine": sp.Rational(25, 3) < 9,
        },
    }


def bad_cases():
    invalid = (
        True,
        False,
        1.0,
        sp.Float(1),
        "1",
        sp.I,
        sp.oo,
        -sp.oo,
        sp.zoo,
        sp.nan,
        None,
        sp.Symbol("ell"),
    )
    rows = [("inexact_time_" + str(i), point, (v,)) for i, v in enumerate(invalid)]
    rows += [
        ("outside_reference_window_" + str(i), point, (v,))
        for i, v in enumerate((-1, 10, Fraction(-1, 2), Fraction(19, 2)))
    ]
    return rows
