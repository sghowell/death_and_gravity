"""Exact actual two-loop pole bounds on every rational subdisc of the unit disc."""

from fractions import Fraction
from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model
from p8_vacuum_light_pole import calibration as one_loop


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational disc radius")
    return sp.Rational(value)


def point(radius=1):
    r = rational(radius)
    if not 0 < r <= 1:
        raise ValueError("Require a positive radius at most one")
    d = data()
    error = r * d["actual_one_plus_two_loop_quadratic_coefficient_upper"]
    return {
        "radius_about_mass_one": r,
        "inverse_propagator_factored_error_upper": error,
        "positive_inverse_factor_modulus_lower": 1 - error,
        "unique_mass_one_pole_through_two_loops": bool(error < 1),
        "light_pole_residue_through_two_loops": sp.S.One,
        "scope": "First-sheet subdisc of the unit disc, for the fixed-order inverse propagator through two loops only. This bound gives no all-orders pole count.",
    }


@cache
def data():
    p = model.data()["actual_parameters"]
    M, g, L = [
        p[k]
        for k in (
            "heavy_mass_squared",
            "cubic_coupling_squared",
            "bare_polynomial_quartic",
        )
    ]
    old = one_loop.point(1)
    B1 = old["quadratic_self_energy_upper_coefficient"]
    alpha = old["finite_kinetic_counterterm_upper"]
    groups = {
        "local_sunset": L**2 / (3 * 20736),
        "mixed_sunset": sp.Rational(33, 2) * L * g / 20736,
        "same_heavy_sunset": sp.Rational(5, 4) * g**2 / 20736,
        "UV_finite_sunset": 28 * g**2 / 20736,
        "nested_decaying_remainder": 2 * g**2 / 20736,
        "nested_constant_alpha": alpha * B1,
    }
    B2 = sum(groups.values())
    combined = B1 + B2
    return {
        "actual_heavy_mass_squared": M,
        "actual_cubic_squared": g,
        "actual_quartic": L,
        "inherited_one_loop_quadratic_coefficient_upper": B1,
        "inherited_asymptotic_multiplier_upper": alpha,
        "actual_two_loop_group_quadratic_coefficient_uppers": groups,
        "actual_two_loop_quadratic_coefficient_upper": B2,
        "actual_one_plus_two_loop_quadratic_coefficient_upper": combined,
        "scope": "Conservative absolute uniform bounds after every required local subtraction. Heavy-mass suppression is intentionally relaxed where convergence is retained. The larger two-loop upper bound does not establish a growing perturbative series; these are unequal conservative estimates, not signed or exact coefficients. Original V/G/B and P8 remain open.",
        "checks": {
            "same_actual_one_loop_bound": B1 - g / (864 * M**2 * (1 - 1 / M)),
            "same_actual_alpha_bound": alpha - g / (288 * M),
            "independent_total_g_squared_prefactor": 28
            + sp.Rational(5, 4)
            + 2
            - sp.Rational(125, 4),
            "actual_six_group_sum": B2
            - alpha * B1
            - (L**2 / 3 + sp.Rational(33, 2) * L * g + sp.Rational(125, 4) * g**2)
            / 20736,
            "rational_two_loop_measure_lower": (16 * 3**2) ** 2 - 20736,
        },
        "bounds": {
            "actual_mass_above_thirty_two": M > 32,
            "actual_cubic_squared_is_small_exact_rational": g
            == sp.Rational(1, 67108864),
            "actual_L_and_g_strictly_positive": L > 0 and g > 0,
            "every_group_bound_positive": all(v > 0 for v in groups.values()),
            "actual_two_loop_quadratic_bound_below_one_e_minus_18": 0
            < B2
            < sp.Rational(1, 10**18),
            "actual_one_plus_two_loop_bound_below_one_e_minus_18": 0
            < combined
            < sp.Rational(1, 10**18),
            "strict_positive_unit_disc_factored_inverse_gap": 1 - combined > 0,
        },
    }


def bad_cases():
    invalid = (
        True,
        False,
        1.0,
        sp.Float(1),
        "1",
        None,
        sp.oo,
        -sp.oo,
        sp.zoo,
        sp.I,
        sp.nan,
        sp.Symbol("r"),
    )
    return [
        ("inexact_radius_" + str(i), point, (v,)) for i, v in enumerate(invalid)
    ] + [
        ("outside_unit_disc_" + str(i), point, (v,))
        for i, v in enumerate((0, -1, Fraction(3, 2), 2))
    ]
