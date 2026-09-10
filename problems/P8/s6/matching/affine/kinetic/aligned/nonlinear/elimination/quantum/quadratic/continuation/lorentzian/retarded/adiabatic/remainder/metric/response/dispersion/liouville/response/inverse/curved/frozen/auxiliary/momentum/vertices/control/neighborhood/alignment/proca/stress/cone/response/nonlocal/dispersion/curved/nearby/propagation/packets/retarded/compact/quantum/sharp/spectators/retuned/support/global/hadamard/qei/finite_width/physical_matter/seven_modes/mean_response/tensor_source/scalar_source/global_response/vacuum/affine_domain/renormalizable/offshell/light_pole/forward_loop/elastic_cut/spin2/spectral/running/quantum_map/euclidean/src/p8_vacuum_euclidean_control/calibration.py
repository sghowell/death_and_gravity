"""Exact all-momentum enclosures and restricted insertion-order bounds."""

from fractions import Fraction
from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model
from p8_vacuum_light_pole import kernel as parent_kernel

from . import dyson, map_domain


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational Euclidean invariant")
    return sp.Rational(value)


def point(invariant=0):
    y = rational(invariant)
    if y < 0:
        raise ValueError("Require a nonnegative Euclidean invariant")
    d = data()
    g, M = d["actual_cubic_squared"], d["actual_heavy_mass_squared"]
    bmin = sp.Rational(3, 16) / M
    low = g * bmin * bmin * (y + 1) / (1024 * (1 + bmin * (y + 1)))
    high = min(
        d["actual_uniform_relative_self_energy_upper"], g * (y + 1) / (864 * M * M)
    )
    return {
        "Euclidean_invariant": y,
        "relative_self_energy_strict_lower": low,
        "relative_self_energy_strict_upper": high,
        "displayed_inverse_relative_lower": 1 - high,
        "single_line_covariance_relative_upper": 1 / (1 - high),
        "scope": "An exact parameter-integral enclosure at this Euclidean invariant. Neither a timelike amplitude nor a full higher-loop self-energy.",
    }


def insertion_tail(order=1):
    if type(order) is not int or not 0 <= order <= 32:
        raise ValueError("Require a native retained insertion order zero through 32")
    return dyson.tails(data()["actual_uniform_relative_self_energy_upper"], order)


@cache
def data():
    p = model.data()["actual_parameters"]
    M, g = p["heavy_mass_squared"], p["cubic_coupling_squared"]
    alpha = g / (288 * M)
    tail = dyson.tails(alpha, 1)
    x = sp.Symbol("unit_Feynman_parameter", real=True)
    Ms = sp.Symbol("heavy_mass_squared_at_least_one", positive=True)
    b = x * (1 - x) / (x * Ms + (1 - x) ** 2)
    return {
        "actual_heavy_mass_squared": M,
        "actual_cubic_squared": g,
        "actual_uniform_relative_self_energy_upper": alpha,
        "actual_single_line_total_relative_shift_upper": alpha / (1 - alpha),
        "actual_two_line_total_relative_shift_upper": alpha
        * (2 - alpha)
        / (1 - alpha) ** 2,
        "actual_single_line_beyond_one_insertion_upper": tail["single_line_tail"],
        "actual_two_line_beyond_one_total_insertion_upper": tail[
            "two_line_total_order_tail"
        ],
        "actual_map_domain": {
            k: v for k, v in map_domain.data().items() if k not in {"checks", "bounds"}
        },
        "checks": {
            "same_fixed_parent_kinetic_majorant": alpha
            - parent_kernel.data()["actual_rational_finite_kinetic_counterterm_upper"],
            "parameter_denominator_below_M_positive_gap": sp.expand(
                Ms - (x * Ms + (1 - x) ** 2) - (1 - x) * (Ms - 1 + x)
            ),
            "parameter_numerator_lower_on_middle_half": sp.expand(
                x * (1 - x)
                - sp.Rational(3, 16)
                - (x - sp.Rational(1, 4)) * (sp.Rational(3, 4) - x)
            ),
            "parameter_weight_above_numerator_over_M": sp.factor(
                b
                - x * (1 - x) / Ms
                - x * (1 - x) ** 2 * (Ms - 1 + x) / (Ms * (x * Ms + (1 - x) ** 2))
            ),
            "positive_subinterval_prefactor": sp.Rational(1, 2)
            * sp.Rational(1, 2)
            / 256
            - sp.Rational(1, 1024),
        },
        "bounds": {
            "actual_mass_above_one_for_continuous_parameter_bounds": M > 1,
            "actual_uniform_relative_self_energy_below_three_e_minus_208": 0
            < alpha
            < sp.Rational(3, 10**208),
            "actual_single_line_total_relative_shift_below_three_e_minus_208": alpha
            / (1 - alpha)
            < sp.Rational(3, 10**208),
            "actual_two_line_total_relative_shift_below_six_e_minus_208": alpha
            * (2 - alpha)
            / (1 - alpha) ** 2
            < sp.Rational(6, 10**208),
            "actual_single_line_beyond_one_insertion_below_one_e_minus_414": 0
            < tail["single_line_tail"]
            < sp.Rational(1, 10**414),
            "actual_two_line_beyond_one_insertion_below_one_e_minus_414": 0
            < tail["two_line_total_order_tail"]
            < sp.Rational(1, 10**414),
            "uniform_relative_kernel_below_one_gives_positive_Euclidean_form": alpha
            < 1,
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
        sp.Symbol("y"),
    )
    rows = [
        ("inexact_Euclidean_invariant_" + str(i), point, (v,))
        for i, v in enumerate(invalid)
    ]
    rows += [
        ("negative_Euclidean_invariant_" + str(i), point, (v,))
        for i, v in enumerate((-1, Fraction(-1, 2)))
    ]
    rows += [
        ("invalid_insertion_order_" + str(i), insertion_tail, (v,))
        for i, v in enumerate((True, False, 1.0, sp.Integer(1), -1, 33, "1", None))
    ]
    rows += [
        ("invalid_multiplier_" + str(i), dyson.tails, (v, 1))
        for i, v in enumerate(
            (
                True,
                False,
                0,
                1,
                1.0,
                sp.Float(1),
                sp.Rational(0),
                sp.Rational(1),
                sp.Rational(-1, 2),
                sp.Rational(3, 2),
                sp.oo,
                sp.I,
                None,
            )
        )
    ]
    return rows
