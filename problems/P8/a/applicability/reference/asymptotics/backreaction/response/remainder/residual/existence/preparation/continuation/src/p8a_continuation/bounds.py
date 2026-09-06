"""Rational pole and resolvent bounds; no numerical Lambert-W assumptions."""

from functools import cache

import sympy as sp
from p8a_preparation import bounds as prior

from .resolvent import named_parameters, positive


def integer_duration_norm_lower(exponent):
    """Lower norm at L=exponent/10^6 using e>8/3.

This is a frozen-operator statement, not an existence domain for A.11.
The returned bound may be negative for very short intervals.
    """
    exponent = positive(exponent)
    if exponent.q != 1:
        raise ValueError("the exponent must be a positive integer")
    return ((sp.Rational(8, 3)**exponent-1)/17-sp.Rational(8, 3))


@cache
def calibration():
    params = named_parameters(delta=prior.calibration()["delta"])
    c = params["c"]
    beta_min, beta_max = -sp.Rational(53, 60), sp.Rational(11, 30)
    lower, upper = sp.Integer(10**6), sp.Integer(2*10**6)
    low_den_upper = lower**2*(6*sp.Rational(5, 2)+beta_max)-c
    high_den_lower = upper**2*(sp.Rational(1, 2)+6*2+beta_min)-c
    w_bound = 2*(1+6*sp.Rational(5, 2)+beta_max)
    weighted_den = high_den_lower/upper**2
    weighted_cap = 2/weighted_den
    at10, at100 = integer_duration_norm_lower(10), integer_duration_norm_lower(100)
    margins = {
        "lower_pole_test": -low_den_upper,
        "upper_pole_test": high_den_lower,
        "W0_less_than_33": 33-w_bound,
        "weighted_denominator": weighted_den,
        "ten_micro_norm_above_1000": at10-1000,
        "hundred_micro_norm_above_10_to_28": at100-10**28,
        "e_greater_than_8_over_3": sum(sp.Rational(1, sp.factorial(j)) for j in range(5))-sp.Rational(8, 3),
        "exp_five_halves_greater_than_10": sum(sp.Rational(5, 2)**j/sp.factorial(j) for j in range(5))-10,
    }
    if any(value <= 0 for value in margins.values()):
        raise ValueError("a rational frozen-resolvent margin failed")
    return {"scale": params["scale"], "delta": params["delta"], "c": c,
            "beta_lower": beta_min, "beta_upper": beta_max,
            "positive_pole_lower": lower, "positive_pole_upper": upper,
            "isolated_log_pole_upper": sp.Integer(3),
            "W0_upper": w_bound, "rounded_W0_upper": sp.Integer(33),
            "positive_step_residue_lower": sp.Rational(1, 17),
            "damped_pair_step_absolute_cap": sp.Rational(8, 3),
            "weighted_sigma": upper, "weighted_norm_upper": weighted_cap,
            "norm_lower_at_10_to_minus_5": at10,
            "norm_lower_at_10_to_minus_4": at100,
            "A11_existing_length": prior.calibration()["length"],
            "strict_margins": margins}
