"""Independent Fraction-only replay of the frozen-operator rational gate."""

from fractions import Fraction as Q
from math import factorial


def serialize(value):
    if isinstance(value, dict):
        return {name: serialize(item) for name, item in value.items()}
    return str(value)


def replay(a11):
    full = a11["derived_constants"]["full_map"]
    delta, length = Q(full["delta"]), Q(full["length"])
    scale = Q(5, 2)
    c = scale*scale/(30*delta)
    # Euler's constant in (0,1), log(5/4) in (0,1/4).
    beta_lower, beta_upper = -Q(19, 30)-Q(1, 4), 1-Q(19, 30)
    lower, upper = Q(10**6), Q(2*10**6)
    low_log_upper = 6*Q(5, 2)
    high_log_lower = Q(1, 2)+6*2
    low_sign_upper = lower**2*(low_log_upper+beta_upper)-c
    high_sign_lower = upper**2*(high_log_lower+beta_lower)-c
    w_upper = 2*(1+low_log_upper+beta_upper)
    weight_den = high_sign_lower/upper**2
    # Use e>8/3 in both integer-exponent duration calibrations.
    at10 = (Q(8, 3)**10-1)/17-Q(8, 3)
    at100 = (Q(8, 3)**100-1)/17-Q(8, 3)
    margins = {
        "lower_pole_test": -low_sign_upper,
        "upper_pole_test": high_sign_lower,
        "W0_less_than_33": 33-w_upper,
        "weighted_denominator": weight_den,
        "ten_micro_norm_above_1000": at10-1000,
        "hundred_micro_norm_above_10_to_28": at100-10**28,
        "e_greater_than_8_over_3": sum(Q(1, factorial(j)) for j in range(5))-Q(8, 3),
        "exp_five_halves_greater_than_10": sum(Q(5, 2)**j/factorial(j) for j in range(5))-10,
    }
    if any(value <= 0 for value in margins.values()):
        raise ValueError("an independent frozen-resolvent rational margin failed")
    return serialize({"scale": scale, "delta": delta, "c": c,
                      "beta_lower": beta_lower, "beta_upper": beta_upper,
                      "positive_pole_lower": lower, "positive_pole_upper": upper,
                      "isolated_log_pole_upper": Q(3), "W0_upper": w_upper,
                      "rounded_W0_upper": Q(33), "positive_step_residue_lower": Q(1, 17),
                      "damped_pair_step_absolute_cap": Q(8, 3),
                      "weighted_sigma": upper, "weighted_norm_upper": 2/weight_den,
                      "norm_lower_at_10_to_minus_5": at10,
                      "norm_lower_at_10_to_minus_4": at100,
                      "A11_existing_length": length, "strict_margins": margins})
