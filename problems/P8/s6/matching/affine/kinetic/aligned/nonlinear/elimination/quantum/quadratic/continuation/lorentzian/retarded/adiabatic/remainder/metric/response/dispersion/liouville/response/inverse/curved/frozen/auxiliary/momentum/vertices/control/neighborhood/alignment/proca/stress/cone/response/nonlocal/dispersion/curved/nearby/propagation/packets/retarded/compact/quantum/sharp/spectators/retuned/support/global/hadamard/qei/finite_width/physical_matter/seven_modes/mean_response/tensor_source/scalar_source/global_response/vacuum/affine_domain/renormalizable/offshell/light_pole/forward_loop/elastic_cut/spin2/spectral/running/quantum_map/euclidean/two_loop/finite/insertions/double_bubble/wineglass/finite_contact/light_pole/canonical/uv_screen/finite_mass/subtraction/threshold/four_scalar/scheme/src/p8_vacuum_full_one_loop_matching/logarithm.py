"""Exact rational enclosure of the interaction-scale logarithm."""

from functools import cache

import sympy as sp
from p8_vacuum_gauge_light_cut_subtraction import logarithm as log_two


def enclosure(terms=8):
    if type(terms) is not int:
        raise TypeError("Require a native integer logarithm truncation")
    if not 1 <= terms <= 64:
        raise ValueError("Require one to 64 retained logarithm terms")
    r = sp.Rational(1, 9)
    low = 2 * sum(r ** (2 * j + 1) / (2 * j + 1) for j in range(terms))
    tail = 2 * r ** (2 * terms + 1) / ((2 * terms + 1) * (1 - r * r))
    two = log_two.enclosure(terms)
    log10lo = 3 * two["lower"] + low
    log10hi = 3 * two["upper"] + low + tail
    return {
        "terms": terms,
        "log_five_fourths_lower": low,
        "log_five_fourths_upper": low + tail,
        "log_ten_lower": log10lo,
        "log_ten_upper": log10hi,
        "scale_log_lower": 400 * log10lo,
        "scale_log_upper": 400 * log10hi,
    }


@cache
def data():
    d = enclosure()
    x = sp.symbols("atanh_ratio", real=True)
    checks = {
        "five_fourths_atanh_endpoint": (1 + sp.Rational(1, 9)) / (1 - sp.Rational(1, 9))
        - sp.Rational(5, 4),
        "positive_real_log_ten_factorization": 2**3 * sp.Rational(5, 4) - 10,
        "reference_mass_squared_logarithm_exponent": 2 * 200 - 400,
        "atanh_log_derivative": sp.factor(
            sp.diff(sp.log((1 + x) / (1 - x)), x) - 2 / (1 - x * x)
        ),
        "atanh_log_zero_anchor": sp.log((1 + x) / (1 - x)).subs(x, 0),
        "scale_log_enclosure_width": d["scale_log_upper"]
        - d["scale_log_lower"]
        - 400 * (d["log_ten_upper"] - d["log_ten_lower"]),
    }
    return {
        "exact_scale_log": "ell=log(mF^2/m_phi^2)=400log(10), with mF=10^200 and physical m_phi=1",
        "rational_enclosure": d,
        "strict_coarse_scale_log_lower": 900,
        "strict_coarse_scale_log_upper": 1000,
        "tail_proof": "log(10)=3log(2)+log(5/4), with positive atanh ratios 1/3 and 1/9. Every omitted odd denominator is bounded below by the first omitted one, then the geometric tail is summed. The complete logarithm, not a free order-one coefficient, is retained.",
        "checks": checks,
    }
