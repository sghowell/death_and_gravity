"""Exact rational enclosure of log(2) from its positive atanh series."""

from functools import cache

import sympy as sp


def enclosure(terms=8):
    if type(terms) is not int:
        raise TypeError("Require a native integer series truncation")
    if not 1 <= terms <= 64:
        raise ValueError("Require between one and 64 retained terms")
    r = sp.Rational(1, 3)
    low = 2 * sum(r ** (2 * j + 1) / (2 * j + 1) for j in range(terms))
    tail = 2 * r ** (2 * terms + 1) / ((2 * terms + 1) * (1 - r * r))
    return {"terms": terms, "lower": low, "upper": low + tail, "tail_upper": tail}


@cache
def data():
    x = sp.symbols("atanh_argument", real=True)
    n = sp.symbols("positive_series_count", integer=True, positive=True)
    d = enclosure()
    checks = {
        "log_two_from_atanh_endpoint": sp.Rational(
            1 + sp.Rational(1, 3), 1 - sp.Rational(1, 3)
        )
        - 2,
        "log_ratio_derivative": sp.factor(
            sp.diff(sp.log((1 + x) / (1 - x)), x) - 2 / (1 - x * x)
        ),
        "log_ratio_zero_anchor": sp.log((1 + x) / (1 - x)).subs(x, 0),
        "tail_geometric_majorant": 2 * x ** (2 * n + 1) / ((2 * n + 1) * (1 - x * x))
        - 2 * x * x ** (2 * n) / ((2 * n + 1) * (1 - x * x)),
        "interval_width": d["upper"] - d["lower"] - d["tail_upper"],
    }
    for j in range(8):
        checks[f"Taylor_coefficient_{j}"] = sp.diff(
            sp.log((1 + x) / (1 - x)), x, 2 * j + 1
        ).subs(x, 0) / sp.factorial(2 * j + 1) - sp.Rational(2, 2 * j + 1)
    return {
        "log_two_enclosure": d,
        "positive_series": "log(2)=2 sum_{j>=0} (1/3)^(2j+1)/(2j+1)",
        "strict_tail_proof": "Every omitted denominator is at least 2N+1; replacing them all by this value and summing a geometric series gives the stated strict positive tail upper bound.",
        "checks": checks,
    }


def bad_cases():
    return [
        (f"invalid_log_truncation_{j}", enclosure, (value,))
        for j, value in enumerate(
            (True, False, 8.0, sp.Integer(8), "8", None, 0, -1, 65)
        )
    ]
