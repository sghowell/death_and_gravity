"""Independent elementary ordered integrals behind the sector bound."""

from functools import cache
from itertools import product
from math import factorial

import sympy as sp

LARGE_PATTERNS = (
    (),
    (-2,),
    (-2, -2),
    (-2, 0),
    (-2, -2, 0),
    (-2, 0, -2),
    (-2, -2, 0, 0),
    (-2, 0, -2, 0),
)


def large_data(powers):
    if (
        type(powers) is not tuple
        or any(type(p) is not int for p in powers)
        or powers not in LARGE_PATTERNS
    ):
        raise ValueError("Unsupported large ordered power pattern")
    return _large_data(powers)


@cache
def _large_exact(powers):
    B = sp.Symbol("heavy_cutoff_at_least_one", positive=True)
    value = sp.Integer(1)
    for p in reversed(powers):
        primitive = sp.integrate(B**p * value, B)
        value = sp.expand(primitive - primitive.subs(B, 1))
    return value


@cache
def _large_data(powers):
    B = sp.Symbol("heavy_cutoff_at_least_one", positive=True)
    value = _large_exact(powers)
    prefix = tuple(sum(p + 1 for p in powers[:j]) for j in range(1, len(powers) + 1))
    zeros = prefix.count(0)
    denominator = factorial(zeros)
    for a in prefix:
        if a < 0:
            denominator *= -a
    z = sp.symbols("spacing0:" + str(len(powers)), nonnegative=True)
    coordinates = [sum(z[j:]) for j in range(len(powers))]
    return {
        "powers": powers,
        "B": B,
        "exact_ordered_integral": value,
        "positive_log_upper": sp.log(B) ** zeros / denominator,
        "log_power": zeros,
        "checks": {
            "anchored_ordered_derivative": sp.simplify(
                sp.diff(value, B)
                - (B ** powers[0] * _large_exact(powers[1:]) if powers else 0)
            ),
            "ordered_integral_at_unit_endpoint": sp.simplify(
                value.subs(B, 1) - (1 if not powers else 0)
            ),
            "log_spacing_exponent_identity": sp.expand(
                sum((p + 1) * u for p, u in zip(powers, coordinates))
                - sum(a * x for a, x in zip(prefix, z))
            ),
        },
    }


@cache
def small_patterns():
    return tuple(
        p
        for n in range(1, 8)
        for p in product((-2, 0), repeat=n)
        if p.count(-2) <= 2 and all(sum(x + 1 for x in p[j:]) > 0 for j in range(n))
    )


def small_data(powers):
    if (
        type(powers) is not tuple
        or any(type(p) is not int for p in powers)
        or powers not in small_patterns()
    ):
        raise ValueError("Nonintegrable or unsupported small ordered pattern")
    return _small_data(powers)


@cache
def _small_data(powers):
    b = sp.Symbol("positive_small_upper", positive=True)
    value = sp.Integer(1)
    for p in reversed(powers):
        value = sp.integrate(b**p * value, (b, 0, b))
    tails = tuple(sum(x + 1 for x in powers[j:]) for j in range(len(powers)))
    formula = b ** tails[0] / sp.prod(tails)
    return {
        "powers": powers,
        "exact_small_ordered_integral": value,
        "positive_tail_exponents": tails,
        "closed_formula": formula,
        "checks": {
            "independent_nested_integral_matches_tail_product": sp.simplify(
                value - formula
            )
        },
    }


@cache
def data():
    # Include suffixes for derivative checks, even when not a leading large pattern.
    rows = []
    checks = {}
    for n, p in enumerate(LARGE_PATTERNS):
        d = large_data(p)
        rows.append({k: v for k, v in d.items() if k != "checks"})
        checks.update({"large_" + str(n) + "_" + k: v for k, v in d["checks"].items()})
    small = []
    for n, p in enumerate(small_patterns()):
        d = small_data(p)
        small.append({k: v for k, v in d.items() if k != "checks"})
        checks.update({"small_" + str(n) + "_" + k: v for k, v in d["checks"].items()})
    return {
        "large_ordered_integrals": tuple(rows),
        "small_ordered_integrals": tuple(small),
        "scope": "Elementary exact nested integrals and logarithmic-coordinate bounds, not exact values of full graph integrals.",
        "checks": checks,
    }
