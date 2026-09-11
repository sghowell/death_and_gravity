"""Strict exact-rational assembly without numerical epsilon or scale cutoffs."""

import sympy as s


def nonnegative(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer, s.Rational)):
        raise TypeError("An exact real rational bound is required")
    value = s.Rational(value)
    if value < 0:
        raise ValueError("A bound must be nonnegative")
    return value


def first_variation_upper(k, scalar, fermion):
    k, scalar, fermion = map(nonnegative, (k, scalar, fermion))
    return k * (4 * scalar + 2 * fermion)


def assemble(tree, first, raw, cross, field_relative):
    tree, first, raw, cross, field_relative = map(
        nonnegative, (tree, first, raw, cross, field_relative)
    )
    if tree == 0:
        raise ValueError("The reference tree coefficient must be positive")
    second = raw + cross + tree * field_relative
    return {
        "second_absolute": second,
        "second_relative": second / tree,
        "total_absolute": first + second,
        "total_relative": (first + second) / tree,
        "formal_uniform_lower": tree - first - second,
    }


def data():
    T, E1, R, C, F = s.symbols("T E1 R C F", positive=True)
    h = s.symbols("h", nonnegative=True)
    second = R + C + T * F
    return {
        "second_absolute_formula": second,
        "formal_coefficient_lower_formula": T - E1 - second,
        "homotopy_error_majorant": h * E1 + h**2 * second,
        "checks": {
            "first_plus_second_sum": E1 + second - (E1 + R + C + T * F),
            "uniform_homotopy_difference_factored": s.expand(
                E1 + second - h * E1 - h**2 * second - (1 - h) * (E1 + (1 + h) * second)
            ),
            "second_relative_split": s.expand(second / T - R / T - C / T - F),
            "formal_positive_reference_margin": (T - E1 - second) + (E1 + second) - T,
        },
        "scope": "Absolute bounds on coefficients of a formal two-loop polynomial for 0<=h<=1. This is not a bound on the physical higher-loop remainder.",
    }
