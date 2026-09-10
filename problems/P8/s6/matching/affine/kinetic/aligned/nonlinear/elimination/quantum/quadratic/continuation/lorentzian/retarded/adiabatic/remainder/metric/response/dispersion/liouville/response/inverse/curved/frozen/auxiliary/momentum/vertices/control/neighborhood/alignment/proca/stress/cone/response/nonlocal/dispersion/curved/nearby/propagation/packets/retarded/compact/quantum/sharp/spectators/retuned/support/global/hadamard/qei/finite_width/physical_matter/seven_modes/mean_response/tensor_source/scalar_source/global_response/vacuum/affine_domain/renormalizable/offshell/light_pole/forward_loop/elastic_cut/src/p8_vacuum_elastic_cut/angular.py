"""Exact positive angular integral with an anchored real atanh branch."""

from functools import cache

import sympy as sp


@cache
def data():
    B, k, g = sp.symbols("positive_B positive_k positive_cubic_squared", positive=True)
    C, z = sp.symbols("real_C physical_cosine", real=True)
    A = C + 2 * g * B / (B * B - k * k * z * z)
    J1 = sp.atanh(k * z / B) / (B * k)
    J2 = z / (2 * B * B * (B * B - k * k * z * z)) + sp.atanh(k * z / B) / (
        2 * B**3 * k
    )
    primitive = C * C * z + 4 * C * g * B * J1 + 4 * g * g * B * B * J2
    average = (
        C * C
        + 4 * C * g * sp.atanh(k / B) / k
        + 2 * g * g / (B * B - k * k)
        + 2 * g * g * sp.atanh(k / B) / (B * k)
    )
    q = sp.Symbol("atanh_ratio_parameter", positive=True)
    ratio = sp.atanh(q) / q
    return {
        "B": B,
        "k": k,
        "g": g,
        "C": C,
        "cosine": z,
        "full_tree_vertex": A,
        "anchored_squared_vertex_primitive": primitive,
        "exact_even_squared_vertex_average": average,
        "threshold_average": (C + 2 * g / B) ** 2,
        "branch": "B=M+k>k>=0; atanh(kz/B) is the real anchored branch for 0<=z<=1. The k=0 threshold is assigned by its continuous limit, not evaluated as zero divided by zero.",
        "numerical_boundary": "The closed expression contains enormous cancellations for the actual parameters. Strict physical bounds use the positive amplitude enclosure, not floating evaluation of this difference.",
        "checks": {
            "first_anchored_inverse_denominator_primitive": sp.factor(
                sp.diff(J1, z) - 1 / (B * B - k * k * z * z)
            ),
            "second_anchored_inverse_square_primitive": sp.factor(
                sp.diff(J2, z) - 1 / (B * B - k * k * z * z) ** 2
            ),
            "full_squared_tree_vertex_primitive": sp.factor(
                sp.diff(primitive, z) - A * A
            ),
            "full_primitive_anchor_zero": primitive.subs(z, 0),
            "closed_exact_angular_average": sp.factor(primitive.subs(z, 1) - average),
            "threshold_atanh_ratio_limit": sp.limit(ratio, q, 0, dir="+") - 1,
            "threshold_constant_vertex_average": sp.expand(
                C * C + 4 * C * g / B + 4 * g * g / (B * B) - (C + 2 * g / B) ** 2
            ),
        },
    }
