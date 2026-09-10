"""Full-propagator radius-one bounds for the complete finite MS shift."""

from functools import cache

import sympy as s
from p8_vacuum_two_loop_finite_contact.calibration import rational


def enclosure(L, tree_b2, ell=1200, ellH=600, Q=144):
    L, tree_b2, ell, ellH, Q = map(rational, (L, tree_b2, ell, ellH, Q))
    if (
        L <= 0
        or tree_b2 <= 0
        or not 0 <= ell <= 1200
        or not 0 <= ellH <= 600
        or not 0 < Q <= 144
    ):
        raise ValueError("Need L,tree>0; 0<=ell<=1200; 0<=ellH<=600; 0<Q_lower<=144")
    old = (24 + 32 * ellH) * L**3 / Q**2
    linear = (24 + 16 * ellH) * ell * L**3 / Q**2
    square = 3 * ell**2 * L**2 * tree_b2 / (4 * Q**2)
    return {
        "old_complete_double_bubble_upper": old,
        "full_linear_scale_term_upper": linear,
        "quadratic_scale_b2_absolute_upper": square,
        "quadratic_scale_b2_relative_upper": square / tree_b2,
        "complete_MS_double_bubble_upper": old + linear + square,
        "complete_MS_double_bubble_relative_upper": (old + linear + square) / tree_b2,
    }


@cache
def data():
    L, g, M, ell, ellH, Q = s.symbols("L g M ell ellH Q", positive=True)
    core = 3 * s.Rational(1, 4) * 2 * (ell / Q) * (2 * L) ** 3 * (2 / Q)
    triangles = (
        3 * s.Rational(1, 4) * (ell / Q) * (2 * L) ** 2 * 16 * g * ellH / (Q * M)
    )
    return {
        "external_C_upper": 2 * L,
        "finite_outer_IR_upper": 2 / Q,
        "both_full_heavy_triangle_upper": 16 * g * ellH / (Q * M),
        "linear_scale_full_disc_upper": (24 + 16 * ellH) * ell * L**3 / Q**2,
        "checks": {
            "all_channel_core_linear_prefactor": s.expand(
                core - 24 * ell * L**3 / Q**2
            ),
            "both_triangle_linear_prefactor": s.expand(
                triangles.subs(g, L * M / 3) - 16 * ell * ellH * L**3 / Q**2
            ),
            "complete_linear_sum": s.expand(
                core
                + triangles.subs(g, L * M / 3)
                - (24 + 16 * ellH) * ell * L**3 / Q**2
            ),
            "actual_logarithm_prefactor": 1200 * (24 + 16 * 600) - 11548800,
            "same_old_prefactor": 24 + 32 * 600 - 19224,
            "same_loop_measure_lower": 144**2 - 20736,
        },
        "scope": "The inherited full complex-disc propagator and triangle bounds, not a heavy-mass expansion. Cauchy applies to the complete linear scale function. The quadratic scale piece uses its separately exact forward coefficient and retains the iterated heavy inverses.",
    }
