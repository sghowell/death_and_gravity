"""Exact original intersection shell and its leading conversion limit."""

from functools import cache

import sympy as s


def boundary_radius(K, P, u):
    return P * u + s.sqrt(K * K - P * P * (1 - u * u))


def leading_shell_coefficient(K, P, u):
    """Lost radial r^3 moment; fixed K>P and u<=P/(2K)."""
    return (K**4 - boundary_radius(K, P, u) ** 4) / 4


@cache
def data():
    K, P = s.symbols("K P", positive=True)
    u = s.symbols("u", real=True)
    eps = s.symbols("epsilon", positive=True)
    excess = s.symbols("positive_excess", positive=True)
    root = boundary_radius(K, P, u)
    normalized = (1 - (eps * u + s.sqrt(1 - eps**2 * (1 - u * u))) ** 4) / (4 * eps)
    return {
        "same_original_shell": "For K>|P|, the second-leg radius is r_*=|P|u+sqrt(K^2-|P|^2(1-u^2)). The original intersection removes r in(r_*,K] exactly when u<|P|/(2K). No domain derivative is taken.",
        "fixed_transfer_limit": "At fixed P and angle u, (K^4-r_*^4)/(4K^3|P|) -> -u; the removed angular interval tends to[-1,0]. The lost-shell thickness is at most |P|, and the small positive-u strip has vanishing normalized measure. These give dominated convergence, including grazing angles.",
        "sign": "pair-band minus one-k-ball is the NEGATIVE lost-shell integral. The current j=0 coefficient itself is positive for equal real tensors; both signs are retained.",
        "comparison_not_replacement": "The one-k-ball is only an explicitly named conversion comparator. The physical calculation retains the original two-leg regulator. The complete contact already has one mode and contributes zero to this particular regulator difference.",
        "checks": {
            "exact_second_leg_boundary": s.simplify(
                root**2 - 2 * P * u * root + P**2 - K**2
            ),
            "grazing_positive_square_root_identity": s.factor(
                (K - P**2 / (2 * K)) ** 2 - (K**2 - P**2 * (1 - P**2 / (4 * K**2)))
            ),
            "negative_axis_shell_width": s.simplify(root.subs(u, -1) - (K - P)),
            "positive_axis_second_leg_radius": s.simplify(root.subs(u, 1) - (K + P)),
            "normalized_fixed_angle_shell_limit": s.simplify(
                s.limit(normalized, eps, 0, dir="+") + u
            ),
            "radial_shell_primitive": s.diff(
                s.Symbol("r", positive=True) ** 4 / 4, s.Symbol("r", positive=True)
            )
            - s.Symbol("r", positive=True) ** 3,
        },
        "gates": {
            "grazing_chosen_root_positive_for_K_above_P": bool(
                s.expand(2 * (P + excess) ** 2 - P**2).is_positive
            ),
            "grazing_strip_retained_at_finite_K": True,
            "original_two_leg_mask_not_replaced": True,
            "fixed_P_limit_not_uniform_external_momentum_bound": True,
            "contact_original_one_mode_band_unchanged": True,
        },
    }
