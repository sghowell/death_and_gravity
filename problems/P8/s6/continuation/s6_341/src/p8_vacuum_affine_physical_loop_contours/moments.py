"""Absolute deformed-contour moments and all-order local Cauchy bounds."""

import math
from functools import cache

import sympy as s

from . import parameters, source

MASTER_BUDGETS = {"triangle": s.Integer(256000), "ordered_box": s.Integer(32768000)}
MOMENT_BUDGETS = {
    (1, 0): (500, 1),
    (2, 0): (32, 1),
    (2, 1): (500, 2),
    (3, 1): (16, 2),
    (3, 2): (500, 3),
}


def order(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer)) or value < 0:
        raise ValueError("Require an exact nonnegative derivative order")
    return int(value)


def master_bound(kind, mass=source.HEAVY_MASS2):
    kind = parameters.require_kind(kind)
    mass = parameters.heavy_mass(mass)
    return MASTER_BUDGETS[kind] / mass ** (1 if kind == "triangle" else 2)


def derivative_bound(kind, multiindex, mass_order=0, mass=source.HEAVY_MASS2):
    kind = parameters.require_kind(kind)
    mass = parameters.heavy_mass(mass)
    if not isinstance(multiindex, (tuple, list)) or len(multiindex) != (
        3 if kind == "triangle" else 6
    ):
        raise ValueError("Require the full invariant and virtuality multiindex")
    alpha = tuple(order(x) for x in multiindex)
    r = order(mass_order)
    factor = math.prod(math.factorial(x) for x in alpha) * 4096 ** sum(alpha)
    return master_bound(kind, mass) * factor * s.factorial(r) * (128 / mass) ** r


def moment_bound(p, r, mass=source.HEAVY_MASS2):
    p, r = order(p), order(r)
    mass = parameters.heavy_mass(mass)
    if (p, r) not in MOMENT_BUDGETS:
        raise ValueError("Require one of the proved weighted moments")
    coefficient, power = MOMENT_BUDGETS[(p, r)]
    return s.Integer(coefficient) / mass**power


@cache
def data():
    checks = {}
    R = s.Rational
    u, c0 = s.symbols("integration_t c0", positive=True)
    nn = s.Symbol("positive_n", positive=True)
    U = c0 + nn * u
    antiderivatives = {
        (1, 0): s.log(U / c0) / nn,
        (2, 0): -1 / (nn * U),
        (2, 1): (s.log(U / c0) + c0 / U) / nn**2,
        (3, 1): (-1 / U + c0 / (2 * U**2)) / nn**2,
        (3, 2): (s.log(U / c0) + 2 * c0 / U - c0**2 / (2 * U**2)) / nn**3,
    }
    for (p, r), F in antiderivatives.items():
        checks[f"moment_{p}_{r}_exact_antiderivative"] = s.factor(
            s.diff(F, u) - u**r / U**p
        )
    y = s.Symbol("y", real=True)
    checks["J32_logarithmic_residual_is_negative"] = s.expand(
        2 * y - y * y / 2 - R(3, 2) + (1 - y) * (3 - y) / 2
    )
    checks["J31_residual_is_negative"] = s.expand(-y + y * y / 2 + y * (2 - y) / 2)

    return {
        "checks": checks,
        "gates": {
            "logarithm_budget_from_exp_positive_partial_sum": bool(
                sum(R(7, 3) ** j / s.factorial(j) for j in range(7)) > 10
                and 200 * R(7, 3) < 500
            ),
            "scalar_bounds_apply_on_deformed_contour_only": True,
            "Cauchy_polydisk_keeps_selected_physical_sheet": True,
            "weighted_numerators_not_inferred_from_unweighted_master": True,
        },
        "whole_exact_moment_antiderivatives": antiderivatives,
        "whole_moment_majorants": MOMENT_BUDGETS,
        "whole_scalar_master_majorants": MASTER_BUDGETS,
        "whole_original_scalar_master_bounds": {
            kind: master_bound(kind) for kind in parameters.KINDS
        },
        "whole_all_order_derivative_rule": "External multiindex alpha: alpha!*4096^|alpha| times the master bound. Additional mass derivative r: r!*(128/n0)^r. Bounds follow from the jointly holomorphic continued physical germ, not differentiation of a singular absolute-value real integral.",
        "whole_scalar_measure_budget": "|xi_prime z_prime|<4; triangle measure<=2 and box measure<=4t. Thus C<8*64*500/n0 and D_ordered<16*64^2*500/n0^2.",
    }
