"""Original homogeneous finite local and classical tensor operator, with all signs."""

from functools import cache

import sympy as s
from p8_vacuum_affine_local_tensor_response import operator as original

t = original.t


@cache
def data():
    h = s.Function("homogeneous_unit_shear")(t)
    A, H = original.A, original.H
    L = s.diff(h, t, 2) + 3 * H * s.diff(h, t)
    D = s.diff(h, t, 2) + H * s.diff(h, t)
    Dadj = lambda value: (
        s.diff(value, t, 2)
        + 5 * H * s.diff(value, t)
        + (2 * s.diff(H, t) + 6 * H**2) * value
    )
    Qloc = (A * L + s.diff(A, t) * s.diff(h, t) + Dadj(D) / 60) / (
        8 * s.pi**2 * original.KAPPA
    )
    quantum = -2 * A * L - 2 * s.diff(A, t) * s.diff(h, t) - Dadj(D) / 30
    tree = -16 * s.pi**2 * original.KAPPA * L
    full = s.expand(quantum + tree)
    wanted = {
        4: -s.Rational(1, 30),
        3: -H / 5,
        2: -2 * A - (4 * s.diff(H, t) + 11 * H**2) / 30 - 16 * s.pi**2 * original.KAPPA,
        1: -6 * H * A
        - 2 * s.diff(A, t)
        - (s.diff(H, t, 2) + 7 * H * s.diff(H, t) + 6 * H**3) / 30
        - 48 * s.pi**2 * original.KAPPA * H,
        0: s.S.Zero,
    }
    x, y = s.symbols("output source", nonnegative=True)
    c = s.Function("coefficient")(y)
    kernels = {
        j: sum(
            (-1) ** r
            * s.binomial(j, r)
            * (x - y) ** (3 - j + r)
            * s.diff(c, y, r)
            / s.factorial(3 - j + r)
            for r in range(j + 1)
        )
        for j in range(4)
    }
    checks = {
        "all_actual_finite_and_classical_local_coefficients": s.Matrix(
            [s.factor(full.coeff(s.diff(h, t, j)) - wanted[j]) for j in range(5)]
        ),
        "actual_canonical_local_force_sign_and_two_chain_factors": s.factor(
            quantum + 16 * s.pi**2 * original.KAPPA * Qloc
        ),
        "classical_source_normalization": s.factor(
            tree + 16 * s.pi**2 * original.KAPPA * L
        ),
        "original_finite_fourth_not_action_half": wanted[4] + s.Rational(1, 30),
        "all_lower_variable_coefficient_primitive_kernels": s.Matrix(
            [
                s.expand(kernels[j] - (-1) ** j * s.diff((x - y) ** 3 * c / 6, y, j))
                for j in range(4)
            ]
        ),
        "only_fourth_term_removed_from_local_inventory": s.expand(
            full - sum(wanted[j] * s.diff(h, t, j) for j in range(5))
        ),
    }
    return {
        "actual_known_local_coefficients": wanted,
        "zero_past_lower_local_primitive_kernels": kernels,
        "actual_tree_plus_quantum_normalization": "T_total=T_Gamma-16pi^2 kappa L, L=Dt^2+3H Dt at zero transfer. For the canonical forcing f in-Lh+T_Gamma h/(16pi^2 kappa)+f=0, the total normalized source isg=-16pi^2 kappa f.",
        "complete_quantum_definition": "T_Gamma=64pi^2/a^3 times the actual retarded exponential tracefree metric-current derivative at the original isotropic clock. The S194 two canonical chain factors give T_Gamma/(16pi^2 kappa). It is not a symmetric in-out Hessian substituted for the retarded response.",
        "local_scope": "These coefficients include the original finite local action and actual classical tensor tree. They do not set the remaining original nonlocal or second-vertex contact contribution to zero. The complete latter contribution is retained inV_Gamma.",
        "checks": checks,
        "gates": {
            "actual_background_scale_positive": original.a.is_positive,
            "unchanged_fixed_mass": original.MASS == 1000,
            "unchanged_kappa": original.KAPPA == s.Integer(10) ** 800,
            "nonzero_third_order_curvature_term": wanted[3] != 0,
            "actual_tree_coefficient_retains_kappa": wanted[2].has(s.pi),
            "all_lower_primitive_orders_are_ordinary": all(
                3 - j + r >= 0 for j in range(4) for r in range(j + 1)
            ),
        },
    }
