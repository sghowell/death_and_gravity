"""Exact solution and composition of the explicitly truncated one-loop flow."""

from functools import cache

import sympy as sp


@cache
def data():
    L0, g0, M0, q, q2 = sp.symbols(
        "positive_initial_quartic positive_initial_cubic_squared positive_initial_heavy_mass_squared positive_flow_root positive_second_flow_root",
        positive=True,
    )
    ell = sp.Symbol("nonnegative_logarithmic_reference_time", nonnegative=True)
    L = L0 / q**3
    g = g0 / q**2
    M = M0 + g0 / L0 * (1 - q)
    qdot = -L0 / q**2
    margin = L - 3 * g / M
    initial_margin = L0 - 3 * g0 / M0
    invariant_numerator = (L0 * M0 - 3 * g0) + 4 * g0 * (1 - q)
    Mtwice = M + (g / L) * (1 - q2)
    return {
        "q": q,
        "L0": L0,
        "g0": g0,
        "M0": M0,
        "ell": ell,
        "positive_flow_root_equation": sp.Eq(q**3, 1 - 3 * L0 * ell),
        "running_polynomial_quartic": L,
        "running_cubic_squared": g,
        "running_heavy_mass_squared": M,
        "running_completed_square_margin": margin,
        "positive_margin_numerator": invariant_numerator,
        "formal_one_loop_flow_boundary": 1 / (3 * L0),
        "scope": "Exact solution of the displayed one-loop beta ODE only. Its formal singularity is not a proven quantum Landau pole or an exclusion of a UV completion.",
        "checks": {
            "root_differential_equation": 3 * q * q * qdot + 3 * L0,
            "quartic_one_loop_flow": sp.factor(sp.diff(L, q) * qdot - 3 * L * L),
            "cubic_squared_one_loop_flow": sp.factor(sp.diff(g, q) * qdot - 2 * L * g),
            "heavy_mass_squared_one_loop_flow": sp.factor(sp.diff(M, q) * qdot - g),
            "quartic_initial_value": L.subs(q, 1) - L0,
            "cubic_initial_value": g.subs(q, 1) - g0,
            "heavy_mass_initial_value": M.subs(q, 1) - M0,
            "positive_completed_square_margin_numerator": sp.factor(
                (L * M - 3 * g) * q**3 - invariant_numerator
            ),
            "margin_initial_value": sp.factor(margin.subs(q, 1) - initial_margin),
            "quartic_flow_composition": sp.factor(L / q2**3 - L0 / (q * q2) ** 3),
            "cubic_squared_flow_composition": sp.factor(g / q2**2 - g0 / (q * q2) ** 2),
            "heavy_mass_flow_composition": sp.factor(
                Mtwice - (M0 + g0 / L0 * (1 - q * q2))
            ),
        },
    }
