"""The fixed contact is affine, with a holomorphic dimensional lift."""

from functools import cache

import sympy as s
from p8_vacuum_two_loop_finite_contact import radial


@cache
def data():
    L, g, Q = s.symbols("L g Q", positive=True)
    d = radial.data()
    y, M, J1, J2 = (d[k] for k in ("y", "M", "J1", "J2"))
    k = 6 * g * J1 / Q
    c = 6 * g**2 * (J1 / M + J2) / Q
    sigma = 6 * ((-L + g / M) * g * J1 + g**2 * J2) / Q
    e, mu = s.symbols("epsilon mu", positive=True)
    measure = s.exp(s.EulerGamma * e) * mu ** (2 * e) / s.gamma(2 - e)
    j1D = measure * y ** (1 - e) / ((y + 1) ** 2 * (y + M))
    j2D = measure * y ** (1 - e) / ((y + 1) ** 2 * (y + M) ** 2)
    kernel = 6 * ((-L + g / M) * g * j1D + g**2 * j2D) / Q
    A = (L - g / M) / 2
    F = A - g / (y + M)
    split = (s.log(M + 1) + 1 / (M + 1)) / M
    return {
        "symbols": {"L": L, "g": g, "M": M, "Q": Q, "y": y, "epsilon": e, "mu": mu},
        "J1": J1,
        "J2": J2,
        "sigma": sigma,
        "positive_affine_k": k,
        "positive_affine_c": c,
        "dimensional_J1_integrand": j1D,
        "dimensional_J2_integrand": j2D,
        "dimensional_sigma_integrand": kernel,
        "J1_split_upper": split,
        "J1_log_upper": s.log(4 * M) / M,
        "analytic_strip": "-1 < Re epsilon < 2, M>1, mu>0",
        "checks": {
            "same_exact_affine_contact": s.expand(sigma + k * L - c),
            "fixed_parameter_contact_slope": s.diff(sigma, L) + k,
            "contact_second_L_derivative": s.diff(sigma, L, 2),
            "dimensional_kernel_affine": s.expand(
                kernel + 6 * g * j1D * L / Q - 6 * g**2 * (j1D / M + j2D) / Q
            ),
            "four_dimensional_J1": s.simplify(
                j1D.subs(e, 0) - d["one_heavy_factor_integrand"]
            ),
            "four_dimensional_J2": s.simplify(
                j2D.subs(e, 0) - d["two_heavy_factor_integrand"]
            ),
            "same_potential_kernel": s.factor(
                kernel.subs(e, 0) - 6 * y * (F**2 - A**2) / (Q * (y + 1) ** 2)
            ),
            "J1_small_radius_primitive": s.diff(s.log(y + 1) + 1 / (y + 1), y)
            - y / (y + 1) ** 2,
            "J1_small_radius_anchor": s.simplify(
                (s.log(M + 1) + 1 / (M + 1) - 1) / M + 1 / M - split
            ),
            "J1_large_radius_majorant_difference": s.factor(
                1 / y**2
                - d["one_heavy_factor_integrand"]
                - ((M + 2) * y**2 + (2 * M + 1) * y + M)
                / (y**2 * (y + 1) ** 2 * (y + M))
            ),
            "squared_potential_negative_kernel": s.factor(
                F**2 - A**2 + g * (F + A) / (y + M)
            ),
        },
        "scope": "The same finite constant-field contact. The analytic dimensional lift is retained in both legs of a pole product; the cancellation holds for every common analytic lift, so no evanescent convention is silently fixed by its four-dimensional value.",
    }
