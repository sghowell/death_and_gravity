"""Independent exact finite contact integrals and their anchored primitives."""

from functools import cache

import sympy as sp


@cache
def data():
    y, M = sp.symbols("nonnegative_radial positive_heavy_mass_squared", positive=True)
    logdiff = sp.log(y + 1) - sp.log(y + M)
    P1 = M * logdiff / (M - 1) ** 2 + 1 / ((M - 1) * (y + 1))
    P2 = (
        (M + 1) * logdiff / (M - 1) ** 3
        + M / ((M - 1) ** 2 * (y + M))
        + 1 / ((M - 1) ** 2 * (y + 1))
    )
    J1 = M * sp.log(M) / (M - 1) ** 2 - 1 / (M - 1)
    J2 = (M + 1) * sp.log(M) / (M - 1) ** 3 - 2 / (M - 1) ** 2
    f1 = y / ((y + 1) ** 2 * (y + M))
    f2 = y / ((y + 1) ** 2 * (y + M) ** 2)
    return {
        "y": y,
        "M": M,
        "one_heavy_factor_integrand": f1,
        "two_heavy_factor_integrand": f2,
        "one_heavy_factor_primitive": P1,
        "two_heavy_factor_primitive": P2,
        "J1": J1,
        "J2": J2,
        "scope": "All-radius exact positive integrals at M>1, with the conservative estimates used only at actual M>32. The second primitive is checked directly, not only inferred by differentiation under an integral.",
        "checks": {
            "first_radial_primitive_derivative": sp.factor(sp.diff(P1, y) - f1),
            "first_radial_zero_anchor": sp.simplify(-P1.subs(y, 0) - J1),
            "first_radial_infinity_anchor": sp.limit(P1, y, sp.oo),
            "second_radial_primitive_derivative": sp.factor(sp.diff(P2, y) - f2),
            "second_radial_zero_anchor": sp.simplify(-P2.subs(y, 0) - J2),
            "second_radial_infinity_anchor": sp.limit(P2, y, sp.oo),
            "independent_mass_derivative_relation": sp.factor(J2 + sp.diff(J1, M)),
            "second_integrand_positive_upper_remainder": sp.factor(
                f1 / M - f2 - y**2 / (M * (y + 1) ** 2 * (y + M) ** 2)
            ),
            "mass_ratio_majorant_above_thirty_two_margin": 2 * (32 - 1) ** 2
            - 32**2
            - 898,
            "mass_ratio_majorant_polynomial": sp.expand(
                2 * (M - 1) ** 2 - M**2 - (M**2 - 4 * M + 2)
            ),
        },
    }
