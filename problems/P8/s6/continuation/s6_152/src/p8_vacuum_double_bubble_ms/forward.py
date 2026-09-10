"""The exact forward coefficient of the quadratic scale term."""

from functools import cache

import sympy as s


@cache
def data():
    L, g, M, nu, ell, Q, h, deltaM = s.symbols("L g M nu ell Q h delta_M")
    D = M - 2
    Cplus = -L + g / (D - nu)
    Cminus = -L + g / (D + nu)
    Czero = -L + g / M
    cube = Cplus**3 + Cminus**3 + Czero**3
    b2cube = s.diff(cube, nu, 2).subs(nu, 0) / 2
    b2tree = 2 * g / D**3
    relative = 3 * (-L + g / D) * (-L + 2 * g / D)
    r = s.symbols("r", nonnegative=True)
    local = (ell / Q) ** 2 * b2cube / 4
    heavy_product = s.diff(g / (M + h * deltaM - 2), h, 2).subs(h, 0) / 2
    return {
        "full_three_channel_C_cubed": cube,
        "C_cubed_forward_coefficient": b2cube,
        "C_cubed_relative_to_tree": relative,
        "quadratic_scale_b2": local,
        "quadratic_scale_relative_upper": 3 * ell**2 * L**2 / (4 * Q**2),
        "checks": {
            "exact_three_channel_second_coefficient": s.factor(
                b2cube / b2tree - relative
            ),
            "relative_factorization": s.expand(
                relative - 3 * L**2 + 9 * L * g / D - 6 * g**2 / D**2
            ),
            "exact_quadratic_scale_relative": s.factor(
                local / b2tree - ell**2 * relative / (4 * Q**2)
            ),
            "t_channel_constant_no_b2": s.diff(Czero**3, nu, 2),
            "unit_ratio_positive_product": s.expand(
                (1 - r) * (1 - 2 * r) - (1 - 3 * r + 2 * r**2)
            ),
            "upper_product_margin": s.expand(
                1 - (1 - r) * (1 - 2 * r) - r * (3 - 2 * r)
            ),
            "ordinary_disjoint_heavy_mass_product": s.factor(
                heavy_product - g * deltaM**2 / D**3
            ),
            "C_cubed_keeps_heavy_cubic_inverse": s.expand((-L + g / D) ** 3).coeff(g, 3)
            - 1 / D**3,
        },
        "scope": "For the actual 0<g/(M-2)<L/2, the dimensionless product (1-r)(1-2r) lies in (0,1), so |b2[C^3 sum]/b2tree|<=3L^2. The h^3 heavy inverse belongs to ordinary repeated heavy counterterm insertions. This fixed contribution is not a new tunable nonlocal interaction.",
    }
