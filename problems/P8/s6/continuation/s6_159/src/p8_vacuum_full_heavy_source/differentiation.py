"""Mass-differentiated sunset sectors with a threshold split, not an M loss."""

from functools import cache

import sympy as s


@cache
def data():
    a, M, r, t, e = s.symbols("a M r t epsilon", positive=True)
    V = 1 + t + r * t
    Bs = (M + a * r + a * r * t, a + M * r + a * r * t, a + a * r + M * r * t)
    Cs = (r + r * t, 1 + r * t, 1 + r)
    rows = {}
    for j, (B, C) in enumerate(zip(Bs, Cs)):
        F = V ** (e - 2) * B ** (1 - 2 * e)
        Fa = (1 - 2 * e) * C * V ** (e - 2) * B ** (-2 * e)
        rows[f"light_mass_derivative_sector_{j}"] = s.simplify(s.diff(F, a) - Fa)
        rows[f"mass_derivative_coefficient_{j}"] = s.diff(B, a) - C
    rho = s.Rational(1, 16)
    derivative_prefactor = s.Rational(9, 8) * (2 + 2 * (2 + rho) + 8 * rho * M)
    small = s.Rational(10, 3) / (1 - rho)
    large = 5 / rho
    rows.update(
        {
            "two_light_label_permutations_per_heavy_position": 2 * len(Bs) - 6,
            "uniform_Fa_difference_majorant": 5
            - 2 * s.Rational(9, 4)
            - s.Rational(1, 2),
            "r_derivative_majorant_gap": 10 * M
            - derivative_prefactor
            - s.Rational(151, 16) * (M - 2)
            - s.Rational(767, 64),
            "small_r_integral_multiplier": small - s.Rational(32, 9),
            "large_r_integral_multiplier": large - 80,
            "per_sector_below_84": 84 - small - large - s.Rational(4, 9),
            "six_sectors_plus_corners_below_700": 700 - 6 * 84 - 108 - 88,
            "corner_derivative_bound": 6 * s.Rational(9, 8) / rho - 108,
            "full_sunset_derivative_multiplier": 40 * 700 - 28000,
            "normalized_loop_insertion_has_half_derivative": s.Rational(1, 2) * 28000
            - 14000,
        }
    )
    return {
        "three_heavy_positions_each_twice": [{"B": B, "C": C} for B, C in zip(Bs, Cs)],
        "regulated_mass_derivative_integrand": "(1-2epsilon) C V^(epsilon-2) B^(-2epsilon), evaluated at a=1. C is the sum of the two light coefficients, so 0<=C<=2 and 0<=partial_r C<=2.",
        "corner_difference_bound": "|Fa(r)-Fa(0)|<=min(10 M r,5)(3M)^(1/8), on |epsilon|=1/16.",
        "split_point": 1 / (3 * M),
        "continued_mass_derivative_bound": "|partial_a S_D(a,a,M)|_(a=1)<28000 mF^(1/4)(3M)^(3/16)/Q^2.",
        "checks": rows,
        "scope": "Differentiation is at fixed common regulator and fixed other masses/scale. The compact corner-subtracted integrand is differentiated before the finite Laurent coefficient. Splitting r at (3M)^-1 avoids replacing the integral by a spurious full power of M.",
    }
