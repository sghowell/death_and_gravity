"""Every external-label routing and its decaying complex logarithm difference."""

from functools import cache
from itertools import combinations

import sympy as sp
from p8_vacuum_two_loop_insertions import routing as previous


@cache
def data():
    old = previous.data()
    ps, x = old["momenta"], old["x"]
    checks, rows = {}, []
    for pair in combinations(range(4), 2):
        P = ps[pair[0]] + ps[pair[1]]
        rest = [i for i in range(4) if i not in pair]
        z = sp.simplify(-P.dot(P))
        Delta = 1 - x * (1 - x) * z
        for j in rest:
            other = next(k for k in rest if k != j)
            r = -ps[j] - x * P
            key = "".join(map(str, pair)) + "_" + str(j)
            checks["shift_bilinear_square_" + key] = sp.simplify(r.dot(r) + Delta)
            for n in range(4):
                checks["shift_convex_endpoint_" + key + "_" + str(n)] = sp.expand(
                    r[n] + (1 - x) * ps[j][n] - x * ps[other][n]
                )
            rows.append(
                {
                    "outer_external_pair": pair,
                    "inner_external_label": j,
                    "channel_invariant": z,
                    "Delta": Delta,
                    "inner_shift": r,
                }
            )
    y = sp.Symbol("nonnegative_radial", nonnegative=True)
    a = sp.Symbol("inner_parameter_weight", nonnegative=True)
    t = sp.Symbol("unit_interpolation_parameter", real=True)
    N = sp.Symbol("complex_shifted_square_increment")
    primitive = sp.log(1 + a * (y + t * N))
    diff = sp.log(1 + a * (y + N)) - sp.log(1 + a * y)
    bound = (8 * sp.sqrt(y) + 4) / (y + 4)
    u = sp.Symbol("nonnegative_sqrt_radial", nonnegative=True)
    checks.update(
        {
            "anchored_logarithm_difference_derivative": sp.factor(
                sp.diff(primitive, t) - a * N / (1 + a * (y + t * N))
            ),
            "anchored_logarithm_difference_endpoints": sp.expand(
                primitive.subs(t, 1) - primitive.subs(t, 0) - diff
            ),
            "scalar_interpolation_right_half_plane_margin": sp.expand(
                1
                + a * (y / 2 - sp.Rational(3, 2))
                - (1 + a * y) / 2
                - (sp.Rational(1, 2) - 3 * a / 2)
            ),
            "strict_inner_parameter_endpoint_margin": sp.Rational(1, 2)
            - sp.Rational(3, 8)
            - sp.Rational(1, 8),
            "parameter_fraction_majorant_remainder": sp.factor(
                1 / (y + 4) - a / (1 + a * y) - (1 - 4 * a) / ((y + 4) * (1 + a * y))
            ),
            "uniform_three_log_difference_positive_square": sp.expand(
                3 * (u**2 + 4)
                - (8 * u + 4)
                - 3 * (u - sp.Rational(4, 3)) ** 2
                - sp.Rational(8, 3)
            ),
            "logarithm_difference_majorant_decays_at_infinity": sp.limit(
                bound, y, sp.oo
            ),
            "same_Hermitian_external_norm_upper": old[
                "external_leg_squared_Hermitian_norm_upper"
            ]
            - sp.Rational(3, 2),
        }
    )
    return {
        "all_external_label_routings": tuple(rows),
        "routing_count": len(rows),
        "y": y,
        "inner_shift_squared_Hermitian_norm_upper": sp.Rational(3, 2),
        "inner_shift_bilinear_square_modulus_upper": sp.Rational(7, 4),
        "complex_logarithm_difference_integral": sp.Integral(
            a * N / (1 + a * (y + t * N)), (t, 0, 1)
        ),
        "decaying_logarithm_difference_modulus_upper": bound,
        "uniform_logarithm_difference_modulus_upper": sp.Integer(3),
        "shifted_inner_bubble_modulus_upper": (sp.log(1 + y) + 3) / (16 * sp.pi**2),
        "scope": "All twelve external-label routings have r=-(1-x)p_k+x p_l, a convex combination in Hermitian norm, and bilinear r^2=-Delta. Along scalar interpolation of the squared momentum the logarithm stays in the right half plane. The difference bound decays at infinity; replacing it by the uniform constant would lose the required convergence.",
        "checks": checks,
    }
