"""A regulator-first sector decomposition of the equal-mass local reference."""

from functools import cache

import sympy as s


@cache
def data():
    r, t, e, u = s.symbols("r t epsilon u", positive=True)
    v, w = 1 + t + r * t, 1 + r * (1 + t)
    B = w ** (1 - 2 * e) / v ** (2 - e)
    b = (1 + t) ** (-2 + e)
    B0, b0 = B.subs(e, 0), b.subs(e, 0)
    D = s.factor((B0 - b0) / r)
    logder = s.log(v) - 2 * s.log(w)
    R1 = s.log(r) * D + (B0 * logder - b0 * s.log(1 + t)) / r
    P = (1 - 2 ** (e - 1)) / (1 - e)
    p2 = (1 - s.log(2) - s.log(2) ** 2 / 2) / 2
    first_primitive = (
        -(t * t + t + 1) / (t * (t + 1) * (r * t + t + 1))
        - s.log(r * t + t + 1) / (t + 1) ** 2
    )
    first_integrated = (t * t + t + 1) / ((t + 1) ** 2 * (1 + 2 * t)) + s.log(
        (1 + t) / (1 + 2 * t)
    ) / (1 + t) ** 2
    transformed = (
        -1 + s.Rational(1, 2) / (1 - u) + s.Rational(3, 2) / (1 + u) - s.log(1 + u)
    )
    transformed_primitive = 1 - s.log(1 - u) / 2 + (s.Rational(1, 2) - u) * s.log(1 + u)
    j = s.Symbol("finite_sector_constant_j", real=True)
    return {
        "symbols": {"r": r, "t": t, "epsilon": e, "u": u, "j": j},
        "sector_B": B,
        "endpoint_B": b,
        "endpoint_integral_P": P,
        "subtracted_R_integrand": r ** (e - 1) * (B - b),
        "R_first_derivative_integrand": R1,
        "finite_constant": p2 + s.Integral(R1, (r, 0, 1), (t, 0, 1)),
        "finite_constant_absolute_upper": s.Integer(22),
        "A_through_first_order": 3 / e + 3 + 6 * j * e,
        "independent_radial_bubble_hypergeometric_A": 3
        / (e * (1 + e))
        * s.Integral(s.hyper((e, 2 - e), (2 + e,), 1 - t * (1 - t)), (t, 0, 1)),
        "checks": {
            "sector_jacobian_r": s.det(s.Matrix([[1, 0], [t, r]])) - r,
            "three_symmetric_primary_two_ordered_secondary_sectors": 3 * 2 - 6,
            "B_first_epsilon_derivative": s.diff(B, e).subs(e, 0) - B0 * logder,
            "endpoint_first_epsilon_derivative": s.diff(b, e).subs(e, 0)
            - b0 * s.log(1 + t),
            "positive_D_numerator": s.factor(
                D - (1 + t + t * t + t**3 - r * t * t) / ((1 + t) ** 2 * v**2)
            ),
            "endpoint_P_zero": P.subs(e, 0) - s.Rational(1, 2),
            "endpoint_P_first": s.diff(P, e).subs(e, 0) - (1 - s.log(2)) / 2,
            "endpoint_P_second_coefficient": s.diff(P, e, 2).subs(e, 0) / 2 - p2,
            "radial_R_zero_primitive": s.factor(s.diff(first_primitive, r) - D),
            "radial_R_zero_anchors": s.expand_log(
                s.simplify(
                    first_primitive.subs(r, 1)
                    - s.limit(first_primitive, r, 0)
                    - first_integrated
                ),
                force=True,
            ).simplify(),
            "rationalized_t_integral": s.expand_log(
                s.simplify(
                    first_integrated.subs(t, u / (1 - u)) / (1 - u) ** 2 - transformed
                ),
                force=True,
            ).simplify(),
            "transformed_R_zero_primitive": s.factor(
                s.diff(transformed_primitive, u) - transformed
            ),
            "R_zero_exact_endpoint": s.simplify(
                transformed_primitive.subs(u, s.Rational(1, 2))
                - transformed_primitive.subs(u, 0)
                - s.log(2) / 2
            ),
            "A_constant_after_endpoint_and_R_zero": s.simplify(
                6 * ((1 - s.log(2)) / 2 + s.log(2) / 2) - 3
            ),
            "R_first_absolute_majorant": s.Rational(4)
            + 4
            + 3 * s.Rational(9, 2)
            - s.Rational(43, 2),
            "finite_j_absolute_majorant": s.Rational(43, 2) + s.Rational(1, 2) - 22,
        },
        "scope": "Endpoint subtraction gives A(e)=6[P(e)/e+R(e)] with R holomorphic for Re e>-1 near zero. The finite constant is a convergent two-parameter integral bounded analytically, not a fitted decimal. Equal-mass vacuum differentiation is justified first on 1/2<Re e<1 and then meromorphically continued.",
    }
