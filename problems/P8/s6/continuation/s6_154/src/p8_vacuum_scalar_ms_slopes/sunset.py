"""Regulator-first local sunset derivative and convergent simplex moments."""

from functools import cache

import sympy as s


@cache
def data():
    x, y, z, r, t = s.symbols("x y z r t", positive=True)
    e, ell, q, L, u = s.symbols("epsilon ell Q L invariant")
    U = x * y + x * z + y * z
    P = x * y * z / U
    w = x * y * z / U**3
    v, a = 1 + t + r * t, 1 + r * (1 + t)
    substitution = {x: 1 / a, y: r / a, z: r * t / a}
    sector_weight = t / v**3
    sector_U = r * v / a**2
    sector_P = r * t / (a * v)
    jacobian = s.det(s.Matrix([r / a, r * t / a]).jacobian([r, t]))
    T0 = s.Rational(1, 2)
    T1 = 6 * s.Integral(sector_weight * s.log(sector_U), (r, 0, 1), (t, 0, 1))
    C = -6 * s.Integral(sector_weight * s.log(1 - sector_P), (r, 0, 1), (t, 0, 1))
    finite = L**2 / (6 * q**2) * (ell * T0 + T1 / 2 + C)
    radial_primitive = -1 / (2 * v**2)
    radial_integral = 1 / (2 * (1 + t) ** 2) - 1 / (2 * (1 + 2 * t) ** 2)
    dim_factor = s.exp(2 * s.EulerGamma * e + 2 * ell * e) * s.gamma(1 + 2 * e) / 2
    # e Gamma(2e) is regular. Its derivative is ell, not ell/2.
    dim_first = s.diff(dim_factor, e).subs(e, 0)
    D = 1 - u * P
    return {
        "symbols": {"epsilon": e, "ell": ell, "Q": q, "L": L, "invariant": u},
        "simplex_U": U,
        "simplex_P": P,
        "simplex_weight": w,
        "sector_U": sector_U,
        "sector_P": sector_P,
        "sector_weight": sector_weight,
        "slope_T0": T0,
        "slope_T1": T1,
        "on_shell_finite_difference_C": C,
        "full_dimensional_slope_integrand": L**2
        / (6 * q**2)
        * s.exp(2 * s.EulerGamma * e + 2 * ell * e)
        * s.gamma(2 * e)
        * w
        * U**e
        * D ** (-2 * e),
        "UV_slope_pole": L**2 / (24 * q**2 * e),
        "finite_MS_on_shell_slope": finite,
        "checks": {
            "Schwinger_derivative_before_Gamma_recurrence": s.factor(
                s.diff(D ** (1 - 2 * e), u) - (2 * e - 1) * P * D ** (-2 * e)
            ),
            "Gamma_derivative_recurrence": s.expand_func(s.gamma(2 * e))
            - (2 * e - 1) * s.gamma(2 * e - 1),
            "projective_simplex_jacobian": s.factor(jacobian - r / a**3),
            "sector_weight_pullback": s.factor(
                w.subs(substitution) * jacobian - sector_weight
            ),
            "sector_U_pullback": s.factor(U.subs(substitution) - sector_U),
            "sector_P_pullback": s.factor(P.subs(substitution) - sector_P),
            "radial_primitive_derivative": s.factor(
                s.diff(radial_primitive, r) - sector_weight
            ),
            "radial_primitive_endpoints": s.factor(
                radial_primitive.subs(r, 1)
                - radial_primitive.subs(r, 0)
                - radial_integral
            ),
            "exact_simplex_weight_integral": s.integrate(6 * radial_integral, (t, 0, 1))
            - T0,
            "regularized_Gamma_constant": dim_factor.subs(e, 0) - s.Rational(1, 2),
            "regularized_Gamma_first_coefficient": dim_first - ell,
            "full_epsilon_integrand_first_coefficient": s.diff(
                s.exp(e * s.log(U) - 2 * e * s.log(1 - P)), e
            ).subs(e, 0)
            - s.log(U)
            + 2 * s.log(1 - P),
            "slope_pole_normalization": L**2 * T0 / (12 * q**2 * e)
            - L**2 / (24 * q**2 * e),
            "on_shell_Symanzik_factorization": s.expand(
                (x + y + z) * U - x * y * z - (x + y) * (y + z) * (z + x)
            ),
        },
        "scope": "Inverse convention D=s-1+Pi. The sunset is +L^2/6 in Pi. Proper local four-point contractions are momentum independent. The displayed finite MS slope is recorded before outer OS; no exact numerical conjecture for T1/2+C is assumed.",
    }
