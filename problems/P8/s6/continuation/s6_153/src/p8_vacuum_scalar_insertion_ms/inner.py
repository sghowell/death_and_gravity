"""The full dimensional on-shell slope fixes the evanescent coefficient."""

from functools import cache

import sympy as s
from p8_vacuum_light_pole import kernel as previous


@cache
def data():
    x, M, g, Q, e, ell, z = s.symbols("x M g Q epsilon ell z", positive=True)
    a = x * (1 - x)
    Delta = x * M + 1 - x - a * z
    D = Delta.subs(z, 1)
    b = a / D
    raw = s.exp(s.EulerGamma * e + ell * e) * s.gamma(e) * Delta ** (-e)
    slope = s.exp(s.EulerGamma * e + ell * e) * s.gamma(1 + e) * b * D ** (-e)
    old = previous.data()
    mapped = old["normalized_parameter_weight"].subs(
        {old["Feynman_parameter"]: x, old["heavy_mass_squared"]: M}
    )
    A0 = s.Integral(b, (x, 0, 1))
    A1 = s.Integral(b * (ell - s.log(D)), (x, 0, 1))
    return {
        "symbols": {"x": x, "M": M, "g": g, "Q": Q, "epsilon": e, "ell": ell},
        "mixed_bubble_Delta": Delta,
        "on_shell_D": D,
        "normalized_b": b,
        "full_dimensional_mixed_bubble_integrand": g * raw / Q,
        "full_dimensional_slope_integrand": g * slope / Q,
        "alpha_zero": g * A0 / Q,
        "alpha_first_epsilon": g * A1 / Q,
        "checks": {
            "same_actual_parent_normalized_weight": s.factor(mapped - b),
            "on_shell_denominator": s.expand(D - x * M - (1 - x) ** 2),
            "dimensional_derivative_before_Gamma_recurrence": s.factor(
                s.diff(Delta ** (-e), z).subs(z, 1) - e * b * D ** (-e)
            ),
            "Gamma_slope_recurrence": s.expand_func(s.gamma(1 + e)) - e * s.gamma(e),
            "dimensional_slope_first_coefficient": s.diff(
                s.exp(s.EulerGamma * e + ell * e) * s.gamma(1 + e) * D ** (-e), e
            ).subs(e, 0)
            - ell
            + s.log(D),
            "normalized_weight_positive_upper_gap": s.factor(
                (1 - x) / M - b - (1 - x) ** 3 / (M * D)
            ),
            "denominator_above_x_M": s.expand(D - x * M - (1 - x) ** 2),
            "denominator_above_one": s.expand(D - 1 - x * (M - 2 + x)),
            "denominator_below_M": s.expand(M - D - (1 - x) * (M - 1 + x)),
            "integrated_normalized_weight_upper": s.integrate((1 - x) / M, (x, 0, 1))
            - 1 / (2 * M),
        },
        "scope": "The fixed physical inner on-shell subtraction is unchanged. Its asymptotic multiplier is the full regulated slope alpha_D, not just alpha_0. The Gamma recurrence and first epsilon coefficient are taken before the outer pole multiplication.",
    }
