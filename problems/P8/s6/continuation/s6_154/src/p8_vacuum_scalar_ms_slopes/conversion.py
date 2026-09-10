"""Proper interaction-reference change and its differentiated kernels."""

from functools import cache

import sympy as s


@cache
def data():
    x, M, L, g, Q, e, ell, u, v = s.symbols(
        "x M L g Q epsilon ell invariant distance", positive=True
    )
    a = x * (1 - x)
    Delta = x * M + 1 - x - a * u
    D = Delta.subs(u, 1)
    b = a / D
    Bprime = b / Q
    B2prime = x * a / (Q * D**2)
    B2 = x / (Q * (D - a * v))
    B2R = x * b**2 * v**2 / (Q * D * (1 - b * v))
    B = s.Function("mixed_bubble")(u)
    BB = s.Function("squared_heavy_bubble")(u)
    entire = 1 / (Q * e) + ell / Q + e * s.Symbol("next_reference")
    old = L * g * entire * B - g**2 * entire * BB / 2
    new = L * g * B / (Q * e) - g**2 * BB / (2 * Q * e)
    changed = -L * g * ell * B / Q + g**2 * ell * BB / (2 * Q)
    f0, f1 = s.symbols("finite_derivative first_epsilon_derivative")
    # Bubble derivatives are finite; an epsilon reference cannot multiply
    # a derivative pole. The undifferentiated local mass pole is not discarded.
    derivative_limit = s.limit((entire - 1 / (Q * e)) * (f0 + e * f1), e, 0)
    return {
        "on_shell_denominator": D,
        "mixed_B_prime_integrand": Bprime,
        "squared_heavy_B2_prime_integrand": B2prime,
        "MS_minus_old_momentum_dependent_self_energy": changed,
        "exact_squared_heavy_OS_integrand": B2R,
        "MS_slope_conversion_absolute_upper": L * g * ell / (2 * Q**2 * M)
        + g**2 * ell / (4 * Q**2 * M**2),
        "MS_outer_OS_quadratic_conversion_upper": ell
        * (L * g + g**2 / M)
        / (6 * Q**2 * M**2 * (1 - 1 / M)),
        "checks": {
            "finite_reference_change_sign": s.expand(
                s.limit(new - old, e, 0) - changed
            ),
            "mixed_bubble_derivative": s.factor(
                s.diff(-s.log(Delta) / Q, u).subs(u, 1) - Bprime
            ),
            "squared_heavy_derivative": s.factor(-s.diff(Bprime, M) - B2prime),
            "mixed_slope_upper_integral": s.integrate((1 - x) / M, (x, 0, 1))
            - 1 / (2 * M),
            "squared_heavy_slope_upper_integral": s.integrate((1 - x) / M**2, (x, 0, 1))
            - 1 / (2 * M**2),
            "exact_squared_heavy_OS_geometric_remainder": s.factor(
                B2 - B2.subs(v, 0) - v * s.diff(B2, v).subs(v, 0) - B2R
            ),
            "mixed_second_order_geometric_weight": s.integrate(
                (1 - x) ** 2 / (2 * M**2), (x, 0, 1)
            )
            - 1 / (6 * M**2),
            "squared_heavy_second_order_geometric_weight": s.integrate(
                (1 - x) ** 2 / M**3, (x, 0, 1)
            )
            - 1 / (3 * M**3),
            "finite_derivative_has_no_epsilon_times_pole_term": derivative_limit
            - ell * f0 / Q,
            "combined_OS_conversion_coefficient": s.factor(
                L * g * ell / (6 * Q**2 * M**2 * (1 - 1 / M))
                + g**2 * ell / (6 * Q**2 * M**3 * (1 - 1 / M))
                - ell * (L * g + g**2 / M) / (6 * Q**2 * M**2 * (1 - 1 / M))
            ),
        },
        "scope": "Change only proper total canonical interaction references I0 to the pole P. The physical inner OS prescription remains fixed. Local outer mass references have zero derivative, not zero finite mass.",
    }
