"""Full Schwinger normalization, compact derivative and positive layer bounds."""

from functools import cache
from math import factorial

import sympy as sp


@cache
def data():
    U, a, b, c, f0, f1, v, r, q, M, H, L = sp.symbols(
        "positive_U Gaussian_a Gaussian_b Gaussian_c constant_F crossing_F forward_increment positive_scale nonnegative_maximum heavy_mass_squared heavy_parameter_sum light_parameter_sum",
        positive=True,
    )
    k, l = sp.symbols("first_loop_component second_loop_component", real=True)
    gaussian = (sp.pi / sp.sqrt(U) / (2 * sp.pi) ** 2) ** 4
    expo = sp.exp(-(f0 + v * f1) / U)
    checks = {
        "Gaussian_completion_of_square": sp.expand(
            a * k * k
            + 2 * b * k * l
            + c * l * l
            - a * (k + b * l / a) ** 2
            - (c - b * b / a) * l * l
        ),
        "two_loop_four_dimensional_measure_normalization": sp.simplify(
            gaussian * (16 * sp.pi**2) ** 2 * U**2 - 1
        ),
        "second_forward_Taylor_coefficient": sp.simplify(
            sp.diff(expo, v, 2) / 2 - f1 * f1 * expo / (2 * U * U)
        ),
        "scaled_positive_mass_majorant_gap": sp.expand(
            L + (4 * M - 63) * H - (L + 2 * M * H) - (2 * M - 63) * H
        ),
        "positive_exponential_layer_identity": sp.integrate(sp.exp(-r), (r, q, sp.oo))
        - sp.exp(-q),
        "box_parameter_sum_below_seventeen": sp.Rational(17)
        - (4 + sp.Rational(3, 64)) ** 2
        - sp.Rational(2551, 4096),
    }
    by_h = {}
    for h in (1, 2, 3):
        N = 4 + h
        gamma = sp.integrate(r ** (N - 2) * sp.exp(-r), (r, 0, sp.oo))
        pref = sp.Rational(17 * factorial(h + 2) * 16, 41472) * sp.Rational(2, 3) ** h
        by_h[h] = {
            "edge_count": N,
            "fourfold_parameter_scaling": 4 ** (h + 2),
            "positive_layer_factor": gamma,
            "rational_coupling_and_loop_prefactor": pref,
        }
        checks["parameter_scaling_exponent_" + str(h)] = N - 4 + 2 - (h + 2)
        checks["homogeneous_layer_integral_" + str(h)] = gamma - factorial(h + 2)
        checks["heavy_coupling_suppression_" + str(h)] = (
            sp.Rational(4 ** (h + 2), 6**h) - 16 * sp.Rational(2, 3) ** h
        )
    return {
        "Gaussian_measure": gaussian,
        "second_forward_coefficient_integrand": f1 * f1 * expo / (2 * U**4),
        "scaled_mass_lower": L + 2 * M * H,
        "box_moment_upper": sp.Integer(17),
        "per_heavy_edge_count": by_h,
        "scope": "Absolute compact-disc bounds of individually UV-finite bare graphs, using the full heavy kernel. No subtraction-dependent graph or two-loop pole normalization is included.",
        "checks": checks,
    }
