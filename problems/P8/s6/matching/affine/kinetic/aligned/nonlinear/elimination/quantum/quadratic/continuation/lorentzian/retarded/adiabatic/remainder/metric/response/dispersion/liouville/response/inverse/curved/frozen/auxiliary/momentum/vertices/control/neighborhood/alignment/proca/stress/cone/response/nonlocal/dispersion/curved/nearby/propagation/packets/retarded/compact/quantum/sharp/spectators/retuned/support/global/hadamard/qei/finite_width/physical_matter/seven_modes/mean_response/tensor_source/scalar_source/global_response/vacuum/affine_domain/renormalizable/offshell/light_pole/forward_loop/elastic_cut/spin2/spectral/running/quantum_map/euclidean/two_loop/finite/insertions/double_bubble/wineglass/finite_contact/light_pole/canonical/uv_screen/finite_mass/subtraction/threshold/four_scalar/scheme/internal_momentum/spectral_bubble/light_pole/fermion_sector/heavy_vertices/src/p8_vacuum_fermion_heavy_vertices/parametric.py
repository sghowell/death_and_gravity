"""Mass replacement and positive simplex bounds on the full forward disc."""

from functools import cache

import sympy as sp


@cache
def data():
    a0, a1, a2, a3 = sp.symbols("alpha0 alpha1 alpha2 alpha3", nonnegative=True)
    s, v, M = sp.symbols(
        "forward_s spectral_mass_squared heavy_mass_squared", real=True
    )
    H = a2 + a3
    L = a0 + a1
    original = v * a0 + a1 + M * H - s * a0 * a1 - (4 - s) * a2 * a3 - L * H
    grouped = L**2 + M * H - s * a0 * a1 - (4 - s) * a2 * a3 + (v - 1) * a0
    lower = L**2 / 4 + M * H - 3 * H**2 / 4 + (v - 1) * a0
    target = sp.Rational(1, 4) + (M - 1) * H + (v - 1) * a0 + H * (1 - H) / 2
    x, a, k = sp.symbols("spectral_simplex_x base_gap mass_increment", positive=True)
    tri_volume = 1 - x
    box_volume = (1 - x) ** 2 / 2
    log_primitive = sp.log(a + k * x) / k
    inv_primitive = -1 / (k * (a + k * x))
    checks = {
        "unequal_mass_full_box_parameter_polynomial": sp.expand(
            (original - grouped).subs(a3, 1 - a0 - a1 - a2)
        ),
        "quarter_gap_decomposition": sp.factor(
            (lower - target).subs(a3, 1 - a0 - a1 - a2)
        ),
        "light_pair_product_gap": sp.expand(L**2 / 4 - a0 * a1 - (a0 - a1) ** 2 / 4),
        "heavy_pair_product_gap": sp.expand(H**2 / 4 - a2 * a3 - (a2 - a3) ** 2 / 4),
        "triangle_remaining_simplex_volume": sp.integrate(1, (a1, 0, 1 - x))
        - tri_volume,
        "box_remaining_simplex_volume": sp.integrate(1 - x - a1, (a1, 0, 1 - x))
        - box_volume,
        "triangle_volume_below_one": 1 - tri_volume - x,
        "box_volume_below_half": sp.expand(
            sp.Rational(1, 2) - box_volume - x * (2 - x) / 2
        ),
        "triangle_log_primitive": sp.diff(log_primitive, x) - 1 / (a + k * x),
        "box_inverse_primitive": sp.diff(inv_primitive, x) - 1 / (a + k * x) ** 2,
        "triangle_endpoint_log_ratio": sp.expand_log(
            (log_primitive.subs(x, 1) - log_primitive.subs(x, 0))
            - sp.log((a + k) / a) / k,
            force=False,
        ),
        "box_endpoint_integral": sp.factor(
            inv_primitive.subs(x, 1) - inv_primitive.subs(x, 0) - 1 / (a * (a + k))
        ),
        "mass_increment_half_gap": v - 1 - v / 2 - (v / 2 - 1),
        "log_argument_comparison_gap": 4 * v - (1 + 4 * (v - 1)) - 3,
        "Gamma_three_propagator_factor": sp.gamma(1) - 1,
        "Gamma_four_propagator_factor": sp.gamma(2) - 1,
    }
    return {
        "forward_domain": "|s-2|<=1; u=4-s; t=0",
        "full_parameter_polynomial": grouped,
        "strict_real_part_lower": sp.Rational(1, 4) + (M - 1) * H + (v - 1) * a0,
        "triangle_simplex_volume": tri_volume,
        "box_simplex_volume": box_volume,
        "triangle_dimensionless_upper": "2 log(4v)/v",
        "box_dimensionless_upper": "4/v",
        "scope": "Each scalar one-loop master excludes its common 1/Q factor. Repeated heavy routes and all triangle/bubble parameter faces are included. No real-contour translation at complex external momenta is assumed.",
        "checks": checks,
    }
