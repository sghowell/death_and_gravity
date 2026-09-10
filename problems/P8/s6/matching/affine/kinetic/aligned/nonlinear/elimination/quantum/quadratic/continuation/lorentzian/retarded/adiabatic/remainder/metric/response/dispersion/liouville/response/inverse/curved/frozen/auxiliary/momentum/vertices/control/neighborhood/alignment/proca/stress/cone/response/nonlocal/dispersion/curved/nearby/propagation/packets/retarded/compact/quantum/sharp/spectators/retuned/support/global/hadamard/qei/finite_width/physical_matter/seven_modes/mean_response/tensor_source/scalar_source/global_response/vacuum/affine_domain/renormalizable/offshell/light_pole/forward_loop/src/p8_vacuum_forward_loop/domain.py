"""Actual forward kinematics and conservative full-parameter analytic domain."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model


@cache
def data():
    s, x = sp.symbols("forward_s Feynman_x", real=True)
    a = x * (1 - x)
    ell = 1 - a * s
    # Formal bilinear products, not Hermitian products.
    p11 = p22 = -1
    p12 = (2 - s) / 2
    r1_squared = (1 - x) ** 2 * p11 + x * x * p22 - 2 * x * (1 - x) * p12
    L, H, M = sp.symbols(
        "light_parameter_sum heavy_parameter_sum heavy_mass_squared", positive=True
    )
    a0, a1, a2, a3 = sp.symbols(
        "alpha_light0 alpha_light1 alpha_heavy2 alpha_heavy3", nonnegative=True
    )
    original = (
        a0
        + a1
        + M * (a2 + a3)
        - s * a0 * a1
        - (4 - s) * a2 * a3
        - (a0 + a1) * (a2 + a3)
    )
    grouped = (a0 + a1) ** 2 + M * (a2 + a3) - s * a0 * a1 - (4 - s) * a2 * a3
    alpha3 = 1 - a0 - a1 - a2
    lower = L * L / 4 + M * H - 3 * H * H / 4
    target = sp.Rational(1, 4) + (M - 1) * H + H * (1 - H) / 2
    actual = model.data()["actual_parameters"]["heavy_mass_squared"]
    return {
        "forward_s": s,
        "crossed_u": 4 - s,
        "transfer_t": 0,
        "Feynman_x": x,
        "a": a,
        "light_denominator_shift": ell,
        "complex_disc": "|s-2|<=1; u=4-s, t=0; all light external masses equal one",
        "Hermitian_external_momentum_norm_squared_upper": sp.Rational(3, 2),
        "bilinear_shifted_heavy_external_square": -ell,
        "full_box_parameter_Delta": grouped,
        "full_parameter_real_part_lower": sp.Rational(1, 4) + (M - 1) * H,
        "analytic_continuation": "Defined by the Euclidean amplitude and its Feynman-parameter continuation, not by an unshifted real Euclidean contour after blindly substituting complex external momenta.",
        "checks": {
            "actual_crossed_forward_sum": s + (4 - s) - 4,
            "actual_shifted_external_bilinear_square": sp.expand(r1_squared + ell),
            "full_four_propagator_parameter_polynomial": sp.expand(
                (original - grouped).subs(a3, alpha3)
            ),
            "parameter_pair_product_light_gap": sp.expand(
                (a0 + a1) ** 2 / 4 - a0 * a1 - (a0 - a1) ** 2 / 4
            ),
            "parameter_pair_product_heavy_gap": sp.expand(
                (a2 + a3) ** 2 / 4 - a2 * a3 - (a2 - a3) ** 2 / 4
            ),
            "strict_parameter_lower_bound_decomposition": sp.expand(
                lower.subs(L, 1 - H) - target
            ),
            "light_parameter_lower_quarter_gap": sp.expand(
                1 - 3 * a - sp.Rational(1, 4) - 3 * (x - sp.Rational(1, 2)) ** 2
            ),
        },
        "bounds": {
            "actual_M_above_twenty_four": actual > 24,
            "uniform_heavy_external_shift_ratio_below_one_half": 6 / actual
            < sp.Rational(1, 4),
            "actual_Q_majorant_below_L_majorant": sp.Rational(576, 25) < actual,
        },
    }
