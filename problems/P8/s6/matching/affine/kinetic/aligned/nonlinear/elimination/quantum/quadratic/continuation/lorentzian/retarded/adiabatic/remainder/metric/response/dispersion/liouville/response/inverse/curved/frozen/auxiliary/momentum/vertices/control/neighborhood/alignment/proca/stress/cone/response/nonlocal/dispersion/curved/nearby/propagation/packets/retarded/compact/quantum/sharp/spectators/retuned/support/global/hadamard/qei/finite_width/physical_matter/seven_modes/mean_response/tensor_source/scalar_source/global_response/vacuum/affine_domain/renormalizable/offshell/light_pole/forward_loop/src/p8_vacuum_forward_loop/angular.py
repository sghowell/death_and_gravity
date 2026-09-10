"""Uniform external-shift expansion at every radial loop momentum."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model


@cache
def data():
    A, z, y, M, c = sp.symbols("A external_shift radial_y M small_c", positive=True)
    g = c * M * M
    inverse = 1 / A - z / A**2 + z * z / (A * A * (A + z))
    V0 = 20 * c * M * (y + 1) / (M + y)
    L = 20 * g * sp.sqrt(y) / (M + y) ** 2
    Q = 192 * g * y / (M + y) ** 3
    first = 7680 * c * c * M**3 * y * (y + 1) / (M + y) ** 4
    second = 1600 * c * c * M**4 * y / (M + y) ** 4
    radial1 = sp.integrate(y / (M + y) ** 4, (y, 0, sp.oo))
    radial2 = sp.integrate(1 / (M + y) ** 4, (y, 0, sp.oo))
    d = model.data()
    p = d["actual_parameters"]
    Mp = p["heavy_mass_squared"]
    gp = p["cubic_coupling_squared"]
    cp = gp / Mp**2
    return {
        "small_parameter": c,
        "cubic_squared": g,
        "radial_vertex_majorant": V0,
        "odd_linear_vertex_majorant": L,
        "quadratic_vertex_remainder_majorant": Q,
        "full_inverse_external_shift_identity": inverse,
        "angular_product_remainder_majorant": first + second,
        "all_three_channels_uniform_amplitude_and_Cauchy_b2_error_upper": 60
        * c
        * c
        * M,
        "actual_angular_b2_error_upper": 60 * cp * cp * Mp,
        "scope": "Only external momentum shifts are expanded. The entire radial y dependence, mixed propagators and heavy threshold remain; no y/M expansion or artificial loop cutoff is made.",
        "checks": {
            "literal_untruncated_external_shift_inverse": sp.factor(
                1 / (A + z) - inverse
            ),
            "radial_cross_term_majorant_coefficient": sp.factor(2 * V0 * Q - first),
            "four_linear_square_majorant_coefficient": sp.factor(4 * L * L - second),
            "first_positive_radial_gap": sp.factor(
                4 * y
                - y * y * (y + 1) / (y + sp.Rational(1, 4)) ** 2
                - y * (3 * y * y + y + sp.Rational(1, 4)) / (y + sp.Rational(1, 4)) ** 2
            ),
            "second_positive_radial_gap": sp.factor(
                1
                - y * y / (y + sp.Rational(1, 4)) ** 2
                - (y / 2 + sp.Rational(1, 16)) / (y + sp.Rational(1, 4)) ** 2
            ),
            "first_exact_radial_integral": radial1 - 1 / (6 * M * M),
            "second_exact_radial_integral": radial2 - 1 / (3 * M**3),
            "all_three_channel_prefactor_bound": sp.Rational(3, 288)
            * (5120 + sp.Rational(1600, 3))
            - sp.Rational(530, 9),
            "linear_shift_constant_square_gap": sp.Integer(20) ** 2
            - 16**2 * sp.Rational(3, 2)
            - 16,
            "Q_to_L_ratio_squared_maximum": sp.Rational(192, 40) ** 2 / M
            - sp.Rational(576, 25) / M,
        },
        "bounds": {
            "twenty_bounds_sixteen_times_square_root_three_halves": 20**2
            > 16**2 * sp.Rational(3, 2),
            "all_channel_Cauchy_constant_below_sixty": sp.Rational(530, 9) < 60,
            "actual_angular_error_strictly_positive": 60 * cp * cp * Mp > 0,
        },
    }
