"""Uniform complex full-clock source bounds and its mixed spatial variation."""

from functools import cache

import sympy as s
from p8_vacuum_affine_retarded_energy import clock

DELTA = s.Rational(1, 100)
RADIUS = s.Rational(1, 100)
GRAD = s.Integer(400000000)


@cache
def data():
    d = s.Rational(1, 50)
    u = s.Rational(13, 25)
    x = s.Rational(1, 20)
    hd = 1 - u * u
    Rminus = 3 * x / hd**3
    RX = (60 * x + 3) / hd**3
    Ru = 6 * u * s.Rational(2, 5) / hd
    CZ = 1 / (1 - x) ** 2 + 3 * 16 / (2 * s.Rational(1, 2) * (1 - x))
    box = d + 3 * s.Rational(8, 5) * (1 + d) + 3 * d
    Z = (
        (1 + d) ** 2 * d
        + 6 * (1 + d) * d * d
        + 9 * d**3
        + s.Rational(8, 5) * (1 + d) * 3 * d * d
    )
    E = 5 + 3 * 3 * (1 + x)
    Q = 15 / (1 - x) + 3 + 52 * s.Rational(3, 100)
    S = (1 + d) * s.Rational(2, 5) * 21
    second = 2 * 10 / RADIUS**2
    first = 15 * second
    each = 15 * first + 15**2 * second
    gradient = 4 * each
    checks = {
        "complete_source_zero_first_variation": clock.data()["checks"][
            "source_first_clock_variation_zero"
        ],
        "full_complex_gradient_square_defect": 2 * d + 4 * d * d - s.Rational(26, 625),
        "second_jet_Cauchy_coefficient": second - 200000,
        "first_jet_zero_background_segment_bound": first - 3000000,
        "one_component_mixed_spatial_derivative_count": each - 90000000,
        "full_gradient_output_count": gradient - 360000000,
        "unchanged_sharp_source_Frechet_bound": clock.bounds()[
            "uniform_spatial_and_Frechet_L2_constant"
        ]
        - 2048,
    }
    return {
        "class": "The unchanged S6.178 compact CD psi class, all coordinate jets through3<=delta<=1/100, with arbitrary compact smooth eta and its coordinate H3 norm",
        "outer_X_disc_radius": s.Rational(1, 10),
        "inner_X_disc_radius": x,
        "complex_scalar_jet_component_radius": d,
        "Cauchy_polydisc_radius_about_real_small_jet": RADIUS,
        "continuous_complex_bounds": {
            "Rminus": Rminus,
            "RX": RX,
            "Ru": Ru,
            "CZ": CZ,
            "Boxu": box,
            "Z": Z,
            "E": E,
            "Q": Q,
            "source_component": S,
        },
        "source_component_display_upper": s.Integer(10),
        "second_jet_derivative_component_upper": second,
        "first_jet_derivative_component_upper_over_delta": first,
        "full_mixed_spatial_source_derivative_upper_over_delta": GRAD,
        "undifferentiated_source_Frechet_upper_over_delta": s.Integer(2048),
        "checks": checks,
        "gates": {
            "full_X_defect_within_inner_disc": 2 * d + 4 * d * d < x,
            "outer_switch_complement_below_one": 2 * s.Rational(1, 9) ** 1024 < 1,
            "outer_w_exponential_bound": 1024 * s.Rational(11, 10) ** 2 < 2**11
            and s.Rational(4, 5) * 1024 > 512,
            "Rminus_below_two_fifths": Rminus < s.Rational(2, 5),
            "RX_below_sixteen": RX < 16,
            "Ru_below_two": Ru < 2,
            "Href_below_three": 4 * u / hd < 3,
            "CZ_below_fifty_two": CZ < 52,
            "Box_below_five": box < 5,
            "Z_below_three_hundredths": Z < s.Rational(3, 100),
            "E_below_fifteen": E < 15,
            "Q_below_twenty_one": Q < 21,
            "full_source_component_below_ten": S < 10,
            "full_gradient_below_four_e_eight": gradient < GRAD,
            "real_small_ball_plus_Cauchy_disc_fits": DELTA + RADIUS == d,
        },
    }
