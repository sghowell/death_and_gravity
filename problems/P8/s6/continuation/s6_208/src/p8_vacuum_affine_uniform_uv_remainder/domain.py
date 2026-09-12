"""Explicit joint inverse-radius/time domain for the actual normalized W8 modes."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import reference
from p8_vector_state.comparison import AMAX

EPS = s.Rational(1, 100)
MASS = reference.MASS
OUTER = reference.OUTER
RADIUS = reference.RADIUS
FIELD = s.Integer(20)
PAIR = s.Integer(1000)


@cache
def constants():
    ea = 6 * OUTER / (1 - 6 * OUTER)
    eg = 2 * EPS + EPS**2
    ew = ea + (1 + ea) * eg + AMAX**2 * EPS**2
    eta = (EPS + eg / (1 + s.sqrt(1 - eg))) / s.sqrt(1 - eg)
    ratios = {
        kind: sum(
            row[j]["complex_coefficient_upper"]
            * (AMAX * EPS / (s.Rational(98, 100) * MASS)) ** (2 * j)
            for j in range(1, 5)
        )
        for kind, row in reference.data()[
            "actual_W8_complex_coefficient_majorants"
        ].items()
    }
    return {
        "inverse_radius_numerator": EPS,
        "relative_time_scale_defect": ea,
        "scaled_second_leg_square_defect": eg,
        "complete_scaled_squared_frequency_defect": ew,
        "scaled_z_modulus_upper": 1 + AMAX**2 * EPS**2 / (1 - ew),
        "scaled_omega_real_lower": s.Rational(98, 100) / AMAX,
        "full_scaled_W_relative_defects": ratios,
        "scaled_W_real_lower": s.Rational(98, 100) / AMAX - s.Rational(2, 1000),
        "scaled_W_absolute_upper": s.Integer(3),
        "scaled_W_time_derivative_upper": 3 / RADIUS,
        "scaled_canonical_momentum_upper": 3 + (EPS / MASS) * (3 / RADIUS + 5),
        "complex_direction_displacement_upper": eta,
        "complex_frame_rotation_displacement_upper": 2 * eta + 4 * eta**2 / (2 - eta),
        "transverse_magnetic_norm_upper": 2
        * s.sqrt(2)
        * (1 + EPS)
        * s.Rational(101, 100),
        "longitudinal_electric_norm_upper": EPS * 4 * 2 / (s.Rational(98, 100) / AMAX),
        "longitudinal_temporal_multiplier_upper": (1 + EPS)
        * s.Rational(101, 100)
        / (s.Rational(98, 100) / AMAX),
        "transverse_field_norm_squared_upper": 8**2 + 4**2 + (2 * EPS) ** 2,
        "longitudinal_field_norm_squared_upper": 1 + 8**2 + 4**2,
    }


@cache
def data():
    c = constants()
    x = s.symbols("x")
    n = s.Matrix(s.symbols("n0:3"))
    P = s.Matrix(s.symbols("P0:3"))
    v = -n + x * P
    checks = {
        "complete_second_leg_scaled_square": s.expand(
            (v.T * v)[0] - (n.T * n)[0] + 2 * x * (n.T * P)[0] - x * x * (P.T * P)[0]
        ),
        "same_outer_inner_time_discs": OUTER - 2 * RADIUS,
        "same_actual_mass": MASS - 1000,
        "full_sixteen_stress_pair_norm_count": 2 * FIELD**2 - 800,
        "maximum_scaled_mass_radius": MASS * (EPS / MASS) - EPS,
    }
    gates = {
        "full_joint_frequency_defect_below_one_fortieth": c[
            "complete_scaled_squared_frequency_defect"
        ]
        < s.Rational(1, 40),
        "positive_scaled_frequency_branch_lower": s.Rational(98, 100) ** 2
        < 1 - c["complete_scaled_squared_frequency_defect"],
        "scaled_frequency_absolute_upper_two": 1
        + c["complete_scaled_squared_frequency_defect"]
        < 4,
        "unchanged_polynomial_z_domain": c["scaled_z_modulus_upper"]
        < s.Rational(101, 100),
        "scaled_W_defects_below_one_thousandth": all(
            v < s.Rational(1, 1000)
            for v in c["full_scaled_W_relative_defects"].values()
        ),
        "old_coefficient_majorants_dominate_new_ratio": AMAX * EPS * s.Rational(99, 100)
        < s.Rational(98, 100),
        "both_scaled_W_real_parts_above_half": c["scaled_W_real_lower"]
        > s.Rational(1, 2),
        "canonical_momentum_below_four": c["scaled_canonical_momentum_upper"] < 4,
        "complex_direction_within_one_fortieth": c[
            "complex_direction_displacement_upper"
        ]
        < s.Rational(1, 40),
        "complex_frame_vectors_below_two": c[
            "complex_frame_rotation_displacement_upper"
        ]
        < 1,
        "transverse_magnetic_norm_below_four": c["transverse_magnetic_norm_upper"] < 4,
        "longitudinal_electric_norm_below_one": c["longitudinal_electric_norm_upper"]
        < 1,
        "longitudinal_temporal_multiplier_below_two": c[
            "longitudinal_temporal_multiplier_upper"
        ]
        < 2,
        "full_transverse_field_norm_below_twenty": c[
            "transverse_field_norm_squared_upper"
        ]
        < FIELD**2,
        "full_longitudinal_field_norm_below_twenty": c[
            "longitudinal_field_norm_squared_upper"
        ]
        < FIELD**2,
        "complete_pair_norm_below_one_thousand": 2 * FIELD**2 < PAIR,
    }
    return {
        "domain": "For real unit n and arbitrary real external P, set U=m+|P| and rho=epsilon/U, epsilon=1/100. Combine |x|<=rho with the original outer/inner complex time discs. The scaled momenta are n,-n+xP. This is a new joint inverse-radius domain, not an extension of the earlier real-internal-k domain by assumption.",
        "full_frequency": "The exact scaled squared-frequency defect is bounded by ea+(1+ea)(2epsilon+epsilon^2)+Amax^2 epsilon^2<1/40. The positive analytic branches satisfy Re Omega>.98/Amax and |Omega|<2; |z|<1.01. All four actual W8 polynomial bounds remain valid and give relative What defect<.001, Re What>.5 and |What|<3. Hence each normalized inverse summed phase has modulus<1.",
        "full_readouts": "Cauchy on the inner time disc gives |What_time|<=60000, Fhat<=1 and |Pihat|<4 for either Schwarz sign. A local analytic complex frame has genuine Euclidean vector norms<2. Transverse norms obey E<8,B<4,mxA<.02; longitudinal norms obey E<1,|mA0|<8,mAsp<4. Every full ten-field norm is<20 and the full sixteen-component pair norm is<1000.",
        "scope": "These improved constants apply to NORMALIZED analytic modes in the stated inverse-radius domain. They are not substituted for the older all-real unscaled mode bound, and do not change the full W8 reference or actual state.",
        "constants": c,
        "checks": checks,
        "gates": {k: bool(v) for k, v in gates.items()},
    }
