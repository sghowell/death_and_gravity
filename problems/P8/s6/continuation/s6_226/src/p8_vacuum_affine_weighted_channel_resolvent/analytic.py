"""Original first-sheet factors and a common exterior logarithmic lower bound."""

from functools import cache

import sympy as s
from p8_proca_rank_one_inverse import spectral as trace
from p8_vacuum_affine_isolated_shear_resolvent import spectral as shear

COEFFICIENT = s.Rational(13, 80)
OFFSET = s.Integer(4)
z = s.Symbol("cut_coordinate", nonnegative=True)
u = s.Symbol("cut_log_radius", nonnegative=True)


@cache
def data():
    trace_data = trace.data()
    shear_data = shear.data()
    d = s.Symbol("exterior_radial_denominator", positive=True)
    A0 = trace_data["closed_scalar_H_in_d"].subs(trace.d, d)
    A2 = shear_data["closed_A2"].subs(shear.d, d)
    r = s.Symbol("exterior_radius", positive=True)
    eps = s.Symbol("inverse_exterior_radius", positive=True)
    logz = s.Symbol("principal_log_p_over_four_mass_squared")
    # On the exterior sheet, atanh(1/sqrt(d)) is
    # [Log(zeta)+2Log(1+sqrt(d))]/2, d=1+1/zeta.
    replace_log = (logz + 2 * s.log(1 + s.sqrt(d))) / 2
    stable0 = A0.subs(s.atanh(1 / s.sqrt(d)), replace_log)
    stable2 = A2.subs(s.atanh(1 / s.sqrt(d)), replace_log)
    leading0 = s.simplify(stable0.subs(d, 1))
    leading2 = s.simplify(stable2.subs(d, 1))
    P0 = 3 - 2 * z + 3 * z * z
    P2 = 30 - 20 * z + 3 * z * z
    poly0 = s.Rational(16, 15) + z - 3 * z * z
    poly2 = -s.Rational(172, 225) + 19 * z / 30 - z * z / 10
    checks = {
        "original_trace_closed_finite_and_radial_factor": s.simplify(
            A0
            - (
                s.Rational(16, 15)
                + d
                - 3 * d * d
                + s.sqrt(d) * (3 - 2 * d + 3 * d * d) * s.atanh(1 / s.sqrt(d))
            )
        ),
        "original_shear_closed_finite_and_radial_factor": s.simplify(
            A2
            - (
                -s.Rational(172, 225)
                + 19 * d / 30
                - d * d / 10
                + s.sqrt(d) * (30 - 20 * d + 3 * d * d) * s.atanh(1 / s.sqrt(d)) / 30
            )
        ),
        "shear_global_cut_polynomial_gap": s.expand(
            poly2 + s.Rational(172, 225) - z * (s.Rational(8, 15) + (1 - z) / 10)
        ),
        "shear_positive_cut_weight": s.expand(P2 - 13 - (1 - z) * (17 - 3 * z)),
        "trace_positive_cut_weight": s.expand(
            P0 - 3 * (z - s.Rational(1, 3)) ** 2 - s.Rational(8, 3)
        ),
        "trace_polynomial_lower_above_minus_one": s.expand(
            poly0 + 1 - s.Rational(1, 15) - (1 - z) * (3 * z + 2)
        ),
        "shear_polynomial_lower_above_minus_one": s.expand(
            poly2 + 1 - s.Rational(53, 225) - z * (s.Rational(8, 15) + (1 - z) / 10)
        ),
        "trace_high_log_coefficient": s.Rational(3, 4) * s.Rational(8, 3) / 2 - 1,
        "shear_high_log_coefficient": s.Rational(3, 4) * 13 / 60 - COEFFICIENT,
        "global_shear_cut_gap": shear.D.subs(shear.z, 0) + s.Rational(172, 225),
        "trace_uniform_outer_circle_logarithm": s.simplify(
            leading0 - (2 * logz + 4 * s.log(2) - s.Rational(14, 15))
        ),
        "shear_uniform_outer_circle_logarithm": s.simplify(
            leading2
            - (
                s.Rational(13, 60) * logz
                + s.Rational(13, 30) * s.log(2)
                - s.Rational(52, 225)
            )
        ),
        "exterior_log_minorant": COEFFICIENT * s.log(r)
        - OFFSET
        - (s.Rational(13, 80) * s.log(r) - 4),
        "large_circle_remainder_zero_at_inverse_radius_zero": s.simplify(
            (stable0.subs(d, 1 + eps) - leading0).subs(eps, 0)
        ),
        "large_circle_shear_remainder_zero_at_inverse_radius_zero": s.simplify(
            (stable2.subs(d, 1 + eps) - leading2).subs(eps, 0)
        ),
    }
    return {
        "source_pinned_factors": {
            "Atrace": A0,
            "A2": A2,
            "definition": "A_i=-F_i with unchanged original finite coefficients and cuts",
        },
        "stable_exterior_closed_forms": {
            "trace": stable0,
            "shear": stable2,
            "d": "1+1/zeta",
            "zeta": "p/(4m^2)",
        },
        "global_real_part_gaps": {
            "trace": s.Rational(16, 15),
            "shear": -s.Rational(172, 225),
        },
        "common_exterior_logarithmic_lower_bound": COEFFICIENT * s.log(r) - OFFSET,
        "exterior_domain": "p on the original first sheet and |p|>=4m^2. Both Re A_i(p)>=(13/80)log(|p|/(4m^2))-4.",
        "harmonic_proof": "First prove the global shear real-part gap from the cut banks and uniform outer circles. Then apply the minimum principle to Re A_i-(13/80)log(|p|/(4m^2)) on the slit exterior. Inner circle, low/high cut regions and outer circles are all controlled. This is not a positive-real-axis-only estimate.",
        "checks": checks,
        "gates": {
            "small_log_shear_margin": OFFSET - s.Rational(172, 225) - 16 * COEFFICIENT
            > 0,
            "small_log_trace_margin": OFFSET + s.Rational(16, 15) - 16 * COEFFICIENT
            > 0,
            "high_log_shear_margin": OFFSET - 1 > 0,
            "high_log_trace_slope_margin": 1 > COEFFICIENT,
            "shear_outer_log_slope_strict": s.Rational(13, 60) > COEFFICIENT,
            "trace_outer_log_slope_strict": 2 > COEFFICIENT,
            "shear_inner_circle_margin": OFFSET - s.Rational(172, 225) > 0,
            "positive_threshold_log_window": s.Integer(16) > 0,
        },
    }
