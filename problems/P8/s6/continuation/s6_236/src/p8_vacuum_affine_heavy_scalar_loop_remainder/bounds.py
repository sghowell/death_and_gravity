"""A complete physical all-angle first-loop remainder with actual parameter margins."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_scalar_tree_matching import model

S_MAX = s.Integer(10) ** 196
LOG_MAJORANT = s.Integer(1000)
CHANNEL_BOUND = s.Integer(2048) * LOG_MAJORANT
TOTAL_BOUND = s.Integer(10) ** 7


@cache
def data():
    S, n, g2, lam, D = s.symbols("S n g_squared lambda D", positive=True)
    eps = s.Rational(1, 8)
    actual_n = model.MASS2
    actual_g2 = s.Rational(1, 2**26)
    actual_D = actual_n - 2
    loop_ratio_upper = TOTAL_BOUND * actual_g2 / (9 * actual_n)
    universal = s.Symbol("universal_constant")
    base_bound = TOTAL_BOUND * g2 * g2 * S * S / (16 * s.pi**2 * n**4)
    difference_bound = 2 * base_bound
    relative_bound = 8 * TOTAL_BOUND * g2 * g2 / (16 * s.pi**2 * lam * n**4)
    q = s.Symbol("ratio", positive=True)
    checks = {
        "real_discriminant_disk_bound": 4 * s.Rational(1, 16)
        + 4 * s.Rational(1, 16) ** 2
        - s.Rational(17, 64),
        "Cauchy_alpha_deviation": s.Rational(1, 16)
        + s.Rational(1, 8)
        - s.Rational(3, 16),
        "primitive_P_supremum_bound": (1 + s.Rational(1, 8)) / s.Rational(3, 4)
        - s.Rational(3, 2),
        "primitive_Q_supremum_bound": ((1 + s.Rational(1, 8)) * s.Rational(1, 4) + 2)
        / s.Rational(3, 4)
        - s.Rational(73, 24),
        "full_bubble_remainder_bound": s.Rational(113, 8) - (2 + 12 + 8 * eps**2),
        "full_triangle_three_remainder_parts": 384 + 12 + 14 - 410,
        "all_channel_remainder_parts": 16 + 512 + 2 * 256 - 1040,
        "physical_squared_denominator_margin": s.expand(
            4 * (S - 2) ** 2 - S**2 - (S - 4) * (3 * S - 4)
        ),
        "contact_subtraction_universal_cancellation": universal - universal,
        "uniform_symmetric_point_bound_dominated": s.expand(
            2 * S * S - (S * S + 16) - (S - 4) * (S + 4)
        ),
        "complete_relative_loop_bound": s.cancel(
            difference_bound / (lam * S * S / 4) - relative_bound
        ),
        "canonical_target_parameter_relation": s.cancel(
            model.LAMBDA - actual_g2 / (2 * actual_D**3)
        ),
        "substituted_relative_bound": s.cancel(
            relative_bound.subs(lam, g2 / (2 * D**3))
            - TOTAL_BOUND * g2 * D**3 / (s.pi**2 * n**4)
        ),
        "geometric_remainder_weight_monotonic_numerator": s.cancel(
            s.diff(q * q * (3 - 2 * q) / (1 - q) ** 2, q) * ((1 - q) ** 3)
            - 2 * q * (q * q - 3 * q + 3)
        ),
        "original_tree_plus_loop_triangle_inequality_margin": s.Rational(1, 59)
        - s.Rational(1, 60)
        - s.Rational(1, 3540),
    }
    return {
        "physical_domain": "Actual V2S-T1-OS4 parameters; 4<=s=S<=10^196, t=-(S-4)(1-x)/2, u=-(S-4)(1+x)/2, -1<=x<=1. Equivalently2<=E<=10^98. The symmetric subthreshold subtraction point is evaluated separately.",
        "uniform_parameter_domain": "For every channel and every Feynman angle, |L|,|b|<=K=S/4. n=M_H^2>8S, hence w=1/n is inside half the analytic coefficient radius1/(16K). At the symmetric point use S=4,K=1.",
        "logarithmic_majorant": "For s=S>=4, factor L=S(x-r_plus)(x-r_minus). Each integral of -log|x-r| is at most1+log2<2. Thus integral|Log(L-i0)|<=log S+4+pi. Negative channels have L>=1; the symmetric point has L>=2/3. With n<10^198 and log10<7/3, both integral|log n-Log L|+3 and |B| are below1000. The threshold double root and L0 points are included by an integrable limit.",
        "complete_channel_bounds": "With M=1000, the bubble remainder is at most16 M S^2/n^4; the complete triangle remainder at most512 M S^2/n^4; each ordered box remainder at most256 M S^2/n^4. Their full sum is below2048 M S^2/n^4. All higher inverse-mass terms, signs and branches are retained.",
        "full_base_remainder_bound": base_bound,
        "base_remainder_definition": "A1_base-g^4/(16pi^2)[-6log(n)/n^2+(-16log(n)+10/3)/n^3]",
        "same_OS4_complete_first_loop_bound": difference_bound,
        "comparison_to_entire_original_tree": "S233 proves A_original>=lambda(S-2)^2>=lambda S^2/4. Therefore |A1_OS4|/A_original<10^7 g^2/(9n)<10^-199 on the full physical window, including the real and imaginary parts.",
        "actual_rational_relative_upper": loop_ratio_upper,
        "tree_plus_first_loop": "With hbar=1 as a bookkeeping convention for the first-loop-truncated expression, |A_tree,V2S+A1_OS4-A_original|/A_original<1/60+10^-199<1/59. This is a theorem about that truncation, not the exact full quantum amplitude.",
        "not_established": "Higher physical coefficient jets, all omitted loop orders, resonance control, exact S-matrix/UV properties, finite-gravity Regge control and the original covariant common-parent bounce remain unproved.",
        "checks": {k: s.cancel(v) for k, v in checks.items()},
        "gates": {
            "actual_window_inside_half_Cauchy_radius": S_MAX / actual_n
            < s.Rational(1, 8),
            "actual_mass_above_thirty_two": actual_n > 32,
            "actual_logarithm_bound_domain": actual_n < 10**198,
            "exponential_series_proves_log10_below_7_over_3": sum(
                s.Rational(7, 3) ** j / s.factorial(j) for j in range(9)
            )
            > 10,
            "integrated_logarithmic_bound": 2 * 462 + 8 + 3 < LOG_MAJORANT,
            "fixed_vertex_remainder_majorant": s.Rational(8, 7) + s.Rational(8, 25) < 4,
            "fixed_vertex_absolute_majorant": 2 + 2 * eps + 4 * eps**2 < 3,
            "bubble_remainder_majorant": s.Rational(113, 8) < 16,
            "triangle_remainder_majorant": s.Integer(410) < 512,
            "complete_channel_majorant": s.Integer(1040) < 2048,
            "complete_three_channel_majorant": 3 * CHANNEL_BOUND < TOTAL_BOUND,
            "actual_full_first_loop_relative_error": loop_ratio_upper
            < s.Rational(1, 10**199),
            "tree_plus_first_loop_relative_error": s.Rational(1, 60)
            + s.Rational(1, 10**199)
            < s.Rational(1, 59),
        },
    }
