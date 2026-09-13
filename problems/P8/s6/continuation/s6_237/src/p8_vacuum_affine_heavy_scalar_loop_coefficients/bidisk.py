"""Full joint complex bidisk and every omitted inverse-mass term."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_scalar_tree_matching import model

RADIUS = s.Rational(1, 4)
LOG_MAJORANT = s.Integer(500)
REMAINDER_CONSTANT = s.Integer(10) ** 12


@cache
def data():
    w, x, q = s.symbols("w channel ratio")
    actual = -3 * w / (1 - 2 * w) + 2 * w * w / (1 - 2 * w) ** 2 + w / (1 - x * w)
    polynomial = -2 * w + (x - 4) * w * w + (x * x - 4) * w**3 + x**3 * w**4
    tail = w**5 * (x**4 / (1 - x * w) + 16 / (1 - 2 * w) ** 2)
    ac = [0, 2, 7, 13, 27]
    wmax = s.Rational(1, 64)
    anorm = sum(ac[j] * wmax ** (j - 1) for j in range(1, 5))
    rnorm = 81 / (1 - 3 * wmax) + 16 / (1 - 2 * wmax) ** 2
    bnorm = sum(
        s.Rational(ac[i] * ac[j], 2) * wmax ** (i + j - 6)
        for i in range(1, 5)
        for j in range(1, 5)
        if i + j >= 6
    )
    tnorm = 2 * sum(
        ac[i] * 2 * 32 ** (j - 1) * wmax ** (i + j - 6)
        for i in range(1, 5)
        for j in range(1, 5)
        if i + j >= 6
    )
    weight = q**4 * (5 - 4 * q) / (1 - q) ** 2
    checks = {
        "complete_fixed_vertex_remainder_after_four_terms": s.cancel(
            actual - polynomial - tail
        ),
        "kinematic_channel_disk_radius": 2 + RADIUS + RADIUS / 2 - s.Rational(19, 8),
        "light_denominator_re_disk_lower": 1 - s.Rational(19, 32) - s.Rational(13, 32),
        "light_logarithm_complete_series_bound": s.Rational(19, 32)
        / (1 - s.Rational(19, 32))
        - s.Rational(19, 13),
        "uniform_inverse_mass_Cauchy_radius": 1 / (16 * s.Integer(2))
        - s.Rational(1, 32),
        "complete_weighted_tail_identity": s.cancel(
            q**4 / (1 - q) + q * s.diff(q**4 / (1 - q), q) - weight
        ),
        "complete_weighted_tail_endpoint": weight.subs(q, s.Rational(1, 2))
        - s.Rational(3, 4),
        "P_weighted_tail_prefactor": 2 * 12 - 24,
        "Q_weighted_tail_prefactor": 4 * 12 - 48,
        "triangle_full_tail_constant": 4 * 32**4 - 2**22,
        "box_full_tail_constant": 24 * 32**4 - 3 * 2**23,
        "entire_triangle_high_polynomial_constant": tnorm - 1003424,
        "entire_bubble_high_polynomial_constant": bnorm - s.Rational(2286169, 8192),
        "entire_a_polynomial_norm": anorm - s.Rational(553819, 262144),
        "entire_a_remainder_norm": rnorm - s.Rational(5981248, 58621),
        "b20_bidisk_Cauchy_factor": RADIUS**-2 - 16,
        "b21_bidisk_Cauchy_factor": RADIUS**-3 - 64,
        "b40_bidisk_Cauchy_factor": RADIUS**-4 - 256,
    }
    return {
        "joint_complex_kinematics": "v=s+t/2-2, |v|<=1/4, |t|<=1/4, s=2+v-t/2, u=2-v-t/2. This is a holomorphic subthreshold bidisk, not the real physical-angle window of S236.",
        "full_parameter_domain": "For each complex channel, |s|,|u|<=19/8<3, |L-1|<=19/32, Re L>=13/32>0, |L|<2 and |b|<2. The real parts of the original denominators are positive at the actual n. The complete parameter integrals and amplitude are jointly holomorphic on a neighborhood of the closed bidisk.",
        "extension_of_coefficient_disk": "Repeat the complete S236 absolute-value proof for COMPLEX L,b with |L|,|b|<=K=2. It depends only on these moduli. On |w|=1/32, |P|<2,|Q|<4 and all coefficient denominators/logarithms stay on their analytic branches. Only P,Q, not -log w, enter Cauchy at w0.",
        "full_logarithmic_majorant": "The principal Log L obeys |Log L|<=19/13<2. Since log n<462, M=500 bounds |log n-Log L|+3 and |B|. No physical-cut logarithm is crossed in this bidisk.",
        "complete_order_four_tail": "For u=32/n<=1/2, Ptail<=4u^4,Qtail<=8u^4. The full weighted series sum for k>=4 is u^4(5-4u)/(1-u)^2<=12u^4. Its derivative and the explicit -log w derivative give |R_C|<=2^22 M/n^5 and |R_D|<=3*2^23 M/n^6.",
        "full_vertex_four_term_polynomial": polynomial,
        "full_vertex_exact_tail": tail,
        "all_remaining_products": "For n>=64 the polynomial vertex norm plus its exact tail is below3/n and |r_a|<=128/n^5. The complete bubble error is below2048/n^6; the triangle pieces have constants6*2^22 M,1024 M,2^20 M; each of the two boxes has3*2^23 M, all at n^-6. No cross term is omitted.",
        "complete_actual_amplitude_remainder": "The full base first-loop amplitude minus its explicit n^-2 through n^-5 terms has modulus below10^12 g^4/(16pi^2 n^6) everywhere on the closed bidisk. The OS4 contact is constant and does not affect its coefficient derivatives.",
        "coefficient_remainder_factors": {"b20": 16, "b21": 64, "b40": 256},
        "boundary": "This analytic bound applies to the complete one-loop coefficient only. It does not estimate omitted loop orders or claim full exact-S-matrix analyticity.",
        "checks": {k: s.cancel(value) for k, value in checks.items()},
        "gates": {
            "complex_light_denominator_nonzero": s.Rational(13, 32) > 0,
            "complex_parameter_K_bound": 1 + s.Rational(19, 32) < 2,
            "complex_log_bound": s.Rational(19, 13) < 2,
            "uniform_log_majorant": 462 + 2 + 3 < LOG_MAJORANT,
            "actual_mass_in_full_product_domain": model.MASS2 > 64,
            "full_vertex_tail_majorant": rnorm < 128,
            "full_vertex_norm_including_tail": anorm + 128 * wmax**4 < 3,
            "full_bubble_polynomial_majorant": bnorm < 512,
            "full_triangle_polynomial_majorant": tnorm < 2**20,
            "full_bubble_all_parts_majorant": 2 * (512 + 384 + 1) < 2048,
            "full_channel_all_parts_majorant": 6 * 2**22
            + 1024
            + 2**20
            + 2 * 3 * 2**23
            + 2048
            < 2**27,
            "full_three_channel_majorant": 3 * 2**27 * LOG_MAJORANT
            < REMAINDER_CONSTANT,
            "joint_bidisk_not_inferred_from_real_window": True,
        },
    }
