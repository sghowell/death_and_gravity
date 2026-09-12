"""Explicit projection norms and bounded known one-point Ward terms."""

from functools import cache

import sympy as s
from p8_vacuum_affine_matched_spatial_current import assembly

DETECTOR = s.Integer(100)
SOURCE = s.Integer(10) ** 20
WARD = s.Integer(10) ** 35
KNOWN = s.Integer(10) ** 118


def leibniz_bound(order, leading):
    if type(order) is not int or not 0 <= order <= 13:
        raise ValueError("Only the specified native derivative orders zero through13")
    return sum(
        s.binomial(order, j) * leading * s.factorial(j) * 4**j for j in range(order + 1)
    )


@cache
def data():
    h13 = leibniz_bound(13, 7)
    c12 = leibniz_bound(12, 28)
    source_raw = 4 * (3 + 4 * h13 + 2 * c12)
    detector_raw = 1 + 4 * 7 + 2 * (1 + 28)
    source_ward = 6 * 84 + 30 * 4 * 13
    detector_ward = 4 * (3 * 376 + 2 * 4 * 52) + 30 * 4 * 13 * 13
    ward_raw = (source_ward + detector_ward) * s.Integer(10) ** 30
    v = s.Symbol("momentum_squared", nonnegative=True)
    checks = {
        "complex_disk_denominator_lower": 1 - s.Rational(3, 4) ** 2 - s.Rational(7, 16),
        "Cauchy_H_bound_exact": 4 * s.Rational(3, 4) / s.Rational(7, 16)
        - s.Rational(48, 7),
        "source_exact_rounded_predecessor": source_raw
        - s.Integer(62408287126207901212),
        "detector_exact_rounded_predecessor": detector_raw - 87,
        "source_Ward_terms_both_retained": source_ward - 2064,
        "detector_Ward_terms_both_retained": detector_ward - 26456,
        "full_one_point_correction_unrounded": ward_raw - 28520 * s.Integer(10) ** 30,
        "full_known_tracefree_projection_factor": assembly.CURRENT * DETECTOR * SOURCE
        - 2 * s.Integer(10) ** 117,
        "source_two_extra_spatial_derivatives_weight": s.expand(
            (1 + v) ** 8 - (1 + v) ** 6 * v**2
        )
        - s.expand((1 + v) ** 6 * (1 + 2 * v)),
        "detector_safe_four_spatial_derivatives_weight": s.expand(
            (1 + v) ** 4 - (1 + v) * v**2
        )
        - s.expand((1 + v) * (1 + 3 * v + 2 * v * v + v**3)),
    }
    return {
        "input_norms": "V04[D]^2=integral_I ||(1-Delta)^2(n_D,beta_D,Q_D)||L2_x^2 dt. U138[G]^2=sum_r0..13 integral_I ||(1-Delta)^4 partial_t^r(n_G,beta_G,Q_G)||L2_x^2 dt. Component tuples use Euclidean/Frobenius Hilbert norms and the unitary Fourier transform.",
        "time_primitive_bound": "Both Iminus and Iplus have L2_time norm at most1 on the unit slab. For j>=1, partial_t^j I f=partial_t^(j-1)f. This derivative identity, spatial commutation and Leibniz's rule control all source jets without an inverse H or momentum.",
        "Cauchy_coefficients": "On |z-t|<=1/4 for real |t|<=1/2, |1+z^2|>=7/16, |H(z)|<7 and |a(z)^-2|<28. Hence H_j<=7j!4^j and c_j<=28j!4^j. With B_r(b)=sum binom(r,j)b j!4^j, the source coefficient is4[3+4B13(7)+2B12(28)].",
        "projection_unrounded": {
            "detector": detector_raw,
            "source": source_raw,
            "B13_H": h13,
            "B12_c": c12,
        },
        "projection_result": "M[Qsyn_D]<=100 V04[D] and Z136[Qsyn_G]<=1e20 U138[G]. Orthogonal trace/tracefree and scalar projections do not enlarge these bounds. The source has the initial-only extension of the proved tracefree domain described in notes/boundary.md.",
        "one_point_coefficients": "On the real slab a<2, |H|<=2, |H'|<=4, a^-2<=1. The original |rho|,|P| and their first derivatives are<1e30. Thus ||E||F<4e30 and ||E'||F<2e31. Both gauge vectors have L2 norm<=3 times the appropriate input and their full derivative matrices<=4 times it; their ADM gauge directions have norm<=13 times it.",
        "one_point_bound_derivation": "The source Lie-density term is bounded by6*84e30 and its chart term by30*4*13e30. For synchronous G, hG, spatial grad hG and hG' have bounds52,52,324 timesU; use376 for the full derivative norm. The detector Lie term is4(3*376+2*4*52)e30 and its chart is30*4*13^2 e30. The sum28520e30 is below1e35.",
        "known_piece": "The already matched tracefree contribution plus the complete one-point Ward correction is bounded by1e118 V04 U138. This bounds only known pieces, not the missing scalar response.",
        "conditional_full_bound": "If the three ordered scalar forms satisfy |Rab(d,g)|<=Lab M[d]Z136[g] in the same initial-only domain and fixed prescription, then |Rfull(D,G)|<=[1e22(2e95+Ltt+Lt0+L0t)+1e35] V04[D]U138[G]. This is an implication, not a proved finite value of any Lab.",
        "boundary": "These are unreduced metric coordinate norms with derivative loss. No canonical reduced scalar norm, same-space inverse, contraction, finite-amplitude or quantum background/stability result follows.",
        "checks": checks,
        "gates": {
            "complex_H_majorant_strict": s.Rational(48, 7) < 7,
            "complex_inverse_scale_majorant_strict": s.Rational(16, 7) ** 4 < 28,
            "source_projection_strict_rounding": source_raw < SOURCE,
            "detector_projection_strict_rounding": detector_raw < DETECTOR,
            "one_point_Ward_strict_rounding": ward_raw < WARD,
            "known_piece_strict_rounding": assembly.CURRENT * DETECTOR * SOURCE + WARD
            < KNOWN,
            "derivative_loss_not_canonical_reduced_inverse": True,
            "conditional_full_bound_not_scalar_existence": True,
        },
    }
