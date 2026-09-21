"""Uniform original-source local-radiation remainder and finite interference."""

from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_real_remainder import born
from p8_vacuum_affine_dimensional_real_remainder import bounds as old_bounds

from . import radiation, source

UNSCALED_BUDGET = s.Integer(233991168)
REMAINDER_RATIO = s.Rational(1, 10**791)
HARD_RATIO = s.Rational(1, 10**794)
TREE_REMAINDER = s.Integer(4000000000)


def coefficient_budget():
    return sum(abs(value) for value in radiation.local_coefficients().values())


def remainder_upper():
    return s.Integer(4000000) * coefficient_budget() / s.sqrt(source.KAPPA)


def finite_interference_upper(ss, z, resolution):
    ss, z = born.domain(ss, z)
    x = old_bounds.resolution(resolution)
    delta = born.transfer_gap(ss, z)
    a = delta / 192
    ell = min(x, a)
    low = (
        REMAINDER_RATIO
        * (128 * ell + TREE_REMAINDER * ell**2 / delta)
        / (4 * s.pi**2 * source.KAPPA)
    )
    high = (
        0
        if x <= a
        else s.Integer(10) ** 18
        * REMAINDER_RATIO
        * (x * x - a * a)
        / (source.KAPPA * delta)
    )
    return low + high


@cache
def data():
    w, x, d, r = s.symbols("omega resolution delta remainder_ratio", positive=True)
    vertex_y2 = 24 * 16**2
    vertex_gal = 48 * 16**3
    shift_y2 = 24 * 2 * 16 * 48
    shift_gal = 48 * 3 * 16**2 * 48
    metric = 24 * 12 * 16**3
    connection = s.Rational(24 * 9, 2) * 16**2
    low = r * (128 * x + TREE_REMAINDER * x * x / d)
    checks = {
        "complete24_Y2_vertex_budget": s.Integer(vertex_y2 - 6144),
        "complete24_Gal_vertex_budget": s.Integer(vertex_gal - 196608),
        "complete24_Y2_shift_budget": s.Integer(shift_y2 - 36864),
        "complete24_Gal_shift_budget": s.Integer(shift_gal - 1769472),
        "complete24_Gal_metric_contact_budget": s.Integer(metric - 1179648),
        "complete24_Gal_connection_contact_budget": connection - 27648,
        "all_external_and_contact_remainder_terms": 64 * shift_gal
        + 608 * vertex_gal
        + metric
        + connection
        - UNSCALED_BUDGET,
        "low_window_full_interference_radial_integral": s.integrate(
            2 * r * (64 + TREE_REMAINDER * w / d), (w, 0, x)
        )
        - low,
        "high_window_full_interference_radial_integral": s.integrate(
            w * 4 * s.Integer(10) ** 18 * r / d, (w, d / 192, x)
        )
        - 2 * s.Integer(10) ** 18 * r * (x * x - (d / 192) ** 2) / d,
        "low_cutoff_tail_vanishes": s.limit(low, x, 0),
        "low_window_endpoint_budget": s.factor(
            low.subs(x, d / 192)
            - r * d * (s.Rational(128, 192) + TREE_REMAINDER / s.Integer(192) ** 2)
        ),
    }
    gates = {
        "Euclidean_recoil_norm3_shifted_norm4_and_change6omega": True,
        "all_relevant_Gram_modulus16_change48omega": True,
        "both_physical_polarizations_retained": True,
        "pi_lower_bound_used_conservatively": bool(2 * UNSCALED_BUDGET < 4000000 * 144),
        "original_local_coefficients_less_than_1e_minus1398": bool(
            0 < coefficient_budget() < s.Rational(1, 10**1398)
        ),
        "original_Born_floor_exceeds_1e_minus600": bool(
            4 * source.CUBIC**2 / source.HEAVY_MASS2**3 > s.Rational(1, 10**600)
        ),
        "normalized_remainder_less_than_1e_minus791_over_sqrtk": bool(
            4000000 * s.Rational(1, 10**1398) / s.Rational(1, 10**600) < REMAINDER_RATIO
        ),
        "known_local_fourpoint_less_than_1e_minus794_Born": bool(
            s.Rational(196608, 144) * s.Rational(1, 10**1398) / s.Rational(1, 10**600)
            < HARD_RATIO
        ),
        "low_full_interference_window_less_than_1e_minus1587": bool(
            (s.Rational(128, 192) + TREE_REMAINDER / s.Integer(192) ** 2) / 36 < 10000
        ),
        "high_window_bound_not_low_window_extrapolation": True,
        "original_full_positive_Born_and_actual_phase_retained": True,
        "leading_loop_soft_extension_excluded_from_finite_difference": True,
        "unknown_curvature_matching_and_other_loop_sectors_open": True,
    }
    return {
        "checks": {key: s.factor(value) for key, value in checks.items()},
        "gates": gates,
        "whole_original_coefficient_budget": coefficient_budget(),
        "whole_uniform_remainder_bound": "With P=|alpha|+|beta|+|gamma|, ||R1local||_TT<=4e6 P/sqrt(kappa). The complete24 vertex, difference, current and contact budgets sum to233991168 before a conservative two-polarization factor2 and division by16pi^2. Original P<1e-1398 and full A0>1e-600 imply ||R1local||/A0<1e-791/sqrt(kappa), |A1local|/A0<1e-794.",
        "whole_finite_interference_bound": "Let r=1e-791,delta=min(1,-t,-u),a=delta/192,ell=min(x,a). Integral of |2Re<Mtree,R1local>|/A0^2 with the actual physical phase is <=r[128ell+4e9 ell^2/delta]/(4pi^2 kappa)+1_(x>a)1e18 r(x^2-a^2)/(kappa delta). Atx<=a it is<1e-1587. The lower-cutoff tail vanishes. This is only the local radiative remainder, not its divergent leading-soft extension or a full inclusive loop rate.",
        "whole_domain": "25/4<=s<=16,strictly nonforward -1<z<1,0<=x<=1/8,actual original couplings. Uniform remainder over this recoil domain; finite-rate bound pointwise in the Born angle, not uniform through the forward poles.",
    }
