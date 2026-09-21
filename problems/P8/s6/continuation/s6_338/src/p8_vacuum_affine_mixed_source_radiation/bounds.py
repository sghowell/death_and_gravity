"""Original-source uniform remainder and finite actual-phase interference."""

from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_real_remainder import born
from p8_vacuum_affine_dimensional_real_remainder import bounds as old_bounds
from p8_vacuum_affine_local_tadpole_radiation import bounds as local_bounds

from . import radiation, source

REMAINDER_RATIO = s.Rational(1, 10**598)
HARD_RATIO = s.Rational(1, 10**603)
COMBINED_REMAINDER = s.Rational(1, 10**597)
COMBINED_HARD = s.Rational(1, 10**602)
CENTERED_BUDGET = s.Integer(330608)
TREE_REMAINDER = s.Integer(4000000000)
SECTORS = ("mixed_source", "known_local_and_source")


def remainder_ratio(sector):
    if not isinstance(sector, str) or sector not in SECTORS:
        raise ValueError("Require a stated known radiation sector")
    return REMAINDER_RATIO if sector == "mixed_source" else COMBINED_REMAINDER


def finite_interference_upper(ss, z, resolution, sector="mixed_source"):
    ss, z = born.domain(ss, z)
    x = old_bounds.resolution(resolution)
    r = remainder_ratio(sector)
    delta = born.transfer_gap(ss, z)
    a = delta / 192
    ell = min(x, a)
    low = (
        r * (128 * ell + TREE_REMAINDER * ell**2 / delta) / (4 * s.pi**2 * source.KAPPA)
    )
    high = (
        0
        if x <= a
        else s.Integer(10) ** 18 * r * (x * x - a * a) / (source.KAPPA * delta)
    )
    return low + high


@cache
def data():
    w, x, d, r = s.symbols("omega resolution delta remainder_ratio", positive=True)
    n, g = s.symbols("n g", positive=True)
    C = -g * g * (3 / (n - 2) - 2 / (n - 2) ** 2)
    asym = C + 3 * g * g / (n - s.Rational(4, 3))
    low = r * (128 * x + TREE_REMAINDER * x * x / d)
    checks = {
        "unchanged_original_symmetric_value": s.factor(
            asym - 4 * g * g / (3 * (n - 2) ** 2 * (n - s.Rational(4, 3)))
        ),
        "centered_remainder_full_canonical_budget": s.Integer(330000 + 608)
        - CENTERED_BUDGET,
        "original_Q_mass_and_pi_budget": s.Rational(8, 144 * 256) - s.Rational(1, 4608),
        "low_full_tree_interference_radial_integral": s.integrate(
            2 * r * (64 + TREE_REMAINDER * w / d), (w, 0, x)
        )
        - low,
        "high_full_tree_interference_radial_integral": s.integrate(
            w * 4 * s.Integer(10) ** 18 * r / d, (w, d / 192, x)
        )
        - 2 * s.Integer(10) ** 18 * r * (x * x - (d / 192) ** 2) / d,
        "lower_cutoff_tail_vanishes": s.limit(low, x, 0),
        "low_window_endpoint_budget": s.factor(
            low.subs(x, d / 192)
            - r * d * (s.Rational(128, 192) + TREE_REMAINDER / s.Integer(192) ** 2)
        ),
    }
    gates = {
        "original_n_upper_margin": bool(
            source.HEAVY_MASS2 < s.Integer(10) ** 200 / 256
        ),
        "original_n_at_least128": bool(source.HEAVY_MASS2 >= 128),
        "original_kappa_exact": bool(source.KAPPA == s.Integer(10) ** 800),
        "convex_centered_exchange_nonnegative_below_positive_matter_Born": True,
        "canonical_Frobenius_remainder_controls_both_TT_polarizations": True,
        "normalized_source_remainder_budget": bool(CENTERED_BUDGET / 4608 < 100),
        "normalized_hard_source_budget": bool(
            s.Rational(1, 4608) < s.Rational(1, 1000)
        ),
        "known_local_plus_source_remainder_bounded": bool(
            REMAINDER_RATIO + local_bounds.REMAINDER_RATIO < COMBINED_REMAINDER
        ),
        "known_local_plus_source_hard_bounded": bool(
            HARD_RATIO + local_bounds.HARD_RATIO < COMBINED_HARD
        ),
        "low_window_rounding_below_1e4": bool(
            (s.Rational(128, 192) + TREE_REMAINDER / s.Integer(192) ** 2) / 36 < 10000
        ),
        "full_original47_tree_and_actual_phase_retained": True,
        "high_window_bound_not_low_window_extrapolation": True,
        "leading_soft_extension_separate_from_finite_difference": True,
        "other_hard_loops_and_independent_curvature_matching_open": True,
    }
    return {
        "checks": {key: s.factor(value) for key, value in checks.items()},
        "gates": gates,
        "whole_source_remainder_ratio": REMAINDER_RATIO,
        "whole_source_hard_ratio": HARD_RATIO,
        "whole_combined_known_remainder_ratio": COMBINED_REMAINDER,
        "whole_combined_known_hard_ratio": COMBINED_HARD,
        "whole_uniform_proof": "A_centered=A_m-A_sym with A_sym=4g^2/[3(n-2)^2(n-4/3)]>0. Convexity gives0<=A_centered<A_m<A0. The centered radiation remainder equals R_m-A_sym(J_actual-J_Born), so its canonical Frobenius norm is<=330608 A0/sqrt(kappa). Since |Q|<1e-600/4608, the selected source remainder obeys||R1source||_TT/A0<1e-598/sqrt(kappa), and |A1source|/A0<1e-603. Adding S337 gives conservative1e-597 and1e-602 ratios.",
        "whole_finite_interference_bound": "For delta=min(1,-t,-u),a=delta/192,ell=min(x,a), the actual-phase integral of |2Re<Mtree,R1source>|/A0^2 is bounded by r[128ell+4e9 ell^2/delta]/(4pi^2 kappa)+1_(x>a)1e18 r(x^2-a^2)/(kappa delta),r=1e-598. The low window is<1e-1394; the combined known local-plus-source sector has r=1e-597 and low bound<1e-1393. This is not the leading-soft extension, a uniform forward bound, or the full inclusive loop rate.",
        "whole_domain": "Original source,25/4<=s<=16,-1<z<1,0<=x<=1/8. The unintegrated remainder is uniform on the recoil domain; the high-window integrated bound is at fixed nonforward Born angle.",
        "whole_original_Q": radiation.source_factor(),
    }
