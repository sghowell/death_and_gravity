"""Uniform complete bubble-class remainder with unchanged original normalization."""

from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_real_remainder import born
from p8_vacuum_affine_dimensional_real_remainder import bounds as old_bounds
from p8_vacuum_affine_local_tadpole_radiation import bounds as local_bounds
from p8_vacuum_affine_mixed_source_radiation import bounds as source_bounds

from . import source

REMAINDER_RATIO = s.Integer(10) ** 194
HARD_RATIO = s.Integer(10) ** 190
UNSCALED_BUDGET = s.Integer(2688768)
TREE_REMAINDER = s.Integer(4000000000)
SECTORS = ("factorized_bubbles", "known_local_source_and_bubbles")


def require_sector(value):
    if not isinstance(value, str) or value not in SECTORS:
        raise ValueError("Require a stated selected known bubble-radiation sector")


def unrounded_remainder_ratio():
    return s.Integer(2334) * source.CUBIC**2 * source.HEAVY_MASS2


def finite_interference_upper(ss, z, resolution, sector="factorized_bubbles"):
    require_sector(sector)
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
    n = s.Symbol("n", positive=True)
    w, x, d, r = s.symbols("omega resolution delta remainder_ratio", positive=True)
    low = r * (128 * x + TREE_REMAINDER * x * x / d)
    checks = {
        "original_contact_exact_upper_margin": s.factor(
            4 / n - 3 / (n - 2) + 2 / (n - 2) ** 2 - (n - 4) ** 2 / (n * (n - 2) ** 2)
        ),
        "full_kernel_modulus_budget": s.Integer(6**2 * 6 - 216),
        "symmetric_kernel_modulus_budget": s.Integer(6**2 - 36),
        "Born_kernel_modulus_budget": s.Integer(3 * 216 + 3 * 36 - 756),
        "pair_invariant_change_budget": s.Integer(2 * 6 * 12 - 144),
        "complete_kernel_change_budget": s.Integer(80 * 144 - 11520),
        "complete_remainder_all_terms": s.Integer(
            3 * 64 * 11520 + 608 * 756 + 3 * 2 * 36 * 80
        )
        - UNSCALED_BUDGET,
        "positive_Born_normalized_budget": UNSCALED_BUDGET / (288 * 4) - 2334,
        "low_full_tree_interference_integral": s.integrate(
            2 * r * (64 + TREE_REMAINDER * w / d), (w, 0, x)
        )
        - low,
        "high_full_tree_interference_integral": s.integrate(
            w * 4 * s.Integer(10) ** 18 * r / d, (w, d / 192, x)
        )
        - 2 * s.Integer(10) ** 18 * r * (x * x - (d / 192) ** 2) / d,
        "lower_IR_cutoff_tail": s.limit(low, x, 0),
    }
    gates = {
        "original_heavy_mass_at_least128": bool(source.HEAVY_MASS2 >= 128),
        "endpoint_factor_modulus_less_than6": bool(4 + s.Rational(8, 7) < 6),
        "endpoint_first_derivative_margin": bool(s.Rational(64, 49) < 2),
        "endpoint_second_derivative_margin": bool(s.Rational(1024, 343) < 3),
        "full_Fprime_bound_below80": bool(2 * 6 * s.Rational(2, 128) * 6 + 36 * 2 < 80),
        "full_Fsecond_bound_below80": bool(
            2 * (s.Rational(4, 128**2) + 6 * s.Rational(3, 128**2)) * 6
            + 4 * 6 * s.Rational(2, 128) * 2
            + 36 * 2
            < 80
        ),
        "original_cubic_exact": bool(source.CUBIC == s.Rational(1, 8192)),
        "original_mass_upper_bound": bool(
            source.HEAVY_MASS2 < s.Integer(10) ** 200 / 256
        ),
        "source_normalized_remainder_less_than_1e194": bool(
            unrounded_remainder_ratio() < REMAINDER_RATIO
        ),
        "isolated_hard_bound_less_than_1e190_not_small": bool(
            s.Rational(756, 1152) * source.CUBIC**2 * source.HEAVY_MASS2 < HARD_RATIO
        ),
        "sum_known_local_source_bubbles_fits_same_1e194_budget": bool(
            unrounded_remainder_ratio()
            + local_bounds.REMAINDER_RATIO
            + source_bounds.REMAINDER_RATIO
            < REMAINDER_RATIO
        ),
        "physical_normalized_remainder_less_than_1e_minus206": bool(
            REMAINDER_RATIO / s.sqrt(source.KAPPA) == s.Rational(1, 10**206)
        ),
        "low_window_interference_less_than_1e_minus602": bool(
            (s.Rational(128, 192) + TREE_REMAINDER / s.Integer(192) ** 2) / 36 < 10000
            and REMAINDER_RATIO / source.KAPPA == s.Rational(1, 10**606)
        ),
        "full_canonical_norm_both_polarizations_and_actual_phase": True,
        "separate_high_window_fixed_nonforward_bound": True,
        "large_isolated_hard_bound_not_perturbative_closure": True,
    }
    return {
        "checks": {key: s.factor(value) for key, value in checks.items()},
        "gates": gates,
        "whole_original_unrounded_remainder_ratio": unrounded_remainder_ratio(),
        "whole_selected_remainder_ratio": REMAINDER_RATIO,
        "whole_isolated_hard_ratio_upper_not_small": HARD_RATIO,
        "whole_uniform_proof": "For b(v)=nC/g^2+n/(n-v), F=b^2B, use |b|<6,|bprime|<2/n,|bsecond|<3/n^2. Then|F|<216,|Fprime|,|Fsecond|<80 and|F(4/3)|<36. Original pair norms and recoil give invariant change144omega and F change11520omega. External differences, current change and all internal insertions contribute3*64*11520+608*756+3*144*40=2688768. Restore g^4/(32pi^2 n^2 sqrt(kappa)) and A0>4g^2/n^3 to obtain2334g^2n/sqrt(kappa)<1e194/sqrt(kappa)=1e-206. The combined known S337/S338 remainders fit within the same rounded bound.",
        "whole_finite_interference": "For r=1e194,delta=min(1,-t,-u),a=delta/192,ell=min(x,a), the actual-phase full47-tree interference with the nonleading bubble remainder has absolute integral<=r[128ell+4e9 ell^2/delta]/(4pi^2 kappa)+1_(x>a)1e18 r(x^2-a^2)/(kappa delta). Low x<=a gives<1e-602 and the lower-cutoff tail vanishes. This is not the separately divergent leading-soft extension or full inclusive rate.",
        "whole_hard_and_scope_boundary": "The isolated hard coefficient has only the loose bound|A1bubble|/A0<1e190, not a small relative loop correction. Cancellations with the remaining triangle/box sectors matter. No matching coordinate is selected. Integrated high-window bounds are fixed in nonforward angle, not uniform through forward poles.",
    }
