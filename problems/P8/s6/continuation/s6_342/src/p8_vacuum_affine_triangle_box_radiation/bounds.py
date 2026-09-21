"""Uniform full selected radiative remainder and actual full-tree interference."""

from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_real_remainder import born
from p8_vacuum_affine_dimensional_real_remainder import bounds as tree_bounds
from p8_vacuum_affine_factorized_bubble_radiation import bounds as known
from p8_vacuum_affine_physical_loop_contours import insertions, moments

from . import source

REMAINDER_RATIO = s.Integer(10) ** 204
HARD_RATIO = s.Integer(10) ** 197
HARD_BUDGET = s.Integer(205824000)
CHANGE_BUDGET = s.Integer(263738720256000)
REMAINDER_BUDGET = s.Integer(17515528820736000)
SECTORS = ("triangle_box", "all_current_known_matter_sectors")


def require_sector(value):
    if not isinstance(value, str) or value not in SECTORS:
        raise ValueError("Require a stated selected massive-matter radiative sector")
    return value


def unrounded_remainder_ratio():
    return REMAINDER_BUDGET * source.CUBIC**2 * source.HEAVY_MASS2 / 576


def finite_interference_upper(ss, z, resolution, sector="triangle_box"):
    require_sector(sector)
    ss, z = born.domain(ss, z)
    x = tree_bounds.resolution(resolution)
    delta = born.transfer_gap(ss, z)
    a = delta / 192
    ell = min(x, a)
    low = (
        REMAINDER_RATIO
        * (128 * ell + s.Integer(4 * 10**9) * ell**2 / delta)
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
    checks = {}
    gates = {}
    Bc = s.Integer(256000)
    Bd = s.Integer(32768000)
    hard = 36 * Bc + 6 * Bd
    change = 6 * Bc * (6 * 160 * 4096 + 576) + 6 * Bd * 320 * 4096
    budget = 64 * change + 1216 * hard + 36 * 10**12 + 6 * 10**14 + 1728 * Bc
    gates["frozen_scalar_contour_budgets"] = moments.MASTER_BUDGETS == {
        "triangle": Bc,
        "ordered_box": Bd,
    }
    gates["frozen_internal_triangle_budget"] = (
        2 * insertions.BUDGETS["triangle_light"]
        + insertions.BUDGETS["triangle_heavy"] / source.HEAVY_MASS2
    ) * 10000 < 10**12
    gates["frozen_internal_box_budget"] = (
        2 * insertions.BUDGETS["box_light"]
        + 2 * insertions.BUDGETS["box_heavy"] / source.HEAVY_MASS2
    ) * 10000 < 10**14
    gates["pair_shift_change_margin"] = 10 * (7 + 6) < 144
    gates["singleton_virtuality_change"] = 2 * (2 + s.sqrt(3)) < 8
    gates["hard_rough_bound_not_small"] = (
        hard * source.CUBIC**2 * source.HEAVY_MASS2 / 288 < 10**197
    )
    unrounded = budget * source.CUBIC**2 * source.HEAVY_MASS2 / 576
    gates["original_TD_remainder_ratio_less_than1e204"] = unrounded < 10**204
    gates["combined_known_matter_remainder_fits_same1e204"] = (
        unrounded + known.REMAINDER_RATIO < 10**204
    )
    gates["actual_normalized_finite_remainder_1e_minus196"] = s.Integer(
        10
    ) ** 204 / s.sqrt(source.KAPPA) == s.Rational(1, 10**196)
    gates["low_interference_1e_minus592"] = (
        s.Rational(128, 192) + s.Integer(4 * 10**9) / 192**2
    ) / 36 < 10000 and s.Integer(10) ** 204 / source.KAPPA == s.Rational(1, 10**596)
    checks["complete_hard_scalar_budget"] = hard - 205824000
    # Exact continuous external-emission budgets, retaining each shifted leg.
    R = s.Rational
    gates["individual_TT_currents_sum_norm_bound"] = 4 * 4 * 4 == 64
    gates["triangle_virtuality_and_pair_path_length"] = 144 + 2 * 8 == 160
    gates["box_all_six_coordinate_path_length"] = 2 * 144 + 4 * 8 == 320
    gates["two_heavy_resolvents_outer_pair_budget"] = 24 * R(1, 4) * 2 * 36 * 4 == 1728
    checks["triangle_A_times_C_change_all24_labels"] = 6 * Bc * (
        6 * 160 * 4096 + 576
    ) - 24 * R(1, 4) * Bc * (6 * 160 * 4096 + 576)
    checks["box_change_all24_labels"] = (
        6 * Bd * 320 * 4096 - 24 * R(1, 4) * Bd * 320 * 4096
    )
    checks["all_triangle_internal_lines_and_endpoint_prefactor"] = (
        36 * 10**12 - 24 * R(1, 4) * 6 * 10**12
    )
    checks["all_box_internal_lines_and_endpoint_prefactor"] = (
        6 * 10**14 - 24 * R(1, 4) * 10**14
    )
    ww, ell, delta, rr = s.symbols("omega ell delta remainder_ratio", positive=True)
    checks["full_tree_low_finite_interference_integral"] = s.integrate(
        2 * rr * (64 + 4 * 10**9 * ww / delta), (ww, 0, ell)
    ) - rr * (128 * ell + 4 * 10**9 * ell**2 / delta)
    checks["lower_cutoff_nonleading_tail_vanishes"] = s.limit(
        rr * (128 * ell + 4 * 10**9 * ell**2 / delta), ell, 0
    )

    checks["whole_remainder_budget_exact"] = budget - REMAINDER_BUDGET
    checks["whole_change_budget_exact"] = change - CHANGE_BUDGET
    return {
        "checks": {key: s.factor(value) for key, value in checks.items()},
        "gates": {key: bool(value) for key, value in gates.items()},
        "whole_hard_budget": HARD_BUDGET,
        "whole_change_budget": CHANGE_BUDGET,
        "whole_unscaled_remainder_budget": REMAINDER_BUDGET,
        "whole_unrounded_original_remainder_ratio": unrounded_remainder_ratio(),
        "whole_rounded_original_remainder_ratio": REMAINDER_RATIO,
        "whole_isolated_hard_ratio_upper_not_small": HARD_RATIO,
        "whole_continuous_estimate": "With Bc256000,Bd32768000: M=36Bc+6Bd; L=6Bc*(6*160*4096+576)+6Bd*320*4096. External coefficient differences, current variation, all triangle lines, all box lines and the outer-heavy branch give B=64L+1216M+36e12+6e14+1728Bc=17515528820736000. Restore g^4/(16pi^2 n^2 sqrt(kappa)) and A0>4g^2/n^3: B*g^2*n/576<1e204. Adding the already known S337/S338/S339 remainder bound and zero S340 still fits1e204.",
        "whole_finite_interference": "For delta=min(1,-t,-u),a=delta/192,ell=min(x,a), the actual-phase integral of|2Re<Mtree,R>|/A0^2 is <=r[128ell+4e9ell^2/delta]/(4pi^2 kappa)+1_(x>a)1e18*r*(x^2-a^2)/(kappa*delta),r=1e204. Low x<=a is<1e-592 and its lower-cutoff tail vanishes. The high bound is finite at every fixed nonforward angle, not uniform through a forward pole.",
        "whole_large_hard_boundary": "The isolated hard triangle/box coefficient has only the loose bound<1e197 relative to A0, not perturbative smallness. S236's sharper total flat cancellation result is separate. The independently divergent leading-soft loop extension, virtual hard interference, all unknown curvature/internal-gravity terms and higher orders are excluded from the finite-remainder interference claim.",
    }
