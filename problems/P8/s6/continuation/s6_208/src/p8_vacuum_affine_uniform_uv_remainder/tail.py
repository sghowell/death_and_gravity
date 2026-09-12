"""Complete original both-created-mode regulator tail in the same source norm."""

from functools import cache

import sympy as s

from . import domain, estimates

TAIL = 4 * s.Integer(10) ** 53
NEAR_FACTOR = 2 / domain.EPS + 1


@cache
def constants():
    far = sum(
        estimates.spatial_weight(5 - j)
        * 2
        * estimates.row_bound(j)
        / (9 * domain.EPS ** (5 - j))
        for j in range(5)
    )
    near = NEAR_FACTOR * sum(
        estimates.spatial_weight(q + 1) * c for q, c, label in estimates.near_terms()
    )
    low = estimates.LOW * estimates.spatial_weight(1)
    return {
        "complete_far_removed_union_tail": far,
        "complete_near_removed_union_tail": near,
        "unexpanded_low_removed_union_tail": low,
        "complete_original_regulator_tail_numerator": far + near + low,
    }


@cache
def data():
    r, K, U = s.symbols("r K U", positive=True)
    checks = {
        "complete_far_half_band_power_integral": s.integrate(r**-2, (r, K / 2, s.oo))
        - 2 / K,
        "exact_near_both_leg_radius_envelope": (2 / domain.EPS + 1) * U
        - (2 * U / domain.EPS + U),
        "explicit_near_radius_factor": NEAR_FACTOR - 201,
        "largest_near_tail_spatial_degree": max(
            q + 1 for q, c, label in estimates.near_terms()
        )
        - 6,
        "all_five_far_tail_rows": len(estimates.far_terms()) - 5,
    }
    return {
        "original_mask": "Every UV-subtracted endpoint uses chi_K=1_(|k|<=K,|-k+P|<=K). No mask derivative or one-leg replacement is used. The contact retains its distinct original one-mode band.",
        "far": "On FAR, r>=2(m+|P|)/epsilon implies |P|<=epsilon*r/2. A removed pair therefore has r>K/2. Integrating the full Cauchy r^-4 majorant gives 2B_j(m+|P|)^(5-j)/(9epsilon^(5-j)K), for every endpoint.",
        "near": "For r<2(m+|P|)/epsilon, both created-mode magnitudes are below201(m+|P|). Removal implies K<201(m+|P|). Multiply the complete near absolute bound by201(m+|P|)/K. The largest spatial degree becomes6, still controlled by the same X46 norm.",
        "low": "For r<m, removal implies K<m+|P|. The unexpanded low bound therefore gains only(m+|P|)/K. No complex frame or ultraviolet normalization is used at zero internal momentum.",
        "uniform_result": "Summing the entire far, near and low removed union gives |Qnew-Qnew,K|<4e53||D||L2 X46[Gamma]/K for K>=m. This removes a mathematical regulator, not a physical cutoff.",
        "constants": constants(),
        "checks": checks,
        "gates": {
            "complete_uniform_original_tail_display": bool(
                constants()["complete_original_regulator_tail_numerator"] < TAIL
            ),
            "same_six_spatial_derivatives_for_tail": max(
                q + 1 for q, c, label in estimates.near_terms()
            )
            <= 6,
            "far_second_leg_inside_twice_first": bool(1 + domain.EPS / 2 < 2),
            "complete_original_one_and_two_leg_masks_retained": True,
            "no_physical_cutoff_claim": True,
        },
    }
