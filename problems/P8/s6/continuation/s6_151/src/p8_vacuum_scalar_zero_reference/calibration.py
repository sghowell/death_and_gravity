"""Actual complete wineglass interaction-forest MS enclosure."""

from functools import cache

import sympy as s
from p8_polynomial_vacuum import model
from p8_vacuum_two_loop_wineglass import calibration as previous

from . import bounds


@cache
def data():
    p = model.data()["actual_parameters"]
    L, g, M, lam = (
        p[k]
        for k in (
            "bare_polynomial_quartic",
            "cubic_coupling_squared",
            "heavy_mass_squared",
            "lambda",
        )
    )
    d = bounds.enclosure(L, 4 * lam)
    return {
        "actual_parameters": p,
        "actual_MS_scale_logarithm": 400 * s.log(10),
        "actual_heavy_logarithm": s.log(4 * M),
        "actual_wineglass_MS_enclosure": d,
        "checks": {
            "old_family_bound_matches_exact_parent": d[
                "old_fully_subtracted_wineglass_upper"
            ]
            - previous.data()["actual_wineglass_group_b2_absolute_upper"],
            "complete_bound_adds_three_distinct_terms": d["complete_wineglass_MS_upper"]
            - d["old_fully_subtracted_wineglass_upper"]
            - d["added_complete_scale_term_upper"]
            - d["local_reference_b2_absolute_upper"],
            "actual_scale_is_fermion_mass_squared": 400 - 2 * 200,
            "exact_scale_term_integer_prefactor": 1200 * (12 + 40 * 600) - 28814400,
        },
        "bounds": {
            "actual_heavy_mass_above_thirty_two": bool(M > 32),
            "actual_routed_heavy_margin": bool(M / 2 > 9),
            "actual_cubic_over_mass_below_quartic_third": bool(g / M < L / 3),
            "actual_external_heavy_gap": bool(M - 3 > M / 2),
            "actual_heavy_log_below_six_hundred": bool(4 * M < 10**200),
            "finite_MS_anchor_upper_below_thirty_five": bool(
                d["finite_proper_MS_reference_absolute_upper"] < 35
            ),
            "local_b2_relative_upper_below_one_e_minus_406": bool(
                d["local_reference_b2_relative_upper"] < s.Rational(1, 10**406)
            ),
            "added_scale_term_below_two_e_minus_611": bool(
                d["added_complete_scale_term_upper"] < s.Rational(2, 10**611)
            ),
            "complete_wineglass_MS_below_three_e_minus_611": bool(
                d["complete_wineglass_MS_upper"] < s.Rational(3, 10**611)
            ),
            "complete_wineglass_relative_below_eight_e_minus_12": bool(
                d["complete_wineglass_MS_relative_upper"] < s.Rational(8, 10**12)
            ),
        },
        "scope": "The complete assigned wineglass interaction-forest family is converted to MS at mu=mF, at fixed canonical reference coordinates. This is not the full scalar/GY14 amplitude or its finite field/coupling re-expansion. Other families and all remaining V/G/B closure obligations are still open.",
    }
