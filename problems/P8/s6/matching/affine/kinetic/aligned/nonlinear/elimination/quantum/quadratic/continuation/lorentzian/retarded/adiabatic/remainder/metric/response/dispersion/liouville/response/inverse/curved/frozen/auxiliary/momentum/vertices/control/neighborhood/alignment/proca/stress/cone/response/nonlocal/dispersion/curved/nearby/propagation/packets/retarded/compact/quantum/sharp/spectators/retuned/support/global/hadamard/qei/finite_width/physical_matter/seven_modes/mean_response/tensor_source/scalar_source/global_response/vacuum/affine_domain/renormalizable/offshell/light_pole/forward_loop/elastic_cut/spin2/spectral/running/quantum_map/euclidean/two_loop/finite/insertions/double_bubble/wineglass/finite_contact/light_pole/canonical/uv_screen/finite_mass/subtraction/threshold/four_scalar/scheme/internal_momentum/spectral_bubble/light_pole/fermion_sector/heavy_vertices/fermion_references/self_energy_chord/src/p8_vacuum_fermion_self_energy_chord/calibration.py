"""Actual separate scalar/gauge self-energy-chord coefficient enclosures."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model
from p8_vacuum_fermion_local_matching import calibration as fermion
from p8_vacuum_gauge_yukawa_screen import calibration as candidate

from . import tail


@cache
def data():
    f = fermion.data()
    c = candidate.data()
    p = model.data()["actual_parameters"]
    m = f["fermion_mass"]
    Yhi = f["rational_Yukawa_squared_upper"]
    ahi = 2 * p["bare_polynomial_quartic"] / 3
    result = tail.enclosure(m, Yhi, ahi, 144)
    E = result["combined_self_energy_chord_b2_absolute_upper"]
    tree = 4 * p["lambda"]
    relative = E / tree
    return {
        "actual_reference_parameters": {
            "mF": m,
            "Y_upper": Yhi,
            "a_upper": ahi,
            "actual_Y": f["actual_Yukawa_squared"],
            "actual_a": c["prospective_gauge_squared"],
            "tree_b2": tree,
        },
        "actual_self_energy_chord_enclosure": result,
        "combined_self_energy_chord_upper_relative_to_tree": relative,
        "checks": {
            "same_actual_fermion_mass": m - 10**200,
            "same_actual_tree_reference": tree - 4 * sp.Rational(1, 10**600),
            "scalar_and_gauge_subsets_sum": E
            - result["scalar_self_energy_chord_b2_absolute_upper"]
            - result["gauge_self_energy_chord_b2_absolute_upper"],
            "relative_tree_dictionary": relative * tree - E,
            "actual_Cauchy_radius_minimum": result["minimum_soft_scaling_Cauchy_radius"]
            - m / 360,
            "same_four_remaining_graph_groups_not_completed": 2
            * result["unbounded_words_per_sector"]
            - 2 * (24 + 12),
        },
        "bounds": {
            "actual_mass_satisfies_uniform_soft_tail_domain": bool(m >= 720),
            "actual_Y_below_shared_rational_upper": bool(
                0 < f["actual_Yukawa_squared"] < Yhi
            ),
            "actual_a_below_shared_rational_upper": bool(
                0 < c["prospective_gauge_squared"] < ahi
            ),
            "exact_radial_prefactor_below_rounded_bound": tail.data()[
                "rounded_prefactor_is_conservative"
            ],
            "actual_scalar_subset_bound_strictly_positive": bool(
                result["scalar_self_energy_chord_b2_absolute_upper"] > 0
            ),
            "actual_gauge_subset_bound_strictly_positive": bool(
                result["gauge_self_energy_chord_b2_absolute_upper"] > 0
            ),
            "combined_chord_bound_below_one_e_minus_1400": bool(
                0 < E < sp.Rational(1, 10**1400)
            ),
            "relative_chord_bound_below_one_e_minus_800": bool(
                0 < relative < sp.Rational(1, 10**800)
            ),
        },
        "scope": "Two partial primitive rows, in the shared proper MS scheme. Positivity of a majorant is not a sign statement for the amplitude. Remaining primitive words, matching conversions, full two-loop error and original P8 remain open.",
    }
