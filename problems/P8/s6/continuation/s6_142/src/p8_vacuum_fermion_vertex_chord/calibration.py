"""Actual complete scalar/gauge quartic primitive bounds, with disjoint subsets."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_opposite_chord import calibration as opposite
from p8_vacuum_fermion_self_energy_chord import calibration as self_energy

from . import tail


@cache
def data():
    old = self_energy.data()
    other = opposite.data()
    p = old["actual_reference_parameters"]
    new = tail.enclosure(p["mF"], p["Y_upper"], p["a_upper"], 144)
    self_bound = old["actual_self_energy_chord_enclosure"]
    opp_bound = other["actual_opposite_chord_enclosure"]
    scalar = (
        self_bound["scalar_self_energy_chord_b2_absolute_upper"]
        + opp_bound["scalar_opposite_chord_b2_absolute_upper"]
        + new["scalar_vertex_chord_b2_absolute_upper"]
    )
    gauge = (
        self_bound["gauge_self_energy_chord_b2_absolute_upper"]
        + opp_bound["gauge_opposite_chord_b2_absolute_upper"]
        + new["gauge_vertex_chord_b2_absolute_upper"]
    )
    vertex = new["combined_vertex_chord_b2_absolute_upper"]
    total = scalar + gauge
    return {
        "actual_reference_parameters": p,
        "actual_vertex_chord_enclosure": new,
        "complete_scalar_quartic_primitive_upper": scalar,
        "complete_gauge_quartic_primitive_upper": gauge,
        "both_complete_quartic_primitive_upper": total,
        "both_primitive_upper_relative_to_tree": total / p["tree_b2"],
        "checks": {
            "same_shared_mass": p["mF"] - 10**200,
            "same_tree_coefficient": p["tree_b2"] - 4 * sp.Rational(1, 10**600),
            "same_rounded_coupling_dictionary": vertex
            - 5 * opp_bound["combined_opposite_chord_b2_absolute_upper"],
            "separate_scalar_and_gauge_sum": total - scalar - gauge,
            "three_disjoint_subsets_sum": total
            - vertex
            - opp_bound["combined_opposite_chord_b2_absolute_upper"]
            - self_bound["combined_self_energy_chord_b2_absolute_upper"],
            "complete_word_count_per_sector": 24
            + 12
            + new["newly_bounded_words_per_sector"]
            - 60,
        },
        "bounds": {
            "mass_in_joint_domain": bool(p["mF"] >= 720),
            "shared_positive_Y_upper": bool(0 < p["actual_Y"] < p["Y_upper"]),
            "shared_positive_a_upper": bool(0 < p["actual_a"] < p["a_upper"]),
            "rounded_vertex_prefactor_conservative": tail.data()[
                "rounded_prefactor_is_conservative"
            ],
            "vertex_bound_positive": bool(vertex > 0),
            "scalar_primitive_bound_positive": bool(scalar > 0),
            "gauge_primitive_bound_positive": bool(gauge > 0),
            "both_complete_primitive_bound_below_one_e_minus_1400": bool(
                0 < total < sp.Rational(1, 10**1400)
            ),
            "relative_bound_below_one_e_minus_800": bool(
                0 < total / p["tree_b2"] < sp.Rational(1, 10**800)
            ),
        },
        "scope": "Two complete paired quartic primitive rows only. Other primitive rows, finite field/parameter conversions and the full matched two-loop error remain open.",
    }
