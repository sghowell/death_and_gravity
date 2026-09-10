"""Actual nonlocal quadratic enclosures at the unchanged GY14 reference."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_self_energy_chord import calibration as previous

from . import pole, regions


@cache
def data():
    p = previous.data()["actual_reference_parameters"]
    d = pole.enclosure(p["mF"], p["Y_upper"], p["a_upper"], 144)
    E = d["both_quadratic_primitives_soft_tail_upper"]
    B = d["nonlocal_on_shell_divided_remainder_upper"]
    return {
        "actual_reference_parameters": p,
        "enclosure": d,
        "checks": {
            "same_mass": p["mF"] - 10**200,
            "all_self_and_vertex_words_summed": E
            - d["self_energy_soft_tail_upper"]
            - d["vertex_soft_tail_upper"],
            "equal_conservative_self_and_vertex_bounds": d[
                "self_energy_soft_tail_upper"
            ]
            - d["vertex_soft_tail_upper"],
            "one_third_after_OS_Cauchy": 3 * B - E,
            "minimum_radius": d["minimum_Cauchy_radius"] - p["mF"] / 360,
        },
        "bounds": {
            "actual_mass_in_domain": bool(p["mF"] >= 720),
            "shared_positive_Y_enclosure": bool(0 < p["actual_Y"] < p["Y_upper"]),
            "shared_positive_gauge_enclosure": bool(0 < p["actual_a"] < p["a_upper"]),
            "self_rounded_prefactor_conservative": bool(
                regions.data()["exact_self_prefactor"] < 10**13
            ),
            "vertex_rounded_prefactor_conservative": bool(
                regions.data()["exact_vertex_prefactor"] < 3 * 10**13
            ),
            "both_soft_tails_strictly_positive": bool(E > 0),
            "actual_nonlocal_upper_below_one_e_minus_799": bool(
                0 < B < sp.Rational(1, 10**799)
            ),
        },
        "scope": "Nonlocal OS primitive correction only; no finite MS local references or full matched canonical pole.",
    }
