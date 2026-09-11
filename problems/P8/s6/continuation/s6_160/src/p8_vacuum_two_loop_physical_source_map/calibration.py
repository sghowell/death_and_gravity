"""Actual fixed map and inherited physical-source two-loop enclosures."""

from functools import cache

import sympy as s
from p8_offshell_vacuum import jets, matching
from p8_polynomial_vacuum import model
from p8_vacuum_full_heavy_source import calibration as source
from p8_vacuum_full_phi_normalization import calibration as pole
from p8_vacuum_full_two_loop_amplitude import calibration as amplitude
from p8_vacuum_full_vacuum_reference import calibration as vacuum


@cache
def data():
    p = model.data()["actual_parameters"]
    lam, gamma = p["lambda"], p["gamma"]
    m = matching.data()
    c = m["redundant_c"].subs({m["lambda"]: lam, m["gamma"]: gamma})
    R = jets.data()["cubic_field_redefinition"]
    variables = sorted(R.free_symbols.intersection(jets.BY_SYMBOL), key=str)
    terms = s.Poly(R, *variables).terms()
    degrees = {sum(powers) for powers, _ in terms}
    derivative_degrees = {
        sum(power * sum(jets.BY_SYMBOL[v]) for v, power in zip(variables, powers))
        for powers, _ in terms
    }
    h = s.Symbol("map_scaling")
    scaled = s.expand(R.subs({v: h * v for v in variables}, simultaneous=True))
    a = amplitude.data()["complete_canonical_two_loop_enclosure"]
    P = pole.data()["complete_canonical_one_plus_two_loop_OS_coefficient_upper"]
    V = vacuum.data()["complete_matched_two_loop_vacuum_absolute_upper"]
    J = source.data()["complete_second_H_source_absolute_upper"]
    return {
        "fixed_actual_map_coefficients": {"lambda": lam, "gamma": gamma, "c": c},
        "literal_R_jet_monomial_count": len(terms),
        "literal_R_field_degrees": sorted(degrees),
        "literal_R_total_derivative_degrees": sorted(derivative_degrees),
        "physical_source_b2_one_plus_two_loop_relative_upper": a["total_relative"],
        "physical_source_b2_two_loop_relative_upper": a["second_relative"],
        "physical_source_b2_formal_uniform_lower": a["formal_uniform_lower"],
        "physical_source_unit_disc_Phi_OS_coefficient_upper": P,
        "physical_source_complete_vacuum_reference_upper": V,
        "physical_source_complete_second_H_source_upper": J,
        "checks": {
            "same_fixed_actual_cubic_map_coefficient": s.factor(
                c - 4 * lam * lam / gamma
            ),
            "literal_map_homogeneous_degree_three": s.expand(scaled - h**3 * R),
            "literal_map_vanishes_at_zero_field": R.subs({v: 0 for v in variables}),
            "literal_map_has_only_cubic_field_degree": sum(n != 3 for n in degrees),
            "literal_map_derivative_order_at_most_four": sum(
                n > 4 for n in derivative_degrees
            ),
            "unchanged_parent_physical_source_b2_bound": a["total_relative"]
            - amplitude.data()["complete_canonical_two_loop_enclosure"][
                "total_relative"
            ],
            "unchanged_parent_physical_source_pole_bound": P
            - pole.data()["complete_canonical_one_plus_two_loop_OS_coefficient_upper"],
            "unchanged_parent_physical_source_H_reference": J
            - source.data()["complete_second_H_source_absolute_upper"],
        },
        "bounds": {
            "fixed_actual_map_coefficients_positive": bool(min(lam, gamma, c) > 0),
            "physical_source_complete_b2_relative_below_one_e_minus_6": bool(
                a["total_relative"] < s.Rational(1, 10**6)
            ),
            "physical_source_complete_second_b2_relative_below_one_e_minus_7": bool(
                a["second_relative"] < s.Rational(1, 10**7)
            ),
            "physical_source_formal_b2_strictly_positive": bool(
                a["formal_uniform_lower"] > 0
            ),
            "physical_source_unit_disc_pole_bound_below_one_e_minus_18": bool(
                P < s.Rational(1, 10**18)
            ),
            "physical_source_complete_vacuum_reference_below_one_e_595": bool(
                V < 10**595
            ),
            "physical_source_complete_second_H_source_below_one_e_189": bool(
                J < 10**189
            ),
        },
        "scope": "These are transferred correlators of the matched physical source J F(Psi), not new estimates for an independently minimally subtracted ordinary Psi. No ordinary-Psi second residue, global coordinate norm or omitted-order physical error is bounded here.",
    }
