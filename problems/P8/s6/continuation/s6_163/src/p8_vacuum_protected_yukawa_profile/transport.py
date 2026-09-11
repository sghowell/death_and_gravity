"""Exactly unchanged named low-order observables in the SAT8 EFT branch."""

from functools import cache

import sympy as s
from p8_vacuum_full_analytic_target_match import transport as classical
from p8_vacuum_two_loop_elastic_cut import calibration as cut
from p8_vacuum_two_loop_physical_source_map import calibration as physical

from . import order


def rows():
    return [
        {
            "observable": "physical_Phi_pole_and_residue",
            "source_endpoints": 2,
            "through_loop": 2,
        },
        {
            "observable": "physical_Phi_four_point_amplitude_and_b2",
            "source_endpoints": 4,
            "through_loop": 2,
        },
        {"observable": "vacuum_reference", "source_endpoints": 0, "through_loop": 2},
        {
            "observable": "stationary_H_source_reference",
            "source_endpoints": 1,
            "through_loop": 2,
        },
        {
            "observable": "physical_Phi_first_amplitude_for_two_loop_cut",
            "source_endpoints": 4,
            "through_loop": 1,
        },
    ]


def validate_rows(observations):
    if observations != rows():
        raise ValueError(
            "The SAT8 low-order transport scope differs from the fixed observables"
        )
    return True


@cache
def data():
    p = physical.data()
    c = cut.data()
    t = classical.data()
    preserved = {
        "total_b2_relative_upper": p[
            "physical_source_b2_one_plus_two_loop_relative_upper"
        ],
        "second_b2_relative_upper": p["physical_source_b2_two_loop_relative_upper"],
        "unit_disc_Phi_pole_coefficient_upper": p[
            "physical_source_unit_disc_Phi_OS_coefficient_upper"
        ],
        "second_vacuum_reference_upper": p[
            "physical_source_complete_vacuum_reference_upper"
        ],
        "second_H_source_upper": p["physical_source_complete_second_H_source_upper"],
        "complete_first_physical_amplitude_upper": c[
            "complete_one_loop_physical_amplitude_upper"
        ],
        "second_elastic_cut_relative_upper": c[
            "complete_second_elastic_cut_tree_relative_upper"
        ],
        "complete_cut_subtracted_formal_relative_upper": c[
            "complete_through_two_loop_improved_formal_band"
        ]["total_relative_error"],
        "classical_full_analytic_target_match_coefficient": t[
            "full_analytic_target_stationary_action_match_coefficient"
        ],
    }
    checks = {
        "all_named_changes_first_possible_after_stated_loop_orders": sum(
            order.new_vertex_floor(row["source_endpoints"]) <= row["through_loop"]
            for row in rows()
        ),
        "same_complete_b2_allowance": preserved["total_b2_relative_upper"]
        - p["physical_source_b2_one_plus_two_loop_relative_upper"],
        "same_second_elastic_cut_allowance": preserved[
            "second_elastic_cut_relative_upper"
        ]
        - c["complete_second_elastic_cut_tree_relative_upper"],
        "same_classical_full_target_error": preserved[
            "classical_full_analytic_target_match_coefficient"
        ]
        - t["full_analytic_target_stationary_action_match_coefficient"],
        "both_active_and_all_spectator_flavors": 2 + 12 - 14,
        "active_color_flavor_multiplicity": 2 * 3 - 6,
    }
    return {
        "transported_observables": rows(),
        "identical_named_coefficient_bounds": preserved,
        "fixed_reference": "The same lower-field GY14 MS interaction values, physical Phi mass/residue and vacuum/H references are imposed at the same scale, with the same regulator before finite parts. New higher-field EFT counterterms retain their loop grades. No finite lower-field retuning is used to assert this equality.",
        "physical_source": "Apply the same full local cubic F and source J F(Psi) to the SAT8 action. The new difference f_R(F)-F starts at scalar degree nine and cannot enter these named graphs through two loops. This is not equality for arbitrary numbers of physical sources, an ordinary-Psi MS composite normalization or a global inverse.",
        "classical_zero_fermion_identity": "The action change is proportional to fermion bilinears and vanishes identically when the classical fermions are zero. The full flat classical target match therefore stays exactly as in S6.162, not merely approximately unchanged by the new profile.",
        "bounds": {
            "all_named_observables_protected_through_declared_orders": all(
                order.new_vertex_floor(row["source_endpoints"]) > row["through_loop"]
                for row in rows()
            ),
            "same_total_b2_relative_below_one_e_minus_6": bool(
                preserved["total_b2_relative_upper"] < s.Rational(1, 10**6)
            ),
            "same_second_b2_relative_below_one_e_minus_7": bool(
                preserved["second_b2_relative_upper"] < s.Rational(1, 10**7)
            ),
            "same_unit_disc_Phi_pole_below_one_e_minus_18": bool(
                preserved["unit_disc_Phi_pole_coefficient_upper"]
                < s.Rational(1, 10**18)
            ),
            "same_second_vacuum_reference_below_one_e_595": bool(
                preserved["second_vacuum_reference_upper"] < 10**595
            ),
            "same_second_H_source_below_one_e_189": bool(
                preserved["second_H_source_upper"] < 10**189
            ),
            "same_second_cut_relative_below_one_e_minus_407": bool(
                preserved["second_elastic_cut_relative_upper"] < s.Rational(1, 10**407)
            ),
            "same_improved_formal_band_below_one_e_minus_6": bool(
                preserved["complete_cut_subtracted_formal_relative_upper"]
                < s.Rational(1, 10**6)
            ),
            "same_classical_full_target_match_below_one_e_minus_800": bool(
                preserved["classical_full_analytic_target_match_coefficient"]
                < s.Rational(1, 10**800)
            ),
        },
        "checks": checks,
        "not_inferred": "A real mass floor does not establish a transported rolling state, common-parent B or an interacting cutoff. The new nonpolynomial EFT is not the original globally renormalizable GY14 model; its all-field UV running and high-energy behavior are not inherited from the old marginal ray. Physical remainders and V/G/B remain open.",
    }
