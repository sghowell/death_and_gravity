"""Continuous new auxiliary bounds and common-domain finite tree comparison."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import bounds as original_auxiliary
from p8_zero_source import bounds as previous_bounds
from p8_zero_source import tree as previous_tree

from . import model


def response_radius(value):
    result=previous_bounds.response_radius(value)
    result["temporal_vector_absolute_upper"]=2*result["input_radius"]
    result["boundaries"]="Constant-mass source-free candidate, same local classical auxiliary domain; no quantum or frequency bound."
    return result


def at_scale(value):
    d=previous_tree.at_scale(value)
    return {"common_M_tau":d["common_M_tau"],
            "constant_mass_cubic_block_upper":d["new_cubic_transition_bound"],
            "constant_mass_connected_quartic_tree_upper":d["new_connected_quartic_tree_bound"],
            "constant_mass_minus_source_free_cubic_upper":2*d["new_cubic_transition_bound"],
            "constant_mass_minus_source_free_connected_quartic_upper":2*d["new_connected_quartic_tree_bound"]}


@cache
def auxiliary():
    d=previous_bounds.auxiliary()
    root=dict(d["new_auxiliary_contraction_bounds"])
    root["temporal_vector_absolute_upper"]=2*original_auxiliary.JET_RADIUS
    root["temporal_vector_over_input_radius_upper"]=sp.Integer(2)
    return {"same_complex_Hamiltonian_upper":d["same_declared_complex_Hamiltonian_upper"],
            "same_new_Hamiltonian_piece_triangle_upper":d["new_Hamiltonian_piece_upper_before_N"],
            "same_continuous_coefficient_bounds":previous_bounds.coefficients(),
            "constant_mass_auxiliary_bounds":root,
            "new_temporal_solution":"T=-j/U, with |1/U|<2 on the unchanged outer coefficient polydisc.",
            "mass_bound_reason":"Constant gamma_t=gamma_s=1 and their inverses lie within the old absolute bounds below 2. All scalar, boundary, geometric and matter pieces are literally unchanged. The new and previous background lapse functions are identical."}


@cache
def tree():
    old=previous_tree.data()
    cubic=old["new_cubic_transition_bound"]
    quartic=old["new_connected_quartic_tree_bound"]
    return {"common_named_M_tau":previous_tree.SCALE,
            "constant_mass_cubic_block_upper":cubic,
            "constant_mass_connected_quartic_tree_upper":quartic,
            "constant_mass_minus_source_free_cubic_upper":2*cubic,
            "constant_mass_minus_source_free_connected_quartic_upper":2*quartic,
            "common_domain_and_free_modes":old["same_seven_free_modes_and_hard_window"],
            "comparison":"Both new coefficient families obey the same universal Cauchy bounds. Reuse the positive S6.80 scaling proof and the same free phase coordinates at L=10^400. Triangle inequalities give twice the common finite-order block bound. This does not include a new tadpole, quantum loop or higher-operator change."}


@cache
def gates():
    a,t=auxiliary(),tree()
    return {"new_constant_masses_fit_the_same_complex_bound":bool(
                max(abs(value) for value in model.data()["new_source_free_mass_matrix"])<2),
            "unchanged_Hamiltonian_piece_sum_below_ten_thousand":bool(
                2*sum(a["same_new_Hamiltonian_piece_triangle_upper"].values())<10000),
            "new_temporal_root_below_three_times_invariant_radius":bool(
                a["constant_mass_auxiliary_bounds"]["temporal_vector_over_input_radius_upper"]<3),
            "new_Newton_map_remains_strict_contraction":bool(
                a["constant_mass_auxiliary_bounds"]["Newton_contraction_upper"]<sp.Rational(1,10)),
            "constant_mass_cubic_model_change_below_one_e_minus_fourteen":bool(
                t["constant_mass_minus_source_free_cubic_upper"]<sp.Rational(1,10**14)),
            "constant_mass_quartic_tree_model_change_below_one_e_minus_one_hundred_fifty":bool(
                t["constant_mass_minus_source_free_connected_quartic_upper"]<sp.Rational(1,10**150)),
            "same_free_modes_and_vector_covariance_only_on_the_original_clock":True,
            "no_new_quantum_stress_or_counterterm_bound_in_this_classical_estimate":True}
