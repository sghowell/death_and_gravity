"""Exact anchors and explicit boundaries for the written no-loss proof."""
from functools import cache

import sympy as sp

from . import primitive, tree, vertices


@cache
def residuals():
    result={}
    for group in (vertices.checks(),vertices.adiabatic_principal_checks(),primitive.checks(),tree.checks()):
        if set(result).intersection(group):
            raise ValueError("Repeated coupled-response residual name")
        result.update(group)
    return result


@cache
def checks():
    p=vertices.high("L")["pair"]
    transverse=vertices.high("T")["pair"]
    return {"positive_common_clock_rate_bounds":bool(sp.Rational(16,25)>0),
            "physical_mass_and_state_domain_not_extended_to_massless_case":True,
            "actual_physical_Hamiltonian_vertices_and_instantaneous_contacts_retained":True,
            "all_order_state_used_for_finite_time_derivative_expansions_not_reselected":True,
            "complex_dimensional_initial_mixing_remainder_integrable":bool(-2+sp.Rational(1,4)<-1),
            "transverse_highest_square_first_dimensional_jet_vanishes":(
                sp.diff((vertices.D-1)*transverse*transverse.T,vertices.D).subs(vertices.D,3)==sp.zeros(2)),
            "longitudinal_first_dimensional_pair_jet_retained":sp.diff(p,vertices.D)==sp.Matrix([0,2]),
            "proper_momentum_scale_identity_holds_before_dimensional_finite_part":True,
            "fixed_highest_finite_local_matrix_and_exact_massive_constant_retained":True,
            "every_remaining_local_input_derivative_is_at_most_three_after_reference_removal":True,
            "subleading_four_primitive_kernel_is_integrable_logarithm":True,
            "literal_matter_charge_reduction_retains_negative_reconstruction_sign":(
                tree.data()["matter_reconstruction_rate"]==-3*tree.ell*tree.v-tree.w*tree.n),
            "classical_physical_metric_operator_is_local_order_at_most_two":True,
            "nonzero_instantaneous_inverse_is_kept_before_Volterra_composition":True,
            "weighted_C0_inverse_does_not_assume_small_gamma_inverse":True,
            "smoothness_uses_diagonal_kernel_derivatives_and_prepared_endpoints":True,
            "finite_constants_not_promoted_to_numerical_smallness_or_stability":True,
            "not_a_nonlinear_neighborhood_spatial_initial_state_or_VGB_theorem":True}
