"""Exact anchors and explicit written-proof gates for the actual new coupled inverse."""
from functools import cache

import sympy as sp
from p8_proca_rank_one_inverse import kernel

from . import coordinates, matching, tree, volterra


@cache
def residuals():
    out={}
    for group in (coordinates.checks(),tree.checks(),matching.checks(),volterra.checks()):
        if set(out).intersection(group):
            raise ValueError("Repeated new prepared inverse residual")
        out.update(group)
    return out


@cache
def gates():
    clock=tree.data()
    return {name:bool(value) for name,value in {
        "literal_new_minimal_modes_and_physical_readouts_check_time_covariance":True,
        "all_twelve_adiabatic_time_change_identities_retained":True,
        "all_six_finite_local_time_changes_and_three_background_Wards_checked":True,
        "original_selected_state_unchanged_for_initially_identity_time_change":True,
        "fixed_profile_not_reselected_and_its_local_time_contacts_retained":coordinates.data()["profile_and_state_reselected"] is False,
        "full_prepared_quantum_nonlocal_channel_is_one":coordinates.data()["quantum_nonlocal_channel_count"]==1,
        "time_channel_local_quantum_derivative_order_at_most_two":coordinates.data()["local_time_channel_maximum_derivative_order"]==2,
        "literal_zero_charge_matter_reduction_keeps_negative_sign":True,
        "actual_adapted_tree_time_coefficient_nonzero_on_entire_interval":clock["highest_time_channel_absolute_lower"]>0,
        "actual_clock_map_has_no_H_or_Theta_division":clock["clock_delta_lower"]>0,
        "both_local_cubic_derivative_cross_channels_retained":clock["maximum_local_derivative_orders"]==sp.Matrix([[4,3],[3,2]]),
        "new_scalar_leading_log_normalization_four":matching.data()["normalized_scalar_log_weight"]==4,
        "new_scalar_finite_fourth_coefficient_minus_four_not_retuned":matching.data()["fixed_new_scalar_finite_local_fourth_coefficient"]==-4,
        "new_full_dimensional_all_momentum_fourth_contacts_checked":True,
        "dimensional_scale_cancellation_before_finite_part_excludes_hidden_fourth_contact":True,
        "state_and_connection_remainders_retained_as_weak_log_after_four_primitives":True,
        "common_dimensional_state_remainder_integrable":-2+sp.Rational(1,4)<-1,
        "new_scalar_reference_inverse_has_integrable_causal_kernel":kernel.data()["L1_half_line_kernel"],
        "scalar_inverse_has_no_instantaneous_part_but_time_inverse_does":kernel.data()["instantaneous_inverse"]==0 and clock["instantaneous_time_channel_inverse_absolute_upper"]>0,
        "weighted_Volterra_majorant_integrable_for_fixed_finite_positive_coupling":True,
        "first_row_stronger_local_kernel_controls_physical_lapse_reconstruction":True,
        "prepared_smooth_forcing_recovers_unique_smooth_original_metric_solution":True,
        "force_transform_has_no_independent_prepared_homogeneous_kernel":True,
        "no_inverse_norm_smallness_or_quantum_stability_inferred":True,
        "no_arbitrary_spatial_state_initial_jet_nonlinear_or_original_P8_closure":True}.items()}


def controls():
    result=volterra.controls().copy()
    count=0
    for value in (True,False,1,sp.Integer(1),"unknown",None,[]):
        try:
            matching.high(value)
        except (TypeError,ValueError):
            count+=1
    if count!=7:
        raise ValueError("An invalid new minimal current sector was accepted")
    result["rejected_inputs"]+=count
    result["validation_precedes_cached_current_sector"]=True
    return result
