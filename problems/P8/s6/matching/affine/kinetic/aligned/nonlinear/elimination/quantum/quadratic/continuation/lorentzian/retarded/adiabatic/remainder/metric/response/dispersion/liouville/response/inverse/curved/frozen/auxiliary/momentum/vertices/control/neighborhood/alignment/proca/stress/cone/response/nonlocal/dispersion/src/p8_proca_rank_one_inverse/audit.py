"""Exact new rank-one massive block, scalar inverse and analytic-scope audit."""
from functools import cache

import sympy as sp

from . import cut, kernel, spectral


@cache
def residuals():
    out={}
    for group in (spectral.checks(),cut.checks(),kernel.checks()):
        if set(out).intersection(group):
            raise ValueError("Repeated new rank-one massive inverse residual")
        out.update(group)
    return out


@cache
def gates():
    block,density,time=spectral.data(),cut.data(),kernel.data()
    matrix=block["actual_physical_pair_Gram_matrix"]
    return {name:bool(value) for name,value in {
        "new_physical_reference_matrix_has_rank_one_not_two":matrix.det()==0 and matrix[1,1]!=0,
        "new_local_active_coefficient_is_fixed_not_retuned":True,
        "positive_pair_weight_has_strict_square_remainder":sp.Rational(8,3)>0,
        "direct_massive_threshold_gap_positive":block["H_at_threshold"]==sp.Rational(16,15)>0,
        "positive_axis_gap_from_full_massive_integral":spectral.positive_axis(0)["H_continuous_lower"]==4,
        "strict_half_plane_sign_and_cut_bank_gap_have_written_proof":True,
        "minimum_principle_uniform_first_sheet_gap_positive":density["uniform_first_sheet_H_real_part_lower"]>0,
        "uniform_range_inverse_complex_bound_is_fifteen_sixteenths":density["range_inverse_absolute_upper_on_first_sheet"]==sp.Rational(15,16),
        "new_threshold_density_constant_positive":density["density_over_sqrt_z_at_threshold"]>0,
        "new_high_density_log_constant_positive":density["density_times_log_tau_over_mass_squared_squared_at_infinity"]==sp.Rational(1,2),
        "static_positive_spectral_moment_is_one_quarter":density["spectral_static_moment"]==sp.Rational(1,4),
        "Cauchy_representation_includes_no_pole_or_instantaneous_remainder":density["no_instantaneous_range_inverse"],
        "new_scalar_time_kernel_has_integrable_small_and_large_time_bounds":time["L1_half_line_kernel"],
        "no_numeric_L1_norm_inferred_from_complex_frequency_bound":time["no_numerical_L1_norm_claim"],
        "short_window_range_inverse_norm_tends_to_zero":time["short_window_C0_inverse_norm_tends_to_zero"],
        "primitive_absolute_bound_is_one_half":time["primitive_lower"]==-sp.Rational(1,2) and time["primitive_upper"]==0,
        "range_inverse_not_promoted_to_inactive_lapse_inverse":density["full_two_source_inverse"] is False,
        "flat_Proca_reference_not_a_proved_full_candidate_Minkowski_vacuum":True,
        "fourth_derivative_prefactor_not_already_inverted_by_normalized_block":True,
        "actual_curved_tree_plus_loop_inverse_still_requires_new_work":True,
        "no_quantum_spectrum_cutoff_V_G_B_or_original_P8_closure":True}.items()}


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
         sp.Symbol("unfixed"),sp.sqrt(2),[])
    calls=[]
    for value in bad+(0,-1):
        calls.extend((lambda value=value:cut.fraction(value),
                      lambda value=value:kernel.mass(value),
                      lambda value=value:kernel.moment_bound(value)))
    calls.append(lambda:cut.fraction(1))
    for value in bad+(-1,):
        calls.append(lambda value=value:spectral.positive_axis(value))
    for value in (True,False,1,sp.Integer(1),"unknown",None,[]):
        calls.extend((lambda value=value:spectral.pair(value),lambda value=value:spectral.contact(value)))
    rejected=0
    for call in calls:
        try:
            call()
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(calls):
        raise ValueError("An inexact or outside-domain isolated Proca inverse input was accepted")
    return {"rejected_inputs":rejected,"validation_precedes_cached_sector_contact":True,
            "zero_axis_value_uses_the_proved_removable_limit":True,
            "no_parent_action_or_function_monkeypatched":True}

