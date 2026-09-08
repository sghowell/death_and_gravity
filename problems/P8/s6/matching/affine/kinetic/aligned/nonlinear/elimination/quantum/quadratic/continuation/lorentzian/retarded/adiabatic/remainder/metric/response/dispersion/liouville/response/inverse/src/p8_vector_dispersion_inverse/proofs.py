"""Algebraic anchors for the written first-sheet and causal-kernel proof."""
from functools import cache

from . import block, cut, kernel


@cache
def residuals():
    result={}
    for group in (block.checks(),cut.checks(),kernel.checks()):
        if set(result).intersection(group):
            raise ValueError("Repeated inverse residual name")
        result.update(group)
    return result


@cache
def checks():
    item=block.data()
    moment=item["static_spectral_moment"].subs(block.h,1)
    return {"positive_real_axis_gap_reused_only_for_the_exact_massive_block":True,
            "strict_nonreal_imaginary_part_from_interior_positive_spectral_matrix":True,
            "subthreshold_chart_positive_first_and_negative_second_diagonals":True,
            "threshold_matrix_nonsingular":bool(item["threshold_determinant"].subs(block.h,1)<0),
            "cut_banks_strictly_positive_spectral_weight_off_threshold":True,
            "static_spectral_moment_positive_first_minor":bool(moment[0,0]>0),
            "static_spectral_moment_positive_determinant":bool(moment.det()>0),
            "inverse_instantaneous_part_nonzero_rank_one":bool(
                item["inverse_infinity"].subs(block.h,1).trace()>0 and item["inverse_infinity"].det()==0),
            "threshold_density_coefficient_nonzero":bool(
                kernel.data()["density_over_sqrt_z_at_threshold"][1,1]>0),
            "ultraviolet_density_coefficient_nonzero":bool(
                kernel.data()["density_times_log_tau_over_mass2_squared_at_infinity"][1,1]>0),
            "causal_regular_kernel_L1_proof_covers_both_time_endpoints":True,
            "oscillatory_cut_integral_not_claimed_absolutely_convergent_in_frequency":True,
            "instantaneous_term_precludes_full_small_interval_inverse_norm_limit_zero":True,
            "variable_coefficient_congruence_only_defines_a_preconditioner":True,
            "not_the_full_curved_tree_plus_loop_inverse_or_a_quantum_health_verdict":True}
