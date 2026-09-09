"""Sharper all-frequency transport and exact high-frequency normal-form constants."""
from functools import cache

import sympy as sp
from p8_proca_nearby_packets import domains, finite, normal_form
from p8_proca_nearby_retarded import remainder as old_remainder

from . import basis

B0_UPPER=sp.Integer(250000)
B0_DOT_UPPER=sp.Integer(5)*10**12
Z_UPPER=sp.Integer(3)*10**12
Z_DOT_UPPER=sp.Integer(6)*10**19
W_UPPER=sp.Integer(10)**6
REMAINDER_UPPER=sp.Integer(2)*10**20
MIN_MOMENTUM=sp.Integer(10)**14
SCALAR_ERROR=sp.Integer(10)**16
SCALAR_MULTIPLIER=sp.Integer(100)


def momentum(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require an exact rational high-frequency momentum")
    value=sp.Rational(value)
    if value<MIN_MOMENTUM:
        raise ValueError("Momentum is below the new proved high-frequency domain")
    return value


@cache
def data():
    S,Si,Sd=6,4,sp.Rational(1,50)
    L=domains.MATRIX_NORM
    b0=Si*L*S+Si*Sd
    w=3*Si*L*S
    perturb=b0+w
    delta=domains.T*10**6
    low=old_remainder.low_frequency()["complete_original_generator_coefficient_majorants"]
    lowmat=sum((coefficient*4**j for j,coefficient in enumerate(low)),sp.zeros(4))
    lowB=max(sum(lowmat[i,j] for j in range(4)) for i in range(4))
    zd=B0_DOT_UPPER/domains.GAP+2*sp.Rational(1,100)*B0_UPPER/domains.GAP**2
    rem=2*(B0_UPPER*Z_UPPER+sp.Rational(3,2)*W_UPPER+Z_DOT_UPPER)
    matrix_error=S*Si*(6*Z_UPPER+4*domains.T*REMAINDER_UPPER)
    return {"actual_inner_interval_length":domains.T,
        "actual_complex_disc_and_inner_Cauchy_radius":(domains.T,domains.T/2),
        "actual_B0_triangle_upper":b0,"declared_B0_upper":B0_UPPER,
        "Cauchy_B0_first_derivative_upper":B0_DOT_UPPER,
        "actual_all_inverse_frequency_triangle_upper":w,"declared_W_upper":W_UPPER,
        "all_k_at_least_one_interaction_picture_generator_upper":perturb,
        "all_k_at_least_one_integrated_generator_upper":delta,
        "exact_interaction_picture_transfer_norm_upper":1/(1-delta),
        "declared_all_k_at_least_one_packet_transfer_norm_upper":sp.Integer(48),
        "declared_all_k_at_least_one_scalar_response_coefficient":SCALAR_MULTIPLIER,
        "unit_low_frequency_ball_original_generator_entry_majorants":lowmat,
        "unit_low_frequency_ball_original_generator_norm_upper":lowB,
        "unit_low_frequency_ball_integrated_generator_upper":domains.T*lowB,
        "declared_unit_low_frequency_ball_scalar_response_upper":sp.Integer(4),
        "declared_Z_upper":Z_UPPER,"actual_Z_dot_triangle_upper":zd,
        "declared_Z_dot_upper":Z_DOT_UPPER,
        "actual_normal_form_remainder_triangle_upper":rem,
        "declared_normal_form_remainder_upper":REMAINDER_UPPER,
        "new_minimum_high_frequency_momentum":MIN_MOMENTUM,
        "actual_two_endpoint_packet_error_coefficient":matrix_error,
        "actual_two_endpoint_scalar_error_coefficient":2*matrix_error,
        "declared_two_endpoint_scalar_error_coefficient":SCALAR_ERROR,
        "exact_normal_form_change_matrix":normal_form.data()["Z"],
        "all_frequencies_covered_without_a_finite_q_Legendre_inverse":True}


def error(value=MIN_MOMENTUM):
    value=momentum(value)
    return {"momentum":value,"near_identity_difference_upper":Z_UPPER/value,
        "integrated_normal_form_remainder_upper":domains.T*REMAINDER_UPPER/value,
        "uniform_two_endpoint_scalar_response_error_upper":SCALAR_ERROR/value**2}


@cache
def coefficient_majorants():
    p=finite.parent
    limits={p.a:(0,3),p.m:(0,8),p.g:(0,2),p.r:(0,8),
        p.rn:(sp.Rational(1,10),1),p.h:(1,6),p.ell:(0,sp.Rational(1,8)),
        p.alpha:(0,1),p.beta:(0,1),p.H:(0,1),finite.R:(sp.Rational(1,2),2)}
    return {name:finite.data()[name].applyfunc(lambda v:domains.laurent_majorant(v,limits))
        for name in ("L0","L1","L2","L3")}


@cache
def checks():
    d=data()
    return {"direct_small_basis_B0_triangle_has_both_terms":
        d["actual_B0_triangle_upper"]-(4*10000*6+4*sp.Rational(1,50)),
        "all_three_inverse_frequency_terms_kept_for_k_at_least_one":
        d["actual_all_inverse_frequency_triangle_upper"]-3*4*10000*6,
        "Cauchy_B0_derivative_uses_actual_inner_radius":
        B0_DOT_UPPER-B0_UPPER/(domains.T/2),
        "unit_low_frequency_original_generator_keeps_all_q_orders":
        d["unit_low_frequency_ball_original_generator_norm_upper"]-sp.Rational(28169,64),
        "gap_derivative_contributes_to_new_change_derivative":
        d["actual_Z_dot_triangle_upper"]-B0_DOT_UPPER/domains.GAP
        -2*sp.Rational(1,100)*B0_UPPER/domains.GAP**2,
        "two_time_error_keeps_both_basis_endpoints_and_both_small_parameters":
        d["actual_two_endpoint_packet_error_coefficient"]-24*(6*Z_UPPER+4*domains.T*REMAINDER_UPPER),
        "new_scalar_error_retains_physical_source_and_output":
        d["actual_two_endpoint_scalar_error_coefficient"]-2*d["actual_two_endpoint_packet_error_coefficient"],
        "new_normal_form_smallness_is_one_fifth_at_lower_endpoint":
        error()["integrated_normal_form_remainder_upper"]-sp.Rational(1,5)}


@cache
def gates():
    d,e=data(),error()
    out={"all_exact_Laurent_L_matrix_norms_below_inherited_majorant":
        all(max(sum(matrix[i,j] for j in range(4)) for i in range(4))<10000
            for matrix in coefficient_majorants().values()),
        "actual_B0_below_declared":d["actual_B0_triangle_upper"]<B0_UPPER,
        "actual_W_below_declared":d["actual_all_inverse_frequency_triangle_upper"]<W_UPPER,
        "complete_interaction_picture_generator_below_one_million":
            d["all_k_at_least_one_interaction_picture_generator_upper"]<10**6,
        "integrated_interaction_picture_generator_is_below_one_half":
            d["all_k_at_least_one_integrated_generator_upper"]<sp.Rational(1,2),
        "interaction_picture_transfer_below_two":d["exact_interaction_picture_transfer_norm_upper"]<2,
        "scalar_high_response_bound_keeps_two_endpoint_prefactor":2*48<SCALAR_MULTIPLIER,
        "unit_ball_original_generator_integral_below_one_half":
            d["unit_low_frequency_ball_integrated_generator_upper"]<sp.Rational(1,2),
        "unit_ball_source_factor_below_two":sp.Rational(101,100)**4<2,
        "actual_Z_below_declared":B0_UPPER/domains.GAP<Z_UPPER,
        "actual_Z_dot_below_declared":d["actual_Z_dot_triangle_upper"]<Z_DOT_UPPER,
        "actual_normal_form_remainder_below_declared":
            d["actual_normal_form_remainder_triangle_upper"]<REMAINDER_UPPER,
        "new_change_invertible_on_complete_declared_band":e["near_identity_difference_upper"]<sp.Rational(1,2),
        "new_integrated_remainder_at_most_one_quarter":e["integrated_normal_form_remainder_upper"]<=sp.Rational(1,4),
        "new_two_time_scalar_error_below_declared":
            d["actual_two_endpoint_scalar_error_coefficient"]<SCALAR_ERROR,
        "new_high_frequency_threshold_is_below_old_one":MIN_MOMENTUM<normal_form.MIN_MOMENTUM,
        "actual_complex_basis_bounds_not_only_real_point_bounds":all(basis.gates().values())}
    return {name:bool(value) for name,value in out.items()}


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
         sp.Symbol("unbounded"),sp.sqrt(2),[],0,-1,MIN_MOMENTUM-1)
    count=0
    for value in bad:
        for call in (momentum,error):
            try:
                call(value)
            except (TypeError,ValueError):
                count+=1
            else:
                raise ValueError("An invalid sharpened frequency domain was accepted")
    return {"rejected_inputs":count}
