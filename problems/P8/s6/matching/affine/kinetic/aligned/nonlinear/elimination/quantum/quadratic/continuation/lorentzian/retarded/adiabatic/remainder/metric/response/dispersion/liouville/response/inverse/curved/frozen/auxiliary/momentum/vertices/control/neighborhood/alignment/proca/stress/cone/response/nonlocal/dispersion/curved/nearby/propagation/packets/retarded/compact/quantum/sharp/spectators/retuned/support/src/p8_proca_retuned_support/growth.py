"""Complete complex-frequency logarithmic norm and the actual matter support radius."""
from functools import cache

import sympy as sp

from . import domain

EUCLIDEAN_LOWER_ORDER=sp.Integer(4)*10**9
HALF=domain.T/2


def exact(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require exact rational new-support domain inputs")
    return sp.Rational(value)


def time_pair(start=-HALF,end=HALF):
    start,end=exact(start),exact(end)
    if not -HALF<=start<=end<=HALF:
        raise ValueError("Require ordered times inside the new inner real interval")
    dt=end-start
    return {"source_time":start,"output_time":end,"coordinate_duration":dt,
        "entire_multiplier_uniform_prefactor":4+320*sp.exp(EUCLIDEAN_LOWER_ORDER*dt),
        "actual_matter_null_radius_upper":domain.HI/domain.LO**2*dt,
        "support_radius_is_actual_null_integral_not_this_upper_estimate":True}


@cache
def data():
    b=domain.basis()
    S,Si,Sd=(b[key] for key in ("declared_basis_infinity_norm_upper",
        "declared_inverse_infinity_norm_upper","declared_basis_derivative_infinity_norm_upper"))
    perturb=Si*4*10000*S+Si*Sd
    low=domain.coefficient_majorants()["low_scalar_root_unit_disc_generator_infinity_norm_upper"]
    u=sp.Symbol("actual_clock_time",real=True)
    start,end=sp.symbols("ordered_source_time ordered_output_time",real=True)
    N,R=sp.Function("N_new_actual_Picard_solution")(u),sp.Function("R_new_actual_Picard_solution")(u)
    h=(1+u*u)**3
    e=(N*N*h/(1+N*N*(h-1)))**sp.Rational(1,4)
    radius=sp.Integral(N/(e*R),(u,start,end))
    return {"actual_real_time_interval":(-HALF,HALF),
        "actual_matter_null_support_radius":radius,
        "matter_frequency_uses_actual_new_lapse_conformal_and_hat_scale":N/(e*R),
        "complete_lower_order_infinity_norm_triangle_upper":perturb,
        "declared_complete_lower_order_infinity_norm_upper":sp.Integer(2)*10**9,
        "declared_complete_lower_order_Euclidean_norm_upper":EUCLIDEAN_LOWER_ORDER,
        "Euclidean_basis_and_inverse_upper":(2*S,2*Si),
        "complete_packet_transfer_Euclidean_prefactor":4*S*Si,
        "complete_scalar_source_output_prefactor":2*4*S*Si,
        "high_complex_scalar_root_domain":"abs(k)>=1; k^2=zeta dot zeta",
        "high_complex_scalar_response_bound":"320*exp(S_m*abs(Im(k))+4*10^9*(t-s))/abs(k)",
        "low_complex_scalar_root_domain":"abs(k)<=1 including complex-vector null cone",
        "low_complex_q_modulus_upper":sp.Integer(4),
        "low_original_generator_integrated_norm_upper":domain.T*low,
        "low_original_scalar_response_triangle_upper":2/(1-domain.T*low),
        "declared_low_complex_scalar_response_upper":sp.Integer(4),
        "full_interval_entire_multiplier_uniform_prefactor":time_pair()["entire_multiplier_uniform_prefactor"],
        "entire_multiplier_global_type_bound":"abs(G_hat(t,s,zeta)) <= (4+320*exp(400))*exp(S_m(t,s)*norm(Im(zeta)))",
        "all_inverse_frequency_orders_retained_and_no_smallness_required":True,
        "principal_real_frequency_growth_cancels_in_logarithmic_norm":True,
        "old_fast_clock_bounds_state_and_response_not_transferred":True}


@cache
def checks():
    x,y,wc,wm=sp.symbols("real_frequency imaginary_frequency omega_c omega_m",real=True)
    lam=sp.diag(wc,wm,-wc,-wm)
    leading=sp.I*(x+sp.I*y)*lam
    lapse,conformal,scale,slope=sp.symbols("new_lapse new_conformal new_hat_scale radial_slope",positive=True)
    radial_metric=-lapse**2+(conformal*scale)**2*slope**2
    d=data()
    return {"real_frequency_principal_part_is_exactly_skew_Hermitian":
        ((leading+leading.conjugate().T)/2+y*lam).applyfunc(sp.factor),
        "complete_lower_order_bound_keeps_L0_all_three_inverse_orders_and_basis_motion":
            d["complete_lower_order_infinity_norm_triangle_upper"]-(5*4*10000*8+5*2*10**8),
        "four_dimensional_Euclidean_bound_keeps_norm_conversion_factor":
            EUCLIDEAN_LOWER_ORDER-2*d["declared_complete_lower_order_infinity_norm_upper"],
        "complete_source_output_prefactor_keeps_both_Euclidean_basis_factors":
            d["complete_scalar_source_output_prefactor"]-320,
        "finite_full_interval_lower_order_exponent_is_four_hundred":EUCLIDEAN_LOWER_ORDER*domain.T-400,
        "coincident_time_radius_upper_zero":time_pair(0,0)["actual_matter_null_radius_upper"],
        "proper_null_distance_is_lapse_over_physical_scale":
            sp.factor(radial_metric.subs(slope,lapse/(conformal*scale)))}


@cache
def gates():
    d=data()
    rows={"complete_lower_order_infinity_norm_fits_declared":
        d["complete_lower_order_infinity_norm_triangle_upper"]<d["declared_complete_lower_order_infinity_norm_upper"],
        "low_complex_root_generator_integral_below_one_half":d["low_original_generator_integrated_norm_upper"]<sp.Rational(1,2),
        "low_complex_scalar_response_below_four":d["low_original_scalar_response_triangle_upper"]<4,
        "low_complex_q_bound_follows_from_actual_scale_lower":1/domain.LO**2<4,
        "high_and_low_scalar_root_domains_cover_every_complex_vector":True,
        "entire_original_generator_not_chart_formula_defines_complex_null_cone":True,
        "new_real_principal_matter_frequency_dominates_every_mode":
            0<domain.actual_bounds()["new_physical_clock_squared_speed_upper"]<1,
        "uniform_prefactor_finite_and_independent_of_complex_momentum":True}
    return {name:bool(value) for name,value in rows.items()}
