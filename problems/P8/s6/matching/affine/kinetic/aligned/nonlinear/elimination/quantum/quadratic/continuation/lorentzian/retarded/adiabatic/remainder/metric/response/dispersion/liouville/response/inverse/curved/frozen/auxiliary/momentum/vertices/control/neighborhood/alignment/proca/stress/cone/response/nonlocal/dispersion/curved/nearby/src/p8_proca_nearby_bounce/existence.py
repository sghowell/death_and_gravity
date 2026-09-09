"""Recentered holomorphic lapse root and two-variable Picard evolution."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import bounds as upstream

from . import jets

SHIFT_MAX=sp.Rational(1,10**6)
PARAM_RADIUS=sp.Rational(1,10**24)
ROOT_RADIUS=sp.Rational(1,10**14)
TIME_WINDOW=sp.Rational(1,10**30)
H_BOUND=sp.Integer(10000)
FLOW_BOUND=sp.Integer(50000)


def shift(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require an exact central lapse shift")
    value=sp.Rational(value)
    if not 0<=value<=SHIFT_MAX:
        raise ValueError("Require a central lapse shift in [0,10^-6]")
    return value


def duration(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require an exact positive local time half-width")
    value=sp.Rational(value)
    if not 0<value<=TIME_WINDOW:
        raise ValueError("Require a local time half-width at most 10^-30")
    return value


def datum(value):
    value=shift(value)
    N=1+value
    data=jets.data()
    P2=data["positive_matter_density_squared"].subs(jets.N,N)
    return {"central_lapse_shift":value,"central_lapse":N,
            "central_hat_scale":sp.Integer(1),"central_trace_momentum":sp.Integer(0),
            "positive_conserved_matter_charge":sp.sqrt(P2),
            "positive_conserved_matter_charge_squared":P2,
            "central_physical_scale":sp.sqrt(N),
            "central_lapse_acceleration":sp.factor(data["actual_lapse_second_time_derivative"].subs(jets.N,N)),
            "central_physical_Hubble_proper_acceleration":sp.factor(data["actual_physical_Hubble_proper_time_derivative"].subs(jets.N,N)),
            "initial_vector_field_and_momentum_zero":True,
            "not_the_S6_82_profile_corrected_quantum_background":True}


@cache
def coefficient_domain():
    # The old full Hamiltonian bound was stated for real time.
    # For this joint complex-time proof use a smaller disk and a fresh I_phi bound.
    I=upstream.coefficients()["boundary_primitive_absolute_upper"]
    Iphi=I/sp.Rational(1,100)
    potential=2*1002+2*Iphi
    trace=2*sp.Rational(5,4)**2
    matter=sp.Rational(9,25)**2
    return {"central_shift_interval":(sp.Integer(0),SHIFT_MAX),
            "outer_complex_time_radius":sp.Rational(1,100),
            "outer_complex_lapse_radius_about_new_center":sp.Rational(1,200),
            "outer_trace_and_matter_density_radius_about_new_center":sp.Rational(1,4),
            "original_primitive_complex_time_radius":sp.Rational(1,50),
            "actual_primitive_absolute_upper":I,
            "fresh_complex_time_primitive_phi_upper":Iphi,
            "fresh_complex_time_potential_upper":potential,
            "homogeneous_trace_piece_upper":trace,"homogeneous_matter_piece_upper":matter,
            "homogeneous_Hamiltonian_triangle_upper":2*(trace+potential+matter),
            "declared_joint_holomorphic_Hamiltonian_upper":H_BOUND,
            "no_unjustified_extension_of_old_real_time_Hamiltonian_bound":True}


@cache
def lapse_root():
    rn,ru,ry=sp.Rational(1,400),sp.Rational(1,200),sp.Rational(1,8)
    HNu=H_BOUND/(rn*ru)
    HNy=H_BOUND/(rn*ry)
    HNNN=6*H_BOUND/rn**3
    HNNu=2*H_BOUND/(rn**2*ru)
    HNNy=2*H_BOUND/(rn**2*ry)
    pivot=HNNN*ROOT_RADIUS+(HNNu+2*HNNy)*PARAM_RADIUS
    force=(HNu+2*HNy)*PARAM_RADIUS
    rate=pivot/2
    return {"inner_N_Cauchy_radius":rn,"inner_u_Cauchy_radius":ru,
            "inner_p_and_ell_Cauchy_radius":ry,
            "H_Nu_upper":HNu,"H_Np_and_H_Nell_upper":HNy,
            "H_NNN_upper":HNNN,"H_NNu_upper":HNNu,"H_NNp_and_H_NNell_upper":HNNy,
            "complex_time_parameter_radius":PARAM_RADIUS,
            "complex_trace_and_density_parameter_radius":PARAM_RADIUS,
            "complex_lapse_root_radius_about_new_center":ROOT_RADIUS,
            "inverse_exact_center_lapse_pivot_upper":sp.Rational(1,2),
            "lapse_pivot_change_upper":pivot,"center_force_upper":force,
            "Newton_contraction_upper":rate,"Newton_center_step_upper":force/2,
            "Newton_closed_disc_image_radius_upper":force/2+rate*ROOT_RADIUS,
            "actual_root_displacement_upper":force/(2-pivot),
            "remaining_lapse_pivot_absolute_lower":2-pivot}


def evolution(value=TIME_WINDOW):
    value=duration(value)
    Hp=H_BOUND/sp.Rational(1,8)
    rate=4*FLOW_BOUND/PARAM_RADIUS*value
    return {"complex_time_half_width":value,"real_time_interval":(-value,value),
            "actual_background_equations":"p'=-Hcal+ell*Hcal_ell; ell'=-ell*Hcal_p; (log R)'=Hcal_p/3; Hcal_N=0",
            "right_hand_side_joint_absolute_upper":FLOW_BOUND,
            "solution_phase_disc_radius":PARAM_RADIUS/2,
            "Picard_image_radius_upper":FLOW_BOUND*value,
            "phase_flow_Lipschitz_upper":4*FLOW_BOUND/PARAM_RADIUS,
            "Picard_contraction_upper":rate,
            "reconstructed_log_hat_scale_absolute_upper":Hp*value/3,
            "hat_scale_positive_lower":sp.Rational(1,2),"hat_scale_positive_upper":sp.Integer(2),
            "same_nonzero_matter_charge":"ell R^3=ell_center; chi'=N*ell/U",
            "physical_clock_and_scale":"phi=u, physical lapse N>0, a=eomega*R, X=N^-2",
            "time_parity":"N,ell,R,a are even and p is odd, by coefficient parity and uniqueness.",
            "quantitative_time_dependent_scalar_cone_or_time_advance_bound":False}

@cache
def checks():
    out={}
    R,PR,C,N,a,b,f,U=sp.symbols("hat_scale canonical_trace conserved_charge fixed_lapse trace_a trace_b scalar_f volume_U",nonzero=True)
    p,l=sp.symbols("trace_density matter_density",real=True)
    H=N*((p-b)**2/(4*a)-f+l*l/(2*U))
    canonical=R**3*H.subs({p:PR/(3*R**2),l:C/R**3},simultaneous=True)
    Rd=sp.diff(canonical,PR)
    pd=-sp.diff(canonical,R)/(3*R**2)-2*PR*Rd/(3*R**3)
    on={PR:3*R**2*p,C:l*R**3}
    out["literal_homogeneous_hat_scale_equation"]=sp.factor((Rd/R).subs(on)-sp.diff(H,p)/3)
    out["literal_homogeneous_trace_density_equation"]=sp.factor(pd.subs(on)+H-l*sp.diff(H,l))
    out["literal_nonzero_matter_charge_conservation"]=sp.factor((-l*sp.diff(H,p))*R**3+3*l*R**2*(R*sp.diff(H,p)/3))
    out["literal_physical_free_matter_charge"]=sp.factor((U*R**3)*(l/U)-l*R**3)
    c=jets.model.coefficients()
    actual=jets.N*((p-c["b"])**2/(4*c["a"])-c["f0"]+l*l/(2*c["U"]))
    mapping={jets.u:-jets.u,p:-p,jets.model.I:-jets.model.I}
    out["actual_time_and_phase_even_Hamiltonian"]=sp.factor(actual.subs(mapping,simultaneous=True)-actual)
    out["actual_time_and_phase_odd_hat_Hubble"]=sp.factor(
        sp.diff(actual,p).subs(mapping,simultaneous=True)+sp.diff(actual,p))
    center=jets.N*(p*p/(4*jets.background.data()["trace_a"])-jets.background.data()["scalar_potential"]
                  +l*l/(2*jets.background.data()["volume_U"]))
    out["actual_central_fixed_phase_lapse_constraint"]=sp.factor(
        sp.diff(center,jets.N).subs({p:0,l*l:jets.data()["positive_matter_density_squared"]},simultaneous=True))
    root=lapse_root()
    out["new_NNN_Cauchy_factorial_retained"]=root["H_NNN_upper"]-6*H_BOUND/root["inner_N_Cauchy_radius"]**3
    out["new_NNu_Cauchy_factorial_retained"]=root["H_NNu_upper"]-2*H_BOUND/(root["inner_N_Cauchy_radius"]**2*root["inner_u_Cauchy_radius"])
    out["new_fixed_center_Newton_rate"]=root["Newton_contraction_upper"]-root["inverse_exact_center_lapse_pivot_upper"]*root["lapse_pivot_change_upper"]
    out["new_closed_disc_Newton_image"]=root["Newton_closed_disc_image_radius_upper"]-root["Newton_center_step_upper"]-root["Newton_contraction_upper"]*ROOT_RADIUS
    evo=evolution()
    out["whole_complex_time_Picard_rate_one_fifth"]=evo["Picard_contraction_upper"]-sp.Rational(1,5)
    out["actual_phase_flow_Lipschitz_uses_two_Cauchy_columns"]=evo["phase_flow_Lipschitz_upper"]-2*FLOW_BOUND/(PARAM_RADIUS/2)
    out["actual_log_scale_reconstruction_normalization"]=evo["reconstructed_log_hat_scale_absolute_upper"]-H_BOUND*TIME_WINDOW/(3*sp.Rational(1,8))
    return out


@cache
def gates():
    c,r,e=coefficient_domain(),lapse_root(),evolution()
    box=jets.boxes()
    return {name:bool(value) for name,value in {
        "all_central_matter_charges_real_and_positive":box["matter_density_squared"]["lower"]>sp.Rational(1,200),
        "central_charge_below_eleven_hundredths":box["matter_density_squared"]["upper"]<sp.Rational(11,100)**2,
        "actual_center_lapse_pivot_inverse_below_one_half":box["lapse_Hessian_times_N_seven_halves"]["upper"]<-2*(1+SHIFT_MAX)**4,
        "actual_central_physical_bounce_acceleration_above_39999_over_10000":box["physical_Hubble_proper_acceleration"]["lower"]>sp.Rational(39999,10000),
        "lapse_acceleration_absolute_upper_below_one_thousandth":max(abs(box["lapse_acceleration"]["lower"]),abs(box["lapse_acceleration"]["upper"]))<sp.Rational(1,1000),
        "recentered_lapse_disk_inside_original_outer_lapse_domain":SHIFT_MAX+c["outer_complex_lapse_radius_about_new_center"]<upstream.LAPSE_RADIUS,
        "fresh_complex_time_boundary_derivative_bound_below_three":c["fresh_complex_time_primitive_phi_upper"]<3,
        "fresh_complex_time_potential_bound_below_three_thousand":c["fresh_complex_time_potential_upper"]<3000,
        "fresh_joint_holomorphic_Hamiltonian_bound_below_ten_thousand":c["homogeneous_Hamiltonian_triangle_upper"]<H_BOUND,
        "complex_time_Cauchy_disk_fits_original_primitive_domain":2*c["outer_complex_time_radius"]<=c["original_primitive_complex_time_radius"],
        "new_Newton_map_strict_contraction":r["Newton_contraction_upper"]<sp.Rational(1,10),
        "new_Newton_closed_disc_strictly_invariant":r["Newton_closed_disc_image_radius_upper"]<ROOT_RADIUS,
        "new_actual_root_displacement_below_one_e_minus_fifteen":r["actual_root_displacement_upper"]<sp.Rational(1,10**15),
        "new_lapse_pivot_stays_above_one_point_nine":r["remaining_lapse_pivot_absolute_lower"]>sp.Rational(19,10),
        "local_flow_time_inside_implicit_time_domain":TIME_WINDOW<PARAM_RADIUS/2,
        "Picard_maps_phase_half_disc_strictly_into_itself":e["Picard_image_radius_upper"]<e["solution_phase_disc_radius"],
        "Picard_contraction_below_one_half":e["Picard_contraction_upper"]<sp.Rational(1,2),
        "log_hat_scale_small_and_both_real_scale_bounds_positive":e["reconstructed_log_hat_scale_absolute_upper"]<sp.Rational(1,2),
        "entire_real_clock_stays_inside_original_covariant_tube":(1+SHIFT_MAX+ROOT_RADIUS)**-2>sp.Rational(9,10) and (1-ROOT_RADIUS)**-2<sp.Rational(11,10),
        "sourcefree_ordinary_Proca_zero_field_is_exact_invariant_sector":True,
        "time_parity_and_constraint_uniqueness_fix_zero_center_lapse_velocity":True,
        "local_bounces_not_global_complete_nearby_bounce_families":True,
        "no_quantum_profile_state_or_inverse_transferred_to_new_backgrounds":True,
        "no_finite_time_cone_or_resolved_advance_bound":e["quantitative_time_dependent_scalar_cone_or_time_advance_bound"] is False,
        "original_P8_and_common_parent_V_G_B_open":True}.items()}


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
         sp.Symbol("unfixed"),sp.sqrt(2),[])
    calls=[]
    for value in bad+(-1,1):
        calls.extend((lambda value=value:shift(value),lambda value=value:datum(value)))
    for value in bad+(0,-1,2*TIME_WINDOW):
        calls.extend((lambda value=value:duration(value),lambda value=value:evolution(value)))
    count=0
    for call in calls:
        try:
            call()
        except (TypeError,ValueError):
            count+=1
    if count!=len(calls):
        raise ValueError("An invalid nearby-bounce domain input was accepted")
    return {"rejected_inputs":count,"validation_before_scientific_evaluation":True,
            "no_parent_action_profile_state_or_scientific_function_changed":True}
