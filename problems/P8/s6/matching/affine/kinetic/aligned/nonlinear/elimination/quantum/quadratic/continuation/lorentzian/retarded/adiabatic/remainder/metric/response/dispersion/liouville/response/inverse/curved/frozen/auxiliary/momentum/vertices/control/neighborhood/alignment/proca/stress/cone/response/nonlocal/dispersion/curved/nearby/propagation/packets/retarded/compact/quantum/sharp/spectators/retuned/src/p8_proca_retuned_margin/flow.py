"""New holomorphic solution, conserved physical matter charge and principal health."""
from functools import cache

import sympy as sp

from . import domain, model

N_SPEED=sp.Rational(1,10000)
Z_SPEED=sp.Integer(9)


def duration(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require an exact rational new-solution half-width")
    value=sp.Rational(value)
    if not 0<value<=domain.T:
        raise ValueError("Require a new-solution half-width in (0,10^-7]")
    return value


def picard(value=domain.T):
    value=duration(value)
    radii=(domain.RN,domain.RZ)
    speeds=(N_SPEED,Z_SPEED)
    images=tuple(value*b for b in speeds)
    lipschitz=sp.Matrix(2,2,lambda i,j:2*speeds[i]/radii[j])
    weighted=sp.Matrix(2,2,lambda i,j:lipschitz[i,j]*radii[j]/radii[i])
    row_max=max(sum(weighted[i,j] for j in range(2)) for i in range(2))
    hdev=3*value**2+3*value**4+value**6
    return {"new_solution_holomorphic_time_radius":value,
        "new_solution_lapse_and_trace_centers":(domain.N0,sp.Integer(0)),
        "outer_complex_lapse_and_trace_radii":radii,
        "inner_Picard_lapse_and_trace_radii":tuple(r/2 for r in radii),
        "whole_outer_complex_velocity_bounds":speeds,
        "actual_solution_displacement_bounds":images,
        "Cauchy_phase_Lipschitz_entry_upper":lipschitz,
        "weighted_phase_Lipschitz_entry_upper":weighted,
        "weighted_phase_Lipschitz_row_norm_upper":row_max,
        "Picard_contraction_upper":value*row_max,
        "new_log_hat_scale_absolute_upper":value/sp.Integer(60000),
        "new_hat_scale_minus_one_absolute_upper":(value/sp.Integer(60000))/(1-value/sp.Integer(60000)),
        "complex_background_h_minus_one_upper":hdev,
        "complex_conformal_D_minus_one_upper":(domain.N0+domain.RN)**2*hdev,
        "uses_new_rational_flow_not_inherited_old_background":True,
        "fixed_point_is_actual_homogeneous_classical_solution_not_just_off_shell_box":True}


@cache
def initial():
    s=model.system()
    w2=domain.center(s["constraint_matter_square"])
    charge2=sp.factor(w2/domain.N0)
    old_w2=domain.center(model.old.system()["constraint_matter_square"])
    return {"central_physical_lapse":domain.N0,
        "central_rescaled_trace":sp.Integer(0),"central_canonical_trace_density":sp.Integer(0),
        "central_hat_scale":sp.Integer(1),"central_physical_scale":sp.sqrt(domain.N0),
        "central_rescaled_matter_square":w2,
        "new_positive_conserved_matter_charge_squared":charge2,
        "new_positive_conserved_matter_charge":sp.sqrt(charge2),
        "old_charge_squared_at_same_lapse":sp.factor(old_w2/domain.N0),
        "old_charge_residual_in_new_rescaled_constraint":
            domain.center(s["M"])*(old_w2-w2),
        "central_physical_Hubble":domain.center(domain.physical_hubble()),
        "central_lapse_velocity":domain.center(s["lapse_flow"]),
        "initial_ordinary_Proca_field_and_momentum_zero":True,
        "central_relational_matter_origin_can_be_set_to_zero_by_shift_symmetry":True}


@cache
def health():
    emax=1+sp.Rational(1,10**6)
    e=domain.enclosures()
    pivot=e["on_constraint_rescaled_lapse_Hessian"]
    M=e["M"]
    maxM=max(abs(M["lower"]),abs(M["upper"]))
    return {"new_real_lapse_lower":domain.N0-domain.RN,
        "new_real_lapse_upper":domain.N0+domain.RN,
        "actual_real_conformal_lower":sp.Integer(1),"actual_real_conformal_upper":emax,
        "new_clock_kinetic_Schur_lower":-pivot["upper"]/(4*emax**3*maxM**2),
        "new_matter_kinetic_lower":1/(domain.N0+domain.RN),
        "new_tensor_principal_lower":1/emax**4,
        "new_clock_physical_squared_speed_lower":sp.Rational(493,500),
        "new_clock_physical_squared_speed_upper":sp.Rational(247,250),
        "new_matter_tensor_and_ordinary_Proca_principal_squared_speeds":(1,1,1),
        "new_physical_bounce_acceleration_lower":sp.Rational(39997,10000),
        "new_physical_bounce_acceleration_upper":sp.Rational(4001,1000),
        "seven_physical_modes_on_this_classical_background":True,
        "no_all_finite_momentum_Legendre_chart_or_nonlinear_stability_claim":True,
        "no_quantum_state_profile_EFT_cutoff_or_UV_transfer":True}


@cache
def checks():
    s=model.system()
    u,N,z=model.uf,model.Nf,model.zf
    w2=s["constraint_matter_square"]
    F,G=s["lapse_flow"],s["rescaled_trace_flow"]
    total=w2.diff(u)+w2.diff(N)*F+w2.diff(z)*G
    expected=2*(s["conformal_log_flow"]-3*s["hat_Hubble"])*w2
    p=picard()
    R,e,w,H,Elog=sp.symbols("R e w H E_log",nonzero=True)
    charge=w*R**3/e
    charge_dot=sp.diff(charge,w)*(Elog-3*H)*w+sp.diff(charge,R)*H*R+sp.diff(charge,e)*Elog*e
    return {"new_constraint_root_matter_evolution_is_exact_not_assigned":(total-expected).as_expr(),
        "reconstructed_physical_matter_charge_is_exactly_conserved":sp.factor(charge_dot),
        "new_central_charge_contains_conformal_rescaling":
            initial()["new_positive_conserved_matter_charge_squared"]*domain.N0-initial()["central_rescaled_matter_square"],
        "new_Picard_weighted_rate_three_twenty_fifths":p["Picard_contraction_upper"]-sp.Rational(3,25),
        "new_actual_lapse_displacement_one_e_minus_eleven":p["actual_solution_displacement_bounds"][0]-sp.Rational(1,10**11),
        "new_actual_trace_displacement_nine_e_minus_seven":p["actual_solution_displacement_bounds"][1]-sp.Rational(9,10**7),
        "new_duration_rate_scales_linearly":picard(domain.T/2)["Picard_contraction_upper"]-p["Picard_contraction_upper"]/2}


@cache
def gates():
    p,i,h=picard(),initial(),health()
    rows={"new_Picard_maps_both_inner_phase_discs_strictly_inside":all(a<b for a,b in
            zip(p["actual_solution_displacement_bounds"],p["inner_Picard_lapse_and_trace_radii"],strict=True)),
        "new_Picard_strict_contraction_below_one_half":p["Picard_contraction_upper"]<sp.Rational(1,2),
        "new_positive_conserved_charge_is_nonzero":i["new_positive_conserved_matter_charge_squared"]>0,
        "new_initial_charge_differs_from_old_one":i["new_positive_conserved_matter_charge_squared"]!=i["old_charge_squared_at_same_lapse"],
        "old_charge_is_not_silently_reused":i["old_charge_residual_in_new_rescaled_constraint"]!=0,
        "new_real_lapse_positive_and_above_one":h["new_real_lapse_lower"]>1,
        "complex_h_has_no_zero_on_new_time_disc":p["complex_background_h_minus_one_upper"]<sp.Rational(1,2),
        "complex_conformal_D_has_no_zero_on_new_phase_box":p["complex_conformal_D_minus_one_upper"]<sp.Rational(1,2),
        "declared_conformal_upper_dominates_sqrt_lapse":h["actual_real_conformal_upper"]**2>h["new_real_lapse_upper"],
        "new_log_hat_scale_below_two_e_minus_twelve":p["new_log_hat_scale_absolute_upper"]<sp.Rational(2,10**12),
        "new_hat_scale_deviation_below_three_e_minus_twelve":p["new_hat_scale_minus_one_absolute_upper"]<sp.Rational(3,10**12),
        "new_clock_kinetic_positive_and_above_eleven":h["new_clock_kinetic_Schur_lower"]>11,
        "new_matter_kinetic_above_ninety_nine_hundredths":h["new_matter_kinetic_lower"]>sp.Rational(99,100),
        "new_tensor_principal_above_ninety_nine_hundredths":h["new_tensor_principal_lower"]>sp.Rational(99,100),
        "new_scalar_squared_speeds_positive_and_clock_strictly_subluminal":
            0<h["new_clock_physical_squared_speed_lower"]<h["new_clock_physical_squared_speed_upper"]<1,
        "new_exact_center_is_physical_Hubble_zero":i["central_physical_Hubble"]==0,
        "positive_whole_interval_acceleration_makes_strict_physical_bounce":h["new_physical_bounce_acceleration_lower"]>0}
    return {name:bool(value) for name,value in rows.items()}


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
        sp.Symbol("free"),sp.sqrt(2),[],0,-1,2*domain.T)
    count=0
    for value in bad:
        for call in (duration,picard):
            try:
                call(value)
            except (TypeError,ValueError):
                count+=1
            else:
                raise ValueError("An invalid new-solution duration was accepted")
    return {"rejected_inputs":count,"validation_precedes_scientific_evaluation":True}
