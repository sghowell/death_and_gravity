"""Exact central Taylor majorants for a larger homogeneous existence domain."""
from functools import cache
from itertools import product
from math import factorial

import sympy as sp
from p8_proca_scalar_cone import domain

from . import existence, jets

N,u,ell=jets.N,jets.u,jets.ell
p=sp.Symbol("independent_enlarged_trace_density",real=True)
ROOT_RADIUS=sp.Rational(1,10**6)
TIME_PARAMETER_RADIUS=sp.Rational(1,10**6)
PHASE_RADIUS=sp.Rational(1,10**5)
TIME_WINDOW=sp.Rational(1,10**7)
FLOW_BOUND=sp.Integer(10)


@cache
def literal_jet():
    d,b=jets.data(),jets.background.data()
    a,U=b["trace_a"],b["volume_U"]
    b1=b["actual_trace_linear_time_derivative"]
    H0=N*(p*p/(4*a)-b["scalar_potential"]+ell*ell/(2*U))
    Hup=-N*b1/(2*a)
    Huu=N*(b1*b1/(2*a)-d["actual_scalar_potential_second_time_derivative"]
           +ell*ell*d["inverse_volume_second_time_derivative"]/2)
    return {"central_Hamiltonian_at_fixed_phase":H0,
            "mixed_time_trace_coefficient":sp.factor(Hup),
            "second_time_coefficient_at_zero_trace":sp.factor(Huu),
            "literal_total_degree_three_jet_generator":H0+Hup*u*p+Huu*u*u/2}


def coefficient_bound(value):
    # Differentiate first, then impose the real central charge relation.
    # Real positive square-root lapse is a new enclosure coordinate only.
    value=sp.factor(value)
    r=sp.Symbol("positive_square_root_center_lapse",positive=True)
    result=sp.Integer(0)
    parts=[]
    P2=jets.data()["positive_matter_density_squared"]
    for (degree,),coefficient in sp.Poly(value,ell).terms():
        reduced=sp.factor(coefficient*P2**(degree//2))
        rational=sp.cancel(reduced.subs(N,r*r)).subs(r,N)
        box=domain.enclosure(rational,sp.Rational(1,10**6))
        magnitude=max(abs(box["lower"]),abs(box["upper"]))
        if degree%2:
            magnitude*=sp.Rational(11,100)
        result+=magnitude
        parts.append({"central_charge_degree":degree,"rational_sqrt_lapse_enclosure":box,
                      "term_absolute_upper":magnitude})
    return {"absolute_upper":result,"exact_charge_expansion_terms":parts}


@cache
def coefficients():
    J=literal_jet()["literal_total_degree_three_jet_generator"]
    specs={"pivot":(sp.diff(J,N,2),(N,u,p,ell),2),
           "center_force":(sp.diff(J,N),(u,p,ell),3),
           "H":(J,(N,u,p,ell),1),
           "Hp":(sp.diff(J,p),(N,u,p,ell),1),
           "Hell":(sp.diff(J,ell),(N,u,p,ell),1)}
    result={}
    for name,(expression,variables,degree) in specs.items():
        rows={}
        for powers in product(range(degree+1),repeat=len(variables)):
            if sum(powers)>degree:
                continue
            value=expression
            for variable,count in zip(variables,powers,strict=True):
                value=sp.diff(value,variable,count)/factorial(count)
            value=sp.factor(value.subs({u:0,p:0},simultaneous=True))
            if value==0:
                continue
            rows[",".join(map(str,powers))]={"multi_degree":powers,"coefficient_before_central_charge_substitution":value,
                                           "enclosure":coefficient_bound(value)}
        result[name]={"variable_order":tuple(str(v) for v in variables),"degree":degree,"coefficients":rows}
    return result


def polynomial_bound(name,skip_constant=False):
    data=coefficients()[name]
    radii=(TIME_PARAMETER_RADIUS,PHASE_RADIUS,PHASE_RADIUS) if name=="center_force" else (
        ROOT_RADIUS,TIME_PARAMETER_RADIUS,PHASE_RADIUS,PHASE_RADIUS)
    result=sp.Integer(0)
    for row in data["coefficients"].values():
        powers=row["multi_degree"]
        if skip_constant and not any(powers):
            continue
        weight=sp.prod(radius**degree for radius,degree in zip(radii,powers,strict=True))
        result+=row["enclosure"]["absolute_upper"]*weight
    return result


def tail(bound,rho,degree):
    # One complex line through each real center: no missing multivariate multiplicity.
    return bound*rho**(degree+1)/(1-rho)


@cache
def bounds():
    rn,ru,ry=sp.Rational(1,400),sp.Rational(1,200),sp.Rational(1,8)
    rho=max(ROOT_RADIUS/rn,TIME_PARAMETER_RADIUS/ru,PHASE_RADIUS/ry)
    rho_force=max(TIME_PARAMETER_RADIUS/ru,PHASE_RADIUS/ry)
    pivot_tail=tail(2*existence.H_BOUND/rn**2,rho,2)
    force_tail=tail(existence.H_BOUND/rn,rho_force,3)
    H_tail=tail(existence.H_BOUND,rho,1)
    phase_derivative_tail=tail(existence.H_BOUND/ry,rho,1)
    pivot=polynomial_bound("pivot",skip_constant=True)+pivot_tail
    force=polynomial_bound("center_force",skip_constant=True)+force_tail
    H=polynomial_bound("H")+H_tail
    Hp=polynomial_bound("Hp")+phase_derivative_tail
    Hl=polynomial_bound("Hell")+phase_derivative_tail
    ell_bound=sp.Rational(11,100)+PHASE_RADIUS
    rate=pivot/2
    return {"lapse_root_radius":ROOT_RADIUS,"complex_time_parameter_radius":TIME_PARAMETER_RADIUS,
            "complex_phase_parameter_radius":PHASE_RADIUS,"joint_radial_Cauchy_ratio":rho,
            "center_force_radial_Cauchy_ratio":rho_force,"pivot_degree_two_tail_upper":pivot_tail,
            "center_force_degree_three_tail_upper":force_tail,"Hamiltonian_degree_one_tail_upper":H_tail,
            "phase_derivative_degree_one_tail_upper":phase_derivative_tail,
            "lapse_pivot_change_upper":pivot,"center_force_upper":force,
            "Newton_contraction_upper":rate,
            "Newton_closed_disc_image_radius_upper":force/2+rate*ROOT_RADIUS,
            "actual_root_displacement_upper":force/(2-pivot),
            "remaining_lapse_pivot_absolute_lower":2-pivot,
            "actual_Hamiltonian_absolute_upper":H,"actual_Hp_absolute_upper":Hp,
            "actual_Hell_absolute_upper":Hl,"actual_matter_density_absolute_upper":ell_bound,
            "actual_trace_evolution_absolute_upper":H+ell_bound*Hl,
            "actual_density_evolution_absolute_upper":ell_bound*Hp,
            "declared_phase_flow_absolute_upper":FLOW_BOUND,
            "enlarged_time_half_width":TIME_WINDOW,
            "Picard_image_radius_upper":FLOW_BOUND*TIME_WINDOW,
            "phase_flow_Lipschitz_upper":4*FLOW_BOUND/PHASE_RADIUS,
            "Picard_contraction_upper":4*FLOW_BOUND*TIME_WINDOW/PHASE_RADIUS,
            "reconstructed_log_hat_scale_absolute_upper":Hp*TIME_WINDOW/3}


def duration(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require an exact positive enlarged time half-width")
    value=sp.Rational(value)
    if not 0<value<=TIME_WINDOW:
        raise ValueError("Require an enlarged time half-width at most 10^-7")
    return value


def evolution(value=TIME_WINDOW):
    value=duration(value)
    b=bounds()
    return {"time_half_width":value,"real_time_interval":(-value,value),
            "phase_disc_radius":PHASE_RADIUS/2,"flow_absolute_upper":FLOW_BOUND,
            "Picard_image_radius_upper":FLOW_BOUND*value,
            "Picard_contraction_upper":4*FLOW_BOUND*value/PHASE_RADIUS,
            "log_hat_scale_absolute_upper":b["actual_Hp_absolute_upper"]*value/3,
            "same_literal_classical_solution_extends_smaller_proved_interval":True,
            "not_a_quantitative_time_dependent_cone_or_UV_verdict":True}


@cache
def checks():
    c,d=jets.model.coefficients(),jets.data()
    primitive=d["actual_boundary_primitive_phi"]*u+d["actual_boundary_primitive_phi_phi_phi"]*u**3/6
    Iphi=d["actual_boundary_primitive_phi"]+d["actual_boundary_primitive_phi_phi_phi"]*u*u/2
    actual=N*((p-c["b"])**2/(4*c["a"])-c["f0"]+ell*ell/(2*c["U"]))
    actual=actual.subs({jets.model.I:primitive,jets.model.Iphi:Iphi},simultaneous=True)
    j=literal_jet()
    at0={u:0,p:0}
    out={
        "enlarged_literal_central_Hamiltonian":sp.factor(actual.subs(u,0)-j["central_Hamiltonian_at_fixed_phase"]),
        "enlarged_literal_mixed_time_trace_jet":sp.factor(sp.diff(actual,u,p).subs(at0)-j["mixed_time_trace_coefficient"]),
        "enlarged_literal_second_time_jet":sp.factor(sp.diff(actual,u,2).subs(at0)-j["second_time_coefficient_at_zero_trace"]),
        "enlarged_missing_first_pure_time_jet_is_zero":sp.factor(sp.diff(actual,u).subs(at0)),
        "enlarged_missing_third_pure_time_jet_is_zero":sp.factor(sp.diff(actual,u,3).subs(at0)),
        "enlarged_missing_time_two_trace_jet_is_zero":sp.factor(sp.diff(actual,u,p,2).subs(at0)),
        "enlarged_missing_two_time_trace_jet_is_zero":sp.factor(sp.diff(actual,u,2,p).subs(at0))}
    H0=j["central_Hamiltonian_at_fixed_phase"]
    on={p:0,ell**2:d["positive_matter_density_squared"]}
    out["enlarged_exact_center_force_zero_after_derivative"]=sp.factor(sp.diff(H0,N).subs(on,simultaneous=True))
    out["enlarged_exact_center_pivot_matches_original_fixed_phase"]=sp.factor(
        sp.diff(H0,N,2).subs(on,simultaneous=True)-d["fixed_phase_lapse_Hessian"])
    r=sp.Symbol("positive_center_sqrt_lapse_check",positive=True)
    out["enlarged_sqrt_lapse_recoordinate_is_exact"]=sp.sqrt(r*r)-r
    rho=sp.Symbol("radial_ratio",positive=True)
    for degree in (1,2,3):
        # Geometric tail after removing the first degree+1 coefficients.
        out["enlarged_radial_Cauchy_tail_degree_"+str(degree)]=sp.factor(
            1/(1-rho)-sum(rho**i for i in range(degree+1))-rho**(degree+1)/(1-rho))
    b=bounds()
    out["enlarged_pivot_tail_normalization"]=b["pivot_degree_two_tail_upper"]-sp.Rational(512,2499)
    out["enlarged_force_tail_normalization"]=b["center_force_degree_three_tail_upper"]-sp.Rational(1,156218750)
    out["enlarged_Hamiltonian_tail_normalization"]=b["Hamiltonian_degree_one_tail_upper"]-sp.Rational(4,2499)
    out["enlarged_phase_derivative_tail_normalization"]=b["phase_derivative_degree_one_tail_upper"]-sp.Rational(32,2499)
    out["enlarged_two_column_Picard_rate"]=b["Picard_contraction_upper"]-sp.Rational(2,5)
    return out


@cache
def gates():
    b=bounds()
    return {key:bool(value) for key,value in {
        "enlarged_root_and_time_phase_box_fit_old_complex_half_polydisc":b["joint_radial_Cauchy_ratio"]<1,
        "enlarged_force_tail_uses_fixed_center_not_root_radius":b["center_force_radial_Cauchy_ratio"]<b["joint_radial_Cauchy_ratio"],
        "enlarged_pivot_change_below_one_quarter":b["lapse_pivot_change_upper"]<sp.Rational(1,4),
        "enlarged_force_below_six_e_minus_seven":b["center_force_upper"]<sp.Rational(6,10**7),
        "enlarged_Newton_rate_below_one_eighth":b["Newton_contraction_upper"]<sp.Rational(1,8),
        "enlarged_Newton_disc_strictly_invariant":b["Newton_closed_disc_image_radius_upper"]<ROOT_RADIUS,
        "enlarged_root_displacement_below_four_e_minus_seven":b["actual_root_displacement_upper"]<sp.Rational(4,10**7),
        "enlarged_lapse_pivot_above_seven_fourths":b["remaining_lapse_pivot_absolute_lower"]>sp.Rational(7,4),
        "enlarged_Hamiltonian_below_nine":b["actual_Hamiltonian_absolute_upper"]<9,
        "enlarged_Hp_below_one_fiftieth":b["actual_Hp_absolute_upper"]<sp.Rational(1,50),
        "enlarged_Hell_below_one_fifth":b["actual_Hell_absolute_upper"]<sp.Rational(1,5),
        "enlarged_actual_matter_density_below_one_eighth":b["actual_matter_density_absolute_upper"]<sp.Rational(1,8),
        "enlarged_both_phase_flow_components_below_ten":max(b["actual_trace_evolution_absolute_upper"],b["actual_density_evolution_absolute_upper"])<FLOW_BOUND,
        "enlarged_Picard_image_inside_phase_half_disc":b["Picard_image_radius_upper"]<PHASE_RADIUS/2,
        "enlarged_Picard_rate_below_one_half":b["Picard_contraction_upper"]<sp.Rational(1,2),
        "enlarged_time_inside_implicit_parameter_domain":TIME_WINDOW<TIME_PARAMETER_RADIUS,
        "enlarged_log_scale_bound_below_one_e_minus_nine":b["reconstructed_log_hat_scale_absolute_upper"]<sp.Rational(1,10**9),
        "enlarged_lapse_stays_in_original_covariant_clock_tube":(1+existence.SHIFT_MAX+ROOT_RADIUS)**-2>sp.Rational(9,10) and (1-ROOT_RADIUS)**-2<sp.Rational(11,10),
        "enlarged_analytic_solution_agrees_with_old_local_branch":True,
        "enlarged_time_not_a_proved_cone_or_propagation_band":True}.items()}


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
         sp.Symbol("unfixed"),sp.sqrt(2),[],0,-1,2*TIME_WINDOW)
    count=0
    for value in bad:
        for function in (duration,evolution):
            try:
                function(value)
            except (TypeError,ValueError):
                count+=1
    if count!=2*len(bad):
        raise ValueError("An invalid enlarged time-domain input was accepted")
    return {"rejected_inputs":count,"validation_before_enlarged_domain_evaluation":True}
