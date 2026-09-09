"""Full central clock-tube family and the weaker-margin negative control."""
from functools import cache

import sympy as sp

from . import domain, model

E=sp.Symbol("total_covariant_margin",positive=True)
y=sp.Symbol("central_lapse_squared",positive=True)
Y_LO,Y_HI,Y_CAP=sp.Rational(10,11),sp.Rational(10,9),sp.Rational(501,500)
CONTROL_LAPSE=sp.Rational(10001,10000)


@cache
def derived():
    c=model.old.coefficients()
    N,u=model.N,model.u
    V=c["V0"]+(E-model.OLD_MARGIN)*(N*N-1)**2/(N*c["D"]*c["h"])
    C=sp.factor(-sp.diff(V,N)+c["Ln"]*V+sp.diff(c["I0"],u)-c["Lt"]*c["I0"])
    at=lambda value:sp.factor(value.subs(u,0))
    M,A,a0,Lnu=(at(c[name]) for name in ("M","A","a0","Ln"))
    Ltu,Bu=at(sp.diff(c["Lt"],u)),at(sp.diff(c["B"],u))
    W=sp.factor(-at(C)/M)
    P=sp.factor(at(sp.diff(C,N))+(at(sp.diff(c["M"],N))+2*M*Lnu)*W)
    Z=sp.factor(at(V-sp.diff(c["d0"],u))+at(c["m0"])*W)
    G=sp.factor(-Bu/(3*M)-3*N*N*Z-N*W)
    speed=sp.factor(-G/(4*N**4*P))
    Fu=sp.factor(at(sp.diff(c["M"],u,2))*W+at(sp.diff(C,u,2))+Bu*Z+2*M*Ltu*W)
    Fz=sp.factor(Bu+2*A*Z-4*M*a0*W)
    Ndd=sp.factor(-(Fu+Fz*Z)/P)
    acceleration=sp.factor((-N*Z/2+Ltu+Lnu*Ndd)/N**2)
    return {"rescaled_matter_square":W,"fixed_phase_rescaled_lapse_pivot":P,
        "actual_rescaled_trace_velocity":Z,"actual_rescaled_clock_gradient":G,
        "clock_physical_squared_speed":speed,"actual_lapse_second_time_derivative":Ndd,
        "physical_Hubble_proper_time_acceleration":acceleration}


def exact(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require exact rational central-family inputs")
    return sp.Rational(value)


def diagnostic(margin=model.NEW_MARGIN,lapse=CONTROL_LAPSE):
    margin,lapse=exact(margin),exact(lapse)
    if not 0<margin<=sp.Rational(1,100):
        raise ValueError("Require total margin in (0,1/100]")
    if lapse<=0 or not sp.Rational(9,10)<=lapse**-2<=sp.Rational(11,10):
        raise ValueError("Require positive lapse in the declared physical clock tube")
    rows={name:sp.factor(value.subs({E:margin,model.N:lapse})) for name,value in derived().items()}
    return {"total_margin":margin,"central_lapse":lapse,**rows,
        "positive_matter_charge_possible":bool(rows["rescaled_matter_square"]>0),
        "positive_clock_kinetic_pivot":bool(rows["fixed_phase_rescaled_lapse_pivot"]<0),
        "strict_classical_bounce_acceleration":bool(rows["physical_Hubble_proper_time_acceleration"]>0),
        "fast_physical_clock":bool(rows["clock_physical_squared_speed"]>1),
        "diagnostic_does_not_assign_a_real_state_to_negative_matter_square":True}


def bernstein(poly,lower=Y_LO,upper=Y_CAP):
    t=sp.Symbol("Bernstein_coordinate",real=True)
    p=sp.Poly(sp.expand(poly.subs(y,lower+(upper-lower)*t)),t,domain=sp.QQ)
    degree=p.degree()
    coefficients=tuple(sp.factor(sum(p.nth(j)*sp.binomial(k,j)/sp.binomial(degree,j)
        for j in range(k+1))) for k in range(degree+1))
    reconstructed=sum(coefficients[k]*sp.binomial(degree,k)*t**k*(1-t)**(degree-k)
        for k in range(degree+1))
    return {"degree":degree,"coefficients":coefficients,"lower":min(coefficients),
        "upper":max(coefficients),"basis_reconstruction_residual":sp.expand(reconstructed-p.as_expr())}


@cache
def theorem():
    p=-1120*y*y+2847*y-1728
    D=-3360*y*y+2847*y+1728
    B=-2464*y*y+8847*y-5184
    A=-sp.Rational(224,25)*y*y-60*y+sp.Rational(1728,25)
    pa=-(75264*y**4-845488*y**3+2124921*y**2-1708562*y+313365)
    qa=25*y*(-1120*y*y+949*y+576)
    ba,bb=bernstein(pa),bernstein(qa)
    return {"physical_clock_tube":"9/10<=X=N^-2<=11/10; real N>0",
        "central_lapse_squared_tube":(Y_LO,Y_HI),
        "strict_positive_matter_branch_lapse_squared_upper":Y_CAP,
        "minus_one_hundred_times_matter_square":p,
        "positive_clock_speed_denominator":D,
        "positive_clock_speed_numerator":B,
        "one_hundredth_times_cone_margin_numerator":A,
        "matter_polynomial_derivative_lower_on_full_tube":sp.diff(p,y).subs(y,Y_HI),
        "matter_polynomial_positive_at_branch_cap":p.subs(y,Y_CAP),
        "speed_denominator_derivative_upper_on_allowed_branch":sp.diff(D,y).subs(y,Y_LO),
        "speed_numerator_derivative_lower_on_allowed_branch":sp.diff(B,y).subs(y,Y_CAP),
        "cone_margin_numerator_derivative_upper":sp.diff(A,y).subs(y,Y_LO),
        "speed_denominator_lower":D.subs(y,Y_CAP),"speed_denominator_upper":D.subs(y,Y_LO),
        "speed_numerator_lower":B.subs(y,Y_LO),
        "positive_clock_squared_speed_lower":B.subs(y,Y_LO)/D.subs(y,Y_LO),
        "positive_clock_squared_cone_margin_lower":100*A.subs(y,Y_CAP)/D.subs(y,Y_LO),
        "actual_central_acceleration_numerator":pa,"actual_central_acceleration_denominator":qa,
        "actual_acceleration_numerator_Bernstein_certificate":ba,
        "actual_acceleration_denominator_Bernstein_certificate":bb,
        "positive_actual_central_acceleration_lower":ba["lower"]/bb["upper"],
        "weaker_one_e_minus_four_margin_counterexample":diagnostic(sp.Rational(1,10000)),
        "selected_margin_at_same_counterexample_lapse":diagnostic(),
        "whole_continuum_not_a_point_scan":True,
        "central_family_result_not_whole_off_center_phase_tube_or_global_nearby_completion":True}


@cache
def checks():
    d=derived()
    n=model.N
    W=-(1000*n**4*E-1125*n**4-400*n*n*E+2849*n*n-600*E-1725)/100
    D=3000*n**4*E-3375*n**4-400*n*n*E+2849*n*n+600*E+1725
    B=2200*n**4*E-2475*n**4-400*n*n*E+8849*n*n-1800*E-5175
    A=8*n**4*E-9*n**4-60*n*n+24*E+69
    s=model.system()
    center=lambda value:sp.factor(value.evaluate([(model.uf,0),(model.zf,0)]).as_expr())
    t=theorem()
    mapping={y:n*n}
    H=domain.physical_hubble()
    hp,hq=H.numer.as_expr(),H.denom.as_expr()
    point={model.u:0,model.z:0}
    direct_acceleration=sp.factor((sp.diff(hp,model.u).subs(point)
        +sp.diff(hp,model.z).subs(point)*center(s["rescaled_trace_flow"]))/(hq.subs(point)*n))
    rows={"generic_central_matter_square_from_literal_retuned_constraint":d["rescaled_matter_square"]-W,
        "generic_central_fixed_phase_pivot_before_family_derivative":d["fixed_phase_rescaled_lapse_pivot"]+D/(400*n**3),
        "generic_central_speed_keeps_actual_moving_gradient":d["clock_physical_squared_speed"]-B/D,
        "generic_central_cone_margin_numerator_is_exact":1-d["clock_physical_squared_speed"]-100*A/D,
        "selected_central_matter_polynomial_matches_new_action":
            t["minus_one_hundred_times_matter_square"].subs(mapping)+100*center(s["constraint_matter_square"]),
        "selected_central_speed_matches_full_new_rational_flow":
            (t["positive_clock_speed_numerator"]/t["positive_clock_speed_denominator"]).subs(mapping)
            -center(s["clock_physical_speed_squared"]),
        "selected_central_acceleration_polynomials_match_literal_time_jets":
            (t["actual_central_acceleration_numerator"]/t["actual_central_acceleration_denominator"]).subs(mapping)
            -d["physical_Hubble_proper_time_acceleration"].subs(E,model.NEW_MARGIN),
        "central_second_jet_reduction_agrees_with_complete_new_physical_Hubble_chain_rule":
            direct_acceleration-d["physical_Hubble_proper_time_acceleration"].subs(E,model.NEW_MARGIN),
        "actual_acceleration_numerator_Bernstein_reconstruction":
            t["actual_acceleration_numerator_Bernstein_certificate"]["basis_reconstruction_residual"],
        "actual_acceleration_denominator_Bernstein_reconstruction":
            t["actual_acceleration_denominator_Bernstein_certificate"]["basis_reconstruction_residual"],
        "matter_square_at_original_clock_is_independent_of_margin":d["rescaled_matter_square"].subs(n,1)-sp.Rational(1,100),
        "central_branch_upper_is_stricter_than_full_clock_tube":Y_CAP-sp.Rational(501,500)}
    return {name:sp.factor(value) for name,value in rows.items()}


@cache
def gates():
    t=theorem()
    weak=t["weaker_one_e_minus_four_margin_counterexample"]
    strong=t["selected_margin_at_same_counterexample_lapse"]
    rows={"full_central_clock_tube_contains_positive_branch_cap":Y_LO<1<Y_CAP<Y_HI,
        "matter_polynomial_strictly_increasing_on_entire_central_clock_tube":t["matter_polynomial_derivative_lower_on_full_tube"]>0,
        "positive_matter_implies_strict_branch_cap":t["matter_polynomial_positive_at_branch_cap"]>0,
        "clock_denominator_decreases_and_stays_positive":t["speed_denominator_derivative_upper_on_allowed_branch"]<0
            and t["speed_denominator_lower"]>0,
        "clock_numerator_increases_and_stays_positive":t["speed_numerator_derivative_lower_on_allowed_branch"]>0
            and t["speed_numerator_lower"]>0,
        "clock_margin_numerator_decreases":t["cone_margin_numerator_derivative_upper"]<0,
        "all_positive_matter_central_clock_squared_speeds_above_one_half":t["positive_clock_squared_speed_lower"]>sp.Rational(1,2),
        "all_positive_matter_central_clock_squared_margins_above_one_over_four_thousand":
            t["positive_clock_squared_cone_margin_lower"]>sp.Rational(1,4000),
        "actual_central_acceleration_numerator_has_positive_Bernstein_coefficients":
            t["actual_acceleration_numerator_Bernstein_certificate"]["lower"]>0,
        "actual_central_acceleration_denominator_has_positive_Bernstein_coefficients":
            t["actual_acceleration_denominator_Bernstein_certificate"]["lower"]>0,
        "all_positive_matter_central_bounce_accelerations_above_three":t["positive_actual_central_acceleration_lower"]>3,
        "weaker_margin_counterexample_has_positive_matter_regular_kinetic_and_bounce":
            weak["positive_matter_charge_possible"] and weak["positive_clock_kinetic_pivot"]
            and weak["strict_classical_bounce_acceleration"],
        "weaker_margin_counterexample_really_has_fast_clock":weak["fast_physical_clock"],
        "selected_margin_repairs_same_lapse_with_its_own_positive_charge":
            strong["positive_matter_charge_possible"] and strong["positive_clock_kinetic_pivot"]
            and strong["strict_classical_bounce_acceleration"] and not strong["fast_physical_clock"],
        "changing_margin_changes_the_positive_charge_not_just_the_speed":
            strong["rescaled_matter_square"]!=weak["rescaled_matter_square"]}
    return {name:bool(value) for name,value in rows.items()}


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
        sp.Symbol("free"),sp.sqrt(2),[])
    pairs=[pair for value in bad for pair in ((value,1),(model.NEW_MARGIN,value))]
    pairs.extend(((0,1),(-1,1),(1,1),(model.NEW_MARGIN,0),
        (model.NEW_MARGIN,-1),(model.NEW_MARGIN,2),(model.NEW_MARGIN,sp.Rational(1,2))))
    count=0
    for pair in pairs:
        try:
            diagnostic(*pair)
        except (TypeError,ValueError):
            count+=1
        else:
            raise ValueError("An invalid central-family input was accepted")
    return {"rejected_inputs":count,"negative_matter_square_is_diagnostic_not_real_state_assignment":True}
