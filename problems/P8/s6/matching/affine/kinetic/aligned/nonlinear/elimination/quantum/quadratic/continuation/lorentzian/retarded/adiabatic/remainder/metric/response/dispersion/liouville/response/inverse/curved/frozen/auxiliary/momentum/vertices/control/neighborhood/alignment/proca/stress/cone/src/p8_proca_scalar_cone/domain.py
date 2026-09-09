"""Exact continuous central-slice cone domain and its finite-amplitude boundary."""
from functools import cache

import sympy as sp
from p8_coupled_energy.bounds import exact_expression
from p8_offclock_scalar import datum

from . import background, principal

N=background.N
RADIUS=sp.Rational(1,10**8)
OUTER_RADIUS=sp.Rational(1,10**6)


def rational(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require an exact native or SymPy rational")
    return sp.Rational(value)


def exact_radius(value):
    value=rational(value)
    if not 0<value<=OUTER_RADIUS:
        raise ValueError("Require a positive central radius at most 10^-6")
    return value


def enclosure(expression,radius):
    expression=exact_expression(expression,(N,))
    radius=exact_radius(radius)
    z=sp.Symbol("central_lapse_displacement",real=True)
    num,den=sp.fraction(sp.cancel(expression))
    a,b=(sp.Poly(value.subs(N,1+z),z,domain=sp.QQ) for value in (num,den))
    if b.nth(0)<0:
        a,b=-a,-b
    tail=lambda p:sum(abs(c)*radius**degree[0] for degree,c in p.terms() if degree[0])
    at,bt=tail(a),tail(b)
    lower,upper=b.nth(0)-bt,b.nth(0)+bt
    if lower<=0:
        raise ValueError("No nonzero denominator on this continuous central domain")
    endpoints=((a.nth(0)-at)/lower,(a.nth(0)-at)/upper,
               (a.nth(0)+at)/lower,(a.nth(0)+at)/upper)
    return {"lower":min(endpoints),"upper":max(endpoints),"radius":radius,
            "numerator_center":a.nth(0),"numerator_tail_upper":at,
            "denominator_center":b.nth(0),"denominator_tail_upper":bt}


def domain(radius=RADIUS):
    return _domain(exact_radius(radius))


@cache
def _domain(radius):
    d,f=background.data(),principal.formula()
    expressions={"matter_momentum_squared":d["matter_momentum_squared"],
                 "lapse_Hessian_times_N_seven_halves":sp.factor(d["fixed_phase_lapse_Hessian"]*N**sp.Rational(7,2)),
                 "clock_gradient_Schur_over_sqrt_N":sp.factor(f["clock_gradient_Schur"]/sp.sqrt(N)),
                 "clock_physical_speed_squared":f["clock_physical_speed_squared"],
                 "physical_clock_cone_margin":f["physical_clock_cone_margin"],
                 "physical_clock_norm_above_original_lower_tube":N**-2-sp.Rational(9,10),
                 "physical_clock_norm_below_original_upper_tube":sp.Rational(11,10)-N**-2,
                 "original_affine_complement_rank_factor":2/N**2-1,
                 "old_affine_temporal_update_factor":datum.data()["gamma_t"],
                 "old_affine_spatial_update_factor":datum.data()["gamma_s"],
                 "positive_physical_tensor_kinetic":N**-2}
    boxes={name:enclosure(value,radius) for name,value in expressions.items()}
    gates={"positive_matter_root":boxes["matter_momentum_squared"]["lower"]>sp.Rational(1,200),
           "negative_fixed_phase_lapse_Hessian":boxes["lapse_Hessian_times_N_seven_halves"]["upper"]<-2,
           "positive_scalar_gradient_Schur":boxes["clock_gradient_Schur_over_sqrt_N"]["lower"]>11,
           "positive_clock_speed":boxes["clock_physical_speed_squared"]["lower"]>sp.Rational(9,10),
           "strictly_subluminal_clock_on_this_whole_interval":boxes["physical_clock_cone_margin"]["lower"]>sp.Rational(1,10**6),
           "inside_original_covariant_clock_norm_tube":all(boxes[name]["lower"]>0 for name in
               ("physical_clock_norm_above_original_lower_tube","physical_clock_norm_below_original_upper_tube")),
           "old_complement_and_new_affine_mass_update_invertible":all(boxes[name]["lower"]>sp.Rational(9,10) for name in
               ("original_affine_complement_rank_factor","old_affine_temporal_update_factor","old_affine_spatial_update_factor")),
           "positive_luminal_tensor_block":boxes["positive_physical_tensor_kinetic"]["lower"]>sp.Rational(9,10)}
    return {"central_lapse_radius":radius,"continuous_rational_enclosures":boxes,
            "gates":{name:bool(value) for name,value in gates.items()},
            "not_the_old_nine_invariant_interaction_polydisc":True}


def point(lapse):
    lapse=rational(lapse)
    if abs(lapse-1)>OUTER_RADIUS:
        raise ValueError("Require the declared central interval |N-1|<=10^-6")
    d,f=background.data(),principal.formula()
    values={"lapse":lapse,"positive_matter_momentum_squared":d["matter_momentum_squared"].subs(N,lapse),
            "fixed_phase_lapse_Hessian":d["fixed_phase_lapse_Hessian"].subs(N,lapse),
            "clock_physical_speed_squared":f["clock_physical_speed_squared"].subs(N,lapse),
            "clock_cone_margin":f["physical_clock_cone_margin"].subs(N,lapse),
            "clock_gradient_Schur":f["clock_gradient_Schur"].subs(N,lapse)}
    return {name:sp.factor(value) for name,value in values.items()}


@cache
def boundary():
    f=principal.formula()
    polynomial=1124999*N**4+7500000*N**2-8625003
    lo=1+sp.Rational(2,10**7)
    hi=1+sp.Rational(21,10**8)
    p=point(1+OUTER_RADIUS)
    derivative_box=enclosure(sp.diff(polynomial,N),OUTER_RADIUS)
    return {"exact_clock_luminality_boundary_polynomial":polynomial,
            "root_bracket":(lo,hi),"polynomial_at_lower":polynomial.subs(N,lo),
            "polynomial_at_upper":polynomial.subs(N,hi),
            "positive_polynomial_derivative_enclosure":derivative_box,
            "actual_superluminal_fixed_phase_constraint_datum":p,
            "control_is_ghost_and_gradient_positive":bool(p["positive_matter_momentum_squared"]>0
                and p["fixed_phase_lapse_Hessian"]<0 and p["clock_gradient_Schur"]>0),
            "control_is_strictly_superluminal":bool(p["clock_physical_speed_squared"]>1),
            "same_positive_mass_source_free_action_not_old_q_squared_vector_mixing":True,
            "denominator_on_outer_interval":enclosure(
                16874985*N**4-14244998*N**2-8625003,OUTER_RADIUS),
            "speed_formula":f["clock_physical_speed_squared"],
            "no_finite_frequency_EFT_or_global_UV_verdict":True}


@cache
def checks():
    f=principal.formula()
    numerator=12374989*N**4-44244998*N**2+25875009
    denominator=16874985*N**4-14244998*N**2-8625003
    p=boundary()["exact_clock_luminality_boundary_polynomial"]
    return {"explicit_actual_clock_speed_rational_function":sp.factor(f["clock_physical_speed_squared"]-numerator/denominator),
            "exact_clock_luminality_boundary_polynomial":sp.factor(f["physical_clock_cone_margin"]-4*p/denominator)}


@cache
def gates():
    d,b=domain(),boundary()
    out=dict(d["gates"])
    out.update({"simple_luminality_boundary_bracket":bool(b["polynomial_at_lower"]<0<b["polynomial_at_upper"]
                   and b["positive_polynomial_derivative_enclosure"]["lower"]>0),
                "no_pole_in_boundary_or_superluminal_control_interval":b["denominator_on_outer_interval"]["upper"]<0,
                "actual_control_remains_ghost_and_gradient_positive":b["control_is_ghost_and_gradient_positive"],
                "actual_control_is_superluminal":b["control_is_strictly_superluminal"],
                "only_central_classical_homogeneous_family_not_full_spatial_neighborhood":True,
                "broader_domain_not_assigned_old_finite_tree_or_auxiliary_polydisc_bounds":True,
                "fixed_quantum_profile_and_functional_response_not_in_classical_cone":True,
                "original_P8_and_V_G_B_still_open":True})
    return {name:bool(value) for name,value in out.items()}
