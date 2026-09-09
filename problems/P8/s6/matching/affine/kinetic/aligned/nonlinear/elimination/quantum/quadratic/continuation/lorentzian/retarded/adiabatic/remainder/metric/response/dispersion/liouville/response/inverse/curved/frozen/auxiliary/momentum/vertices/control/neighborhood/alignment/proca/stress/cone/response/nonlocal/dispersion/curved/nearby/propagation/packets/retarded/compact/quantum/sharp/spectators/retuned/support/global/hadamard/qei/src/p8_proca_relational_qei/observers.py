"""Actual physical observer cones and named all-time admissible worldlines."""
from functools import cache

import sympy as sp
from p8_proca_global_hadamard import domains
from p8_proca_global_support import model

SAFE_SPEED=sp.Rational(99,100)
CLOCK_LOWER=sp.Rational(993,1000)
INTER_SPEED=sp.Rational(999,1000)
INTER_RADIUS=sp.Rational(1,100)


def exact_speed(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require an exact rational physical observer speed")
    value=sp.Rational(value)
    if not 0<=value<1:
        raise ValueError("The physical observer must be future timelike")
    return value


def curve(value):
    value=exact_speed(value)
    F=model.background()["physical_matter_null_primitive"]
    return {"physical_speed":value,"coordinate_path_along_first_axis":value*F,
        "proper_time_as_function_of_u":sp.sqrt(1-value**2)*model.u,
        "named_uniformly_subclock_class":bool(value<=SAFE_SPEED)}


@cache
def data():
    d=domains.data()
    upper=sp.factor(d["exact_J_times_h_squared_numerator"].subs(model.u,INTER_RADIUS)/800)
    c2upper=upper/(upper+4*model.EPSILON)
    c,a,k,udot,v=sp.symbols("positive_clock_speed positive_scale positive_covector_radius positive_time_tangent positive_observer_speed",positive=True)
    nx,ny,nz=sp.symbols("direction_x direction_y direction_z",real=True)
    normal=sp.Matrix([c/v,sp.sqrt(1-c*c/(v*v)),0])
    return {"safe_global_physical_speed_upper":SAFE_SPEED,"global_clock_speed_rational_lower":CLOCK_LOWER,
        "uniform_signed_contraction_margin":CLOCK_LOWER-SAFE_SPEED,
        "general_positive_mode_covector_contraction":"u_dot*k/a * (physical_velocity dot direction - mode_speed)",
        "interluminal_named_physical_speed":INTER_SPEED,"interluminal_certified_time_radius":INTER_RADIUS,
        "interluminal_interval_J_h_squared_upper":upper,
        "interluminal_interval_clock_squared_speed_upper":c2upper,
        "explicit_interluminal_clock_conormal_direction":normal,
        "conormal_symbols":{"clock_speed":c,"scale":a,"radius":k,"time_tangent":udot,"observer_speed":v},
        "generic_covector_direction":sp.Matrix([nx,ny,nz]),
        "conormal_failure_is_not_by_itself_a_nonexistence_or_QEI_no_go":True,
        "zero_clock_matter_contraction_direction_exists_when_observer_speed_exceeds_clock":True}


@cache
def checks():
    d=data()
    c,a,k,udot,v=(d["conormal_symbols"][name] for name in
        ("clock_speed","scale","radius","time_tangent","observer_speed"))
    n=d["explicit_interluminal_clock_conormal_direction"]
    vx,vy,vz,nx,ny,nz=sp.symbols("vx vy vz nx ny nz",real=True)
    vel,direction=sp.Matrix([vx,vy,vz]),sp.Matrix([nx,ny,nz])
    speed=sp.Symbol("physical_speed",real=True)
    path=speed*model.background()["physical_matter_null_primitive"]
    scale=model.background()["a"]
    return {
        "real_Cauchy_Schwarz_observer_margin_is_Gram_identity":sp.expand(
            vel.dot(vel)*direction.dot(direction)-vel.dot(direction)**2-vel.cross(direction).dot(vel.cross(direction))),
        "named_constant_physical_speed_curve_has_actual_metric_norm":sp.factor(
            -1+scale**2*sp.diff(path,model.u)**2+1-speed**2),
        "named_safe_signed_contraction_margin_is_three_over_one_thousand":d["uniform_signed_contraction_margin"]-sp.Rational(3,1000),
        "interluminal_clock_conormal_direction_is_unit":sp.simplify(n.dot(n)-1),
        "interluminal_clock_covector_annihilates_timelike_curve":sp.factor(
            udot*k/a*(-c+v*n[0])),
        "actual_global_clock_squared_speed_lower_replayed":domains.data()["global_clock_squared_speed_lower"]-sp.Rational(1199,1215),
    }


@cache
def gates():
    d=data()
    return {name:bool(value) for name,value in {
        "named_uniform_observers_are_physically_timelike":0<=SAFE_SPEED<1,
        "actual_global_clock_speed_above_named_rational_lower":CLOCK_LOWER**2<sp.Rational(1199,1215),
        "named_uniform_observers_are_strictly_slower_than_clock_everywhere":SAFE_SPEED<CLOCK_LOWER,
        "named_interluminal_observer_is_physically_timelike":0<INTER_SPEED<1,
        "whole_named_interluminal_interval_is_inside_gamma_slab":0<INTER_RADIUS<sp.Rational(1,4),
        "whole_named_interluminal_interval_clock_slower_than_observer":d["interluminal_interval_clock_squared_speed_upper"]<INTER_SPEED**2,
        "actual_interluminal_bound_is_strictly_positive":d["interluminal_interval_clock_squared_speed_upper"]>0,
    }.items()}


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
        sp.Symbol("free"),sp.sqrt(2),[],1,2,-1)
    rejected=0
    for value in bad:
        for call in (exact_speed,curve):
            try:
                call(value)
            except (TypeError,ValueError):
                rejected+=1
    if rejected!=2*len(bad):
        raise ValueError("An inadmissible physical worldline speed was accepted")
    return {"rejected_inputs":rejected}
