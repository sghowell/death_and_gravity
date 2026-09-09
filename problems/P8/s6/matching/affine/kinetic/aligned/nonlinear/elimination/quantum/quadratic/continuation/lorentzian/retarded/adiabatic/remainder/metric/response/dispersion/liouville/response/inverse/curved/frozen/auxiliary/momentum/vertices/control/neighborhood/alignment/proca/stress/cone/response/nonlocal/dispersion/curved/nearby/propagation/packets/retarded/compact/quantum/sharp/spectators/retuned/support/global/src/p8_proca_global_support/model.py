"""The actual global clock and universal complete scalar coefficient calculus."""
from functools import cache

import sympy as sp
from p8_aligned_margin import dynamics as margin
from p8_proca_retuned_margin import model as actual

old=margin.old
u,q,a,H=old.u,old.q,old.a,old.H
k=sp.Symbol("nonzero_comoving_scalar_root",nonzero=True)
EPSILON=actual.NEW_MARGIN
COEFFICIENTS=(old.H,old.theta,old.lam,old.ell,old.w,old.J,margin.delta_J)
JETS={value:sp.Symbol(str(value)+"_actual_time_derivative",real=True) for value in COEFFICIENTS}


def dt(value):
    if isinstance(value,sp.MatrixBase):
        return value.applyfunc(dt)
    return sp.diff(value,a)*H*a-2*H*q*sp.diff(value,q)+sum(
        sp.diff(value,key)*jet for key,jet in JETS.items())


@cache
def background():
    d=dict(old.background())
    d["new_total_margin"]=EPSILON
    d["delta_J"]=4*EPSILON/d["h"]**2
    d["J_new"]=d["J"]+d["delta_J"]
    d["clock_physical_squared_speed"]=d["J"]/d["J_new"]
    d["physical_matter_null_primitive"]=sp.atan(u)/2+u/(2*(1+u*u))
    return d


@cache
def substitution():
    d=background()
    sub={getattr(old,name):d[name] for name in ("a","H","theta","lam","ell","w","J")}
    sub[margin.delta_J]=d["delta_J"]
    sub.update({jet:sp.diff(sub[key],u) for key,jet in JETS.items()})
    return sub


@cache
def checks():
    d=background()
    return {"retuned_global_margin_is_the_new_literal_total_not_old_margin":EPSILON-sp.Rational(1,200),
        "global_clock_scale_is_the_original_complete_bounce":d["a"]-(1+u*u)**2,
        "global_clock_Theta_has_only_the_zero_at_the_bounce":sp.factor(d["theta"]-u*(4-1/d["h"])/(1+u*u)),
        "global_clock_Lambda_has_the_explicit_regular_gamma_strip_form":sp.factor(d["lam"]-1+3/(2*d["h"])),
        "global_clock_matter_density_has_the_all_time_upper_one_tenth":sp.factor(d["ell"]-1/(10*(1+u*u)**6)),
        "actual_global_physical_matter_null_primitive_derivative":sp.factor(
            sp.diff(d["physical_matter_null_primitive"],u)-1/d["a"]),
        "global_clock_matter_mixing_identity":sp.factor(d["w"]+d["ell"]*d["lam"]),
        "global_clock_full_moving_gradient_identity":sp.factor(d["J"]+d["w"]**2/2
            -d["theta"]*(d["H"]*d["lam"]+sp.diff(d["lam"],u))
            +d["lam"]*sp.diff(d["theta"],u)+d["theta"]**2)}
