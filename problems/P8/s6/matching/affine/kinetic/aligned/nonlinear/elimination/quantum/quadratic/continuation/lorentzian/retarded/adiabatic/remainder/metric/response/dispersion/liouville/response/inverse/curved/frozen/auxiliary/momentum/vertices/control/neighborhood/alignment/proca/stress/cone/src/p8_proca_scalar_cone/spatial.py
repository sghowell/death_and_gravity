"""Independent full spatial Hessian to first order in central background time."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model
from p8_coupled_momentum import constraints
from p8_offclock_scalar import spatial as old
from p8_physical import geometry
from p8_physical import jets as j
from sympy.polys.domains import QQ_I

from . import background

t,k,P=old.t,old.k,old.P
eta=sp.Symbol("formal_central_time",real=True)
pd=background.pd
CHANNELS=("curvature","scalar_p","matter","matter_p")


def clean(value):
    return old.root_relation(sp.factor(value))


def at_lapse(value):
    return sp.factor(value.subs(model.N,t*t))


def context():
    ctx=constraints.context(((1,0,0),(-1,0,0)),symbolic_scale=True)
    ctx.domain=QQ_I.frac_field(t,P,k,eta,pd)
    ctx.zero_coefficient=ctx.domain.zero
    ctx.one_coefficient=ctx.domain.one
    return ctx


def pair(left,right):
    if type(left) is not str or type(right) is not str or left not in CHANNELS or right not in CHANNELS:
        raise ValueError("Require one of the four declared scalar phase channels")
    return _pair(left,right)


@cache
def _pair(left,right):
    ctx=context()
    fields={name:ctx.jet() for name in CHANNELS}
    for label,name in enumerate((left,right)):
        fields[name]+=ctx.leg(label)
    zero=j.zeros(ctx)
    vectors=[ctx.jet() for _ in range(3)]
    geo=geometry.derive(ctx,fields["curvature"],zero)
    # The canonical background coefficient is -p0/2, not a presumed
    # equality to the actual hat Hubble rate away from the clock.
    p0=eta*pd
    solved=constraints.derive(ctx,geo,fields["curvature"],zero,fields["scalar_p"],zero,
                             fields["matter"],fields["matter_p"],vectors,vectors,-p0/2,P)
    momentum=solved["momentum"]
    volume=geo["volume"]
    mixed=j.mscale(j.mmul(momentum,geo["metric"]),volume.power(-1))
    trace=j.trace(mixed)
    invariant_p=sp.Rational(2,3)*trace
    matter=(P+fields["matter_p"])/volume
    shear=j.trace(j.mmul(mixed,mixed))-trace*trace/3
    gradient=sum(geo["inverse"][a][b]*fields["matter"].derivative(a)*fields["matter"].derivative(b)
                 for a in range(3) for b in range(3))
    d=background.data()
    a,U,f,b=at_lapse(d["trace_a"]),at_lapse(d["volume_U"]),at_lapse(d["scalar_potential"]),eta*at_lapse(d["actual_trace_linear_time_derivative"])
    N=t*t
    on=volume*N*((invariant_p-b)**2/(4*a)-f-shear/(sp.Rational(3,2)*a)
                 +matter*matter/(2*U)+t*gradient/2-geo["curvature"]/(2*t**3))
    # Differentiate at fixed pd and P before imposing their background
    # equations. The full first-order temporal Hamiltonian is used.
    force=j.Jet(ctx,{mask:ctx.convert(sp.diff(ctx.domain.to_sympy(coefficient),t)/(2*t))
                     for mask,coefficient in on.data.items()}).homogeneous(1)
    h=at_lapse(d["fixed_phase_lapse_Hessian"])
    # Two labelled legs give [force^2]_12=2*force_1*force_2.
    # Project the exact first time jet before rational division; omitted
    # time powers cannot enter this derivative. This avoids a needless
    # high-dimensional polynomial GCD in the scalar coefficient domain.
    Hhat=eta*at_lapse(d["hat_Hubble_time_derivative"])
    moving=-2*Hhat*j.contract(momentum,geo["metric"])
    moving+=pd*(6*fields["curvature"]+3*fields["curvature"]**2)/2
    moving-=P*fields["matter_p"]/t
    raw=(on+moving).coefficient(ctx.full)
    left_force,right_force=force.coefficient(1),force.coefficient(2)
    F0,G0=left_force.subs(eta,0),right_force.subs(eta,0)
    F1,G1=sp.diff(left_force,eta).subs(eta,0),sp.diff(right_force,eta).subs(eta,0)
    value0=raw.subs(eta,0)-F0*G0/h
    value1=sp.diff(raw,eta).subs(eta,0)-(F1*G0+F0*G1)/h
    # The actual hat Hubble and canonical trace coefficient are distinct.
    # Their moving boundaries are retained before the principal limit.
    actual_pd=at_lapse(d["trace_momentum_time_derivative"])
    value=clean(value0.subs(pd,actual_pd))+eta*clean(value1.subs(pd,actual_pd))
    return {"kernel_first_time_jet":value,"spatial_checks":solved["checks"],
            "moving_boundary_kernel":sp.factor(moving.coefficient(ctx.full).subs(pd,actual_pd)),
            "fixed_phase_derivative_precedes_background_equations":True}


@cache
def phase_matrix():
    result=sp.zeros(4)
    for a,left in enumerate(CHANNELS):
        for b,right in enumerate(CHANNELS[a:],start=a):
            value=pair(left,right)["kernel_first_time_jet"]
            result[a,b]=value
            result[b,a]=value.subs(k,-k)
    return sp.ImmutableMatrix(result)
