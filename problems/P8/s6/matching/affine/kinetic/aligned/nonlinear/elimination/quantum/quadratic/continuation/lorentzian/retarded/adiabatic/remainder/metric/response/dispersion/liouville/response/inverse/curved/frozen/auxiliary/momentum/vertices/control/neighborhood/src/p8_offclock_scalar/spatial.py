"""Full spatially reduced central-slice scalar/longitudinal phase Hessian."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import center, model
from p8_coupled_momentum import constraints
from p8_physical import geometry
from p8_physical import jets as j
from sympy.polys.domains import QQ_I

from . import datum

t=sp.Symbol("positive_lapse_square_root",positive=True)
k=constraints.wave_scale
P=datum.P
CHANNELS=("curvature","scalar_p","matter","matter_p","longitudinal","longitudinal_p")


def at_lapse(value):
    return sp.factor(value.subs(model.N,t*t))


def root_relation(value):
    """Reduce polynomial matter powers only after fixed-phase derivatives."""
    value=sp.cancel(value)
    numerator,denominator=sp.fraction(value)
    relation=P**2-at_lapse(datum.data()["matter_momentum_squared"])
    a=sp.rem(numerator,relation,P)
    b=sp.rem(denominator,relation,P)
    if b==0:
        raise ValueError("An on-family denominator vanished")
    return sp.factor(a/b)


def context():
    ctx=constraints.context(((1,0,0),(-1,0,0)),symbolic_scale=True)
    ctx.domain=QQ_I.frac_field(t,P,k)
    ctx.zero_coefficient=ctx.domain.zero
    ctx.one_coefficient=ctx.domain.one
    return ctx


def fields(ctx,left,right):
    f={name:ctx.jet() for name in CHANNELS}
    for index,name in enumerate((left,right)):
        f[name]+=ctx.leg(index)
    vector=[f["longitudinal"].derivative(a) for a in range(3)]
    vector_p=[]
    for a in range(3):
        value=ctx.jet()
        for mask,coefficient in f["longitudinal_p"].data.items():
            wave=ctx.wave[mask]
            q=sum(x*x for x in wave)
            value+=j.Jet(ctx,{mask:coefficient*ctx.convert(sp.I*wave[a]/q)})
        vector_p.append(value)
    return f,vector,vector_p


def invariant_polynomial(ctx,expression,values):
    result=ctx.jet()
    for powers,coefficient in sp.Poly(sp.expand(expression),*model.VARIABLES).terms():
        term=ctx.jet(at_lapse(coefficient))
        for value,power in zip(values,powers):
            if power:
                term*=value**power
        result+=term
    return result


def derivative_N(value):
    ctx=value.context
    return j.Jet(ctx,{mask:ctx.convert(root_relation(sp.diff(ctx.domain.to_sympy(coefficient),t)/(2*t)))
                     for mask,coefficient in value.data.items()})


def pair(left,right):
    if type(left) is not str or type(right) is not str or left not in CHANNELS or right not in CHANNELS:
        raise ValueError("Require declared coupled scalar phase channels")
    return _pair(left,right)


@cache
def _pair(left,right):
    ctx=context()
    f,W,Pi=fields(ctx,left,right)
    zero=j.zeros(ctx)
    geo=geometry.derive(ctx,f["curvature"],zero)
    solved=constraints.derive(ctx,geo,f["curvature"],zero,f["scalar_p"],zero,
                              f["matter"],f["matter_p"],W,Pi,0,P)
    momentum=solved["momentum"]
    volume=geo["volume"]
    mixed=j.mscale(j.mmul(momentum,geo["metric"]),volume.power(-1))
    trace=j.trace(mixed)
    dp=sp.Rational(2,3)*trace
    dc=(P+f["matter_p"])/volume-sp.Rational(1,10)
    divergence=sum(Pi[a].derivative(a) for a in range(3))/volume
    shear=j.trace(j.mmul(mixed,mixed))-trace*trace/3
    zeta=sp.Rational(1,10**6)
    electric=sum(geo["metric"][a][b]*Pi[a]*Pi[b] for a in range(3) for b in range(3))/zeta/volume**2
    F=[[W[b].derivative(a)-W[a].derivative(b) for b in range(3)] for a in range(3)]
    magnetic=zeta*sum(geo["inverse"][a][c]*geo["inverse"][b][d]*F[a][b]*F[c][d]
                      for a in range(3) for b in range(3) for c in range(3) for d in range(3))
    vector=sum(geo["inverse"][a][b]*W[a]*W[b] for a in range(3) for b in range(3))
    gradient=sum(geo["inverse"][a][b]*f["matter"].derivative(a)*f["matter"].derivative(b)
                 for a in range(3) for b in range(3))
    values=(dp,dc,divergence,shear,electric,magnetic,vector,gradient,geo["curvature"])
    on=volume*invariant_polynomial(ctx,center.data()["full_closed_bounce_Hamiltonian"],values)
    force=derivative_N(on).homogeneous(1)
    hessian=at_lapse(datum.data()["fixed_phase_lapse_Hessian"])
    reduced=on-force*force/(2*hessian)
    return {"kernel":root_relation(reduced.coefficient(ctx.full)),
            "unreduced_kernel":root_relation(on.coefficient(ctx.full)),
            "lapse_force_left":root_relation(force.coefficient(1)),
            "lapse_force_right":root_relation(force.coefficient(2)),
            "spatial_checks":solved["checks"],
            "matter_density_relation_imposed_after_N_derivative":True}


@cache
def phase_matrix():
    result=sp.zeros(6)
    for a,left in enumerate(CHANNELS):
        for b,right in enumerate(CHANNELS[a:],start=a):
            value=pair(left,right)["kernel"]
            result[a,b]=value
            result[b,a]=value.subs(k,-k)
    return sp.ImmutableMatrix(result)
