"""Exact flat-York recursion with both physical matter and vector generators."""
from fractions import Fraction

import sympy as sp
from p8_auxiliary_neighborhood import model
from p8_physical import jets as j
from p8_physical import momentum as old
from sympy.polys.domains import QQ_I

ZeroMomentumConstraint=old.ZeroMomentumConstraint
wave_scale=sp.Symbol("positive_wave_scale",positive=True)


def exact(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational,Fraction)):
        raise TypeError("Require an exact rational wavevector component")
    if isinstance(value,Fraction):
        return sp.Rational(value.numerator,value.denominator)
    return sp.Rational(value)


def context(momenta,*,symbolic_time=False,symbolic_scale=False):
    if type(symbolic_time) is not bool or type(symbolic_scale) is not bool:
        raise TypeError("Require native boolean parameter flags")
    if type(momenta) not in (tuple,list) or len(momenta) not in (2,3,4):
        raise ValueError("Require two through four labelled external wavevectors")
    if any(type(v) not in (tuple,list) or len(v)!=3 for v in momenta):
        raise ValueError("Require three components in every wavevector")
    vectors=tuple(tuple(exact(x) for x in v) for v in momenta)
    if any(sum(v[a] for v in vectors)!=0 for a in range(3)):
        raise ValueError("Require zero total external wavevector")
    n=len(vectors)
    for mask in range(1,(1<<n)-1):
        wave=tuple(sum(vectors[i][a] for i in range(n) if mask&(1<<i)) for a in range(3))
        if sum(x*x for x in wave)==0:
            raise ZeroMomentumConstraint("A proper-subset wavevector vanishes")
    result=j.Context(vectors)
    parameters=([model.u] if symbolic_time else [])+([wave_scale] if symbolic_scale else [])
    if parameters:
        result.domain=QQ_I.frac_field(*parameters)
        result.zero_coefficient=result.domain.zero
        result.one_coefficient=result.domain.one
    if symbolic_scale:
        result.wave={mask:tuple(wave_scale*x for x in wave) for mask,wave in result.wave.items()}
    return result


def physical_current(chi,chi_p,vector,vector_p,background_matter_density):
    divergence=sum(vector_p[a].derivative(a) for a in range(3))
    current=[(background_matter_density+chi_p)*chi.derivative(a)
             +sum(vector_p[b]*(vector[b].derivative(a)-vector[a].derivative(b)) for b in range(3))
             -vector[a]*divergence for a in range(3)]
    return current


def residual(momentum,geo,current):
    gravity=old.residual(momentum,geo["christoffel"])
    return [gravity[a]-sum(geo["inverse"][a][b]*current[b] for b in range(3))/2 for a in range(3)]


def derive(ctx,geo,zeta,tensor,scalar_p,tensor_p,chi,chi_p,vector,vector_p,H,ell,*,order=None):
    order=ctx.n-1 if order is None else order
    if type(order) is not int or not 1<=order<=min(3,ctx.n-1):
        raise ValueError("Require native solved order between one and the proper-subset degree")
    unit=j.identity(ctx)
    free=j.madd(j.mscale(unit,-H*(1+zeta)+scalar_p/6),
                j.madd(tensor_p,j.mscale(tensor,H)))
    total=free
    current=physical_current(chi,chi_p,vector,vector_p,ell)
    vector_orders=[]
    for degree in range(1,order+1):
        source=[value.homogeneous(degree) for value in residual(total,geo,current)]
        correction=old.solve_vector(ctx,source)
        total=j.madd(total,old.york(correction))
        vector_orders.append(correction)
    remaining=residual(total,geo,current)
    checks={f"constraint_{degree}_{a}":remaining[a].homogeneous(degree).is_zero()
            for degree in range(1,order+1) for a in range(3)}
    if not all(checks.values()):
        raise ValueError("The full coupled spatial momentum recursion failed")
    return {"momentum":total,"vector_orders":vector_orders,"physical_lower_current":current,"checks":checks}


def linear_target(ctx,scalar_p,chi,ell):
    source=scalar_p-3*ell*chi
    result=[]
    for a in range(3):
        values={}
        for mask,value in source.data.items():
            wave=ctx.wave[mask]
            k2=sum(x*x for x in wave)
            if k2==0:
                if value!=ctx.zero_coefficient:
                    raise ZeroMomentumConstraint("Nonzero homogeneous linear momentum source")
                continue
            values[mask]=value*ctx.convert(sp.I*wave[a]/(8*k2))
        result.append(j.Jet(ctx,values))
    return result
