"""Full geometric scalar reduction on an arbitrary homogeneous solution."""
from functools import cache

import sympy as sp
from p8_coupled_momentum import constraints
from p8_physical import geometry
from p8_physical import jets as j
from sympy.polys.domains import QQ_I

a,c,f,m,g,r,an,cn,mn,rn,h,p0,ell=sp.symbols("trace_square trace_linear potential matter_square matter_gradient curvature trace_square_N trace_linear_N matter_square_N curvature_N lapse_Hessian background_trace background_matter",nonzero=True)
k=constraints.wave_scale
CHANNELS=("v","p","chi","P")
v,p,chi,P=sp.symbols("v p chi P",real=True)
H=(2*a*p0+c)/3
pd=-a*p0*p0-c*p0-f+m*ell*ell
alpha=(2*an*p0+cn)/3
beta=2*mn*ell


@cache
def target():
    F=alpha*p+beta*(P-3*ell*v)+4*rn*k*k*v
    base=m*(P-3*ell*v)**2+sp.Rational(2,3)*a*ell*p*chi+(g*k*k-a*ell*ell)*chi*chi+2*r*k*k*v*v
    return sp.factor(base-F*F/(2*h))


@cache
def _pair(left,right):
    ctx=constraints.context(((1,0,0),(-1,0,0)),symbolic_scale=True)
    ctx.domain=QQ_I.frac_field(a,c,f,m,g,r,an,cn,mn,rn,h,p0,ell,k)
    ctx.zero_coefficient,ctx.one_coefficient=ctx.domain.zero,ctx.domain.one
    fields={name:ctx.jet() for name in CHANNELS}
    for i,name in enumerate((left,right)):
        fields[name]+=ctx.leg(i)
    zero=j.zeros(ctx)
    W=[ctx.jet() for _ in range(3)]
    geo=geometry.derive(ctx,fields["v"],zero)
    solved=constraints.derive(ctx,geo,fields["v"],zero,fields["p"],zero,
                             fields["chi"],fields["P"],W,W,-p0/2,ell)
    mom=solved["momentum"]
    vol=geo["volume"]
    mixed=j.mscale(j.mmul(mom,geo["metric"]),vol.power(-1))
    tr=j.trace(mixed)
    pinv=sp.Rational(2,3)*tr
    linv=(ell+fields["P"])/vol
    shear=j.trace(j.mmul(mixed,mixed))-tr*tr/3
    grad=sum(geo["inverse"][i][z]*fields["chi"].derivative(i)*fields["chi"].derivative(z)
             for i in range(3) for z in range(3))
    on=vol*(a*pinv*pinv+c*pinv+f+m*linv*linv-sp.Rational(8,3)*a*shear+g*grad+r*geo["curvature"])
    fn=-an*p0*p0-cn*p0-mn*ell*ell
    force=(vol*(an*pinv*pinv+cn*pinv+fn+mn*linv*linv-sp.Rational(8,3)*an*shear
                +rn*geo["curvature"])).homogeneous(1)
    moving=-2*H*j.contract(mom,geo["metric"])
    moving+=(pd+3*H*p0)*(6*fields["v"]+3*fields["v"]**2)/2
    moving-=2*m*ell*fields["P"]
    actual=(on+moving-force*force/(2*h)).coefficient(ctx.full)
    variables={"v":v,"p":p,"chi":chi,"P":P}
    expected=sp.diff(target(),variables[left],variables[right])
    return {"actual":sp.factor(actual),"expected":sp.factor(expected),"residual":sp.factor(actual-expected)}


def pair(left,right):
    if not isinstance(left,str) or not isinstance(right,str):
        raise TypeError("Scalar phase channels must be strings")
    if left not in CHANNELS or right not in CHANNELS:
        raise ValueError("Unknown scalar phase channel")
    return _pair(left,right)


@cache
def checks():
    return {"full_moving_geometric_phase_pair_"+left+"_"+right:pair(left,right)["residual"]
            for i,left in enumerate(CHANNELS) for right in CHANNELS[i:]}


def controls():
    calls=[(value,"v") for value in (None,1,True,sp.Symbol("v"),[],sp.nan)]
    calls += [("v",value) for value in ("","tensor","unknown",None,0,False,[])]
    rejected=0
    for left,right in calls:
        try:
            pair(left,right)
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(calls):
        raise ValueError("An invalid geometric phase pair was accepted")
    return {"rejected_inputs":rejected}
