"""Full spatial-Hamiltonian quadratic replay in the actual canonical variables."""
from functools import cache

import sympy as sp
from p8_aligned_margin import dynamics as margin
from p8_auxiliary_neighborhood import model
from p8_physical import geometry
from p8_physical import jets as j

from . import constraints

CHANNELS=("curvature","scalar_p","matter","matter_p","tensor","tensor_p","tensor2","tensor2_p",
          "Wx","Wy","Wz","Px","Py","Pz")


@cache
def lapse_coefficients():
    actual=model.coefficients()
    generic=model.generic()
    bg=actual["background"]
    I2=model.adm.clock_lapse()["I_NN"]*(model.N-1)**2/2
    Q2=-3*sp.diff(bg["h"],model.u)*(model.N-1)**2/(2*bg["h"]**3)
    out={}
    for field in (model.dp,model.dc,model.j,model.curvature):
        row=sp.diff(generic["Hamiltonian_over_hat_volume"],field).subs(dict.fromkeys(model.VARIABLES,0))
        row=row.subs(model.specialization(),simultaneous=True).subs({
            model.I:I2,model.Iphi:sp.diff(I2,model.u),model.Q:Q2},simultaneous=True)
        out[str(field)]=sp.factor(sp.diff(row,model.N).subs(model.N,1))
    return out


@cache
def target_scalar():
    old=margin.old
    v,p,chi,pm=sp.symbols("v p chi pm",real=True)
    bg=model.coefficients()["background"]
    wave=constraints.wave_scale
    mapping={old.v:v,old.shift:-(p-3*bg["ell"]*chi)/(2*wave**2),
             old.matter:chi,old.pm:pm-3*bg["ell"]*v,old.pi:0,old.sigma:0,
             old.q:wave**2,old.J:bg["J"],margin.delta_J:4*model.MARGIN/bg["h"]**2}
    mapping.update({symbol:bg[key] for symbol,key in (
        (old.theta,"theta"),(old.lam,"lam"),(old.ell,"ell"),(old.w,"w"),(old.H,"H"))})
    value=margin.scalar()["H"].subs(mapping,simultaneous=True)
    return {"fields":(v,p,chi,pm),"Hamiltonian":sp.factor(value)}


@cache
def compact_scalar():
    target=target_scalar()
    v,p,chi,pm=target["fields"]
    bg=model.coefficients()["background"]
    q=constraints.wave_scale**2
    ell,w=bg["ell"],bg["w"]
    R=-bg["theta"]*p/2+(bg["lam"]*q-sp.Rational(3,2)*ell*w)*v+w*pm/2
    Je=bg["J"]+4*model.MARGIN/bg["h"]**2
    value=(pm-3*ell*v)**2/2-ell*p*chi/2+(q/2+sp.Rational(3,4)*ell**2)*chi**2-q*v**2+R**2/Je
    return {"R":R,"J_e":Je,"Hamiltonian":value,
            "prior_regular_scalar_replay":sp.factor(value-target["Hamiltonian"])}


def fields(ctx,left,right):
    values={name:ctx.jet() for name in ("zeta","scalar_p","chi","chi_p")}
    tensor,tensor_p=j.zeros(ctx),j.zeros(ctx)
    vector,vector_p=([ctx.jet() for _ in range(3)] for _ in range(2))
    for leg,channel in enumerate((left,right)):
        amp=ctx.leg(leg)
        if channel in ("curvature","scalar_p","matter","matter_p"):
            name={"curvature":"zeta","scalar_p":"scalar_p","matter":"chi","matter_p":"chi_p"}[channel]
            values[name]+=amp
        elif channel in ("tensor","tensor_p","tensor2","tensor2_p"):
            E=sp.Matrix([[0,0,0],[0,0,1],[0,1,0]]) if channel.startswith("tensor2") else sp.diag(0,1,-1)
            dest=tensor_p if channel.endswith("_p") else tensor
            coefficient=sp.Rational(1,2) if channel.endswith("_p") else 1
            for a in range(3):
                for b in range(3):
                    dest[a][b]+=amp*E[a,b]*coefficient
        else:
            dest=vector if channel[0]=="W" else vector_p
            dest["xyz".index(channel[1])]+=amp
    return values,tensor,tensor_p,vector,vector_p


def pair(left,right):
    if type(left) is not str or type(right) is not str or left not in CHANNELS or right not in CHANNELS:
        raise ValueError("Require declared physical quadratic phase channels")
    return _pair(left,right)


@cache
def _pair(left,right):
    ctx=constraints.context(((1,0,0),(-1,0,0)),symbolic_time=True,symbolic_scale=True)
    f,t,tp,W,Pi=fields(ctx,left,right)
    geo=geometry.derive(ctx,f["zeta"],t)
    bg=model.coefficients()["background"]
    H,ell=bg["H"],bg["ell"]
    solved=constraints.derive(ctx,geo,f["zeta"],t,f["scalar_p"],tp,f["chi"],f["chi_p"],W,Pi,H,ell)
    momentum=solved["momentum"]
    volume=geo["volume"]
    mixed=j.mscale(j.mmul(momentum,geo["metric"]),volume.power(-1))
    trace=j.trace(mixed)
    dp=sp.Rational(2,3)*trace+2*H
    dc=(ell+f["chi_p"])/volume-ell
    divergence=sum(Pi[a].derivative(a) for a in range(3))/volume
    shear=j.trace(j.mmul(mixed,mixed))-trace*trace/3
    zeta=sp.Rational(1,10**6)
    electric=sum(geo["metric"][a][b]*Pi[a]*Pi[b] for a in range(3) for b in range(3))/zeta/volume**2
    F=[[W[b].derivative(a)-W[a].derivative(b) for b in range(3)] for a in range(3)]
    magnetic=zeta*sum(geo["inverse"][a][c]*geo["inverse"][b][d]*F[a][b]*F[c][d]
                      for a in range(3) for b in range(3) for c in range(3) for d in range(3))
    vector=sum(geo["inverse"][a][b]*W[a]*W[b] for a in range(3) for b in range(3))
    gradient=sum(geo["inverse"][a][b]*f["chi"].derivative(a)*f["chi"].derivative(b)
                 for a in range(3) for b in range(3))
    Fclock=sp.factor(model.adm.coefficients()["F"].subs(model.adm.X,1))
    on=(-sp.Rational(3,4)*(dp-2*H)**2-Fclock+divergence**2/2+2*shear
        +electric/2+magnetic/4+vector/2+(ell+dc)**2/2+gradient/2-geo["curvature"]/2)
    force=(3*bg["theta"]*dp-bg["w"]*dc-bg["lam"]*geo["curvature"]/2).homogeneous(1)
    Je=bg["J"]+4*model.MARGIN/bg["h"]**2
    density=volume*(on+force*force/(4*Je))
    Adot=sp.diff(H,model.u)+3*H*H
    density-=2*H*j.contract(momentum,geo["metric"])
    density-=Adot*(6*f["zeta"]+3*f["zeta"]**2-j.contract(t,t)/2)
    density-=ell*f["chi_p"]
    target=target_scalar()
    poly=sp.Poly(target["Hamiltonian"],*target["fields"])
    expected=ctx.jet()
    variables=(f["zeta"],f["scalar_p"],f["chi"],f["chi_p"])
    for powers,coefficient in poly.terms():
        if any(power and field.is_zero() for field,power in zip(variables,powers)):
            continue
        term=ctx.jet(coefficient)
        for field,power in zip(variables,powers):
            term*=field**power
        expected+=term
    # Tensor E:E=2, momentum is the dual E/2, so this is canonical.
    expected+=2*j.contract(tp,tp)
    expected+=sum(j.contract([[t[a][b].derivative(c) for b in range(3)] for a in range(3)],
                              [[t[a][b].derivative(c) for b in range(3)] for a in range(3)])/8 for c in range(3))
    expected+=sum(Pi[a]*Pi[a] for a in range(3))/(2*zeta)
    expected+=sum(Pi[a].derivative(a) for a in range(3))**2/2
    expected+=zeta*sum(F[a][b]*F[a][b] for a in range(3) for b in range(3))/4
    expected+=sum(W[a]*W[a] for a in range(3))/2
    return {"actual":sp.factor(density.coefficient(ctx.full)),
            "expected":sp.factor(expected.coefficient(ctx.full)),
            "residual":sp.factor((density-expected).coefficient(ctx.full))}


@cache
def checks():
    bg=model.coefficients()["background"]
    coefficients=lapse_coefficients()
    out={"actual_linear_lapse_trace_invariant":sp.factor(coefficients[str(model.dp)]-3*bg["theta"]),
         "actual_linear_lapse_matter_invariant":sp.factor(coefficients[str(model.dc)]+bg["w"]),
         "actual_linear_lapse_vector_divergence_vanishes":coefficients[str(model.j)],
         "actual_linear_lapse_spatial_curvature_invariant":sp.factor(coefficients[str(model.curvature)]+bg["lam"]/2)}
    out["compact_actual_canonical_scalar_Hamiltonian"]=compact_scalar()["prior_regular_scalar_replay"]
    for i,left in enumerate(CHANNELS):
        for right in CHANNELS[i:]:
            out["full_quadratic_phase_pair_"+left+"_"+right]=pair(left,right)["residual"]
    return out


@cache
def phase_matrix():
    k=constraints.wave_scale
    result=sp.zeros(len(CHANNELS))
    for i,left in enumerate(CHANNELS):
        for h,right in enumerate(CHANNELS[i:],start=i):
            value=pair(left,right)["actual"]
            result[i,h]=value
            # The real integrated quadratic Hamiltonian has
            # H_ji(k)=H_ij(-k), without assuming a pointwise symmetric
            # coefficient matrix before integrations by parts.
            if h!=i:
                result[h,i]=value.subs(k,-k)
    return sp.ImmutableMatrix(result)
