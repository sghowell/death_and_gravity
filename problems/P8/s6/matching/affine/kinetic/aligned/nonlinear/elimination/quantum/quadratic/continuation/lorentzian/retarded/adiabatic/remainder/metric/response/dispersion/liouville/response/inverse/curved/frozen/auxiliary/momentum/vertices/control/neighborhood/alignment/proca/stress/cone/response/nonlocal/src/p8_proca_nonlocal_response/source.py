"""New constant-physical-mass prepared two-source mode and readout vertices."""
from functools import cache

import sympy as sp
from p8_vector_hadamard import series
from p8_vector_metric_local import canonical, variation
from p8_vector_metric_local import jets as generic
from p8_vector_state import wkb

u,z=wkb.u,wkb.z
n=sp.symbols("ordinary_prepared_lapse_source0:12",real=True)
v=sp.symbols("ordinary_prepared_logscale_source0:12",real=True)
H,lam=wkb.background()["H"],wkb.background()["lambda"]
kind=generic.kind
ZERO={field:sp.Integer(0) for row in (generic.alpha,generic.beta,generic.alpha2,generic.beta2) for field in row}


def output(value):
    if type(value) is not str or value not in ("energy","pressure"):
        raise ValueError("Require physical energy or pressure readout")
    return value


def clean(value):
    fields=n+v
    polynomial=sp.Poly(sp.expand(value),*fields)
    return sp.Add(*(sp.factor(coefficient)*sp.prod(field**power for field,power in zip(fields,powers))
                    for powers,coefficient in polynomial.terms()))


def time(value):
    return clean(sp.diff(value,u)+wkb.background()["z_prime"]*sp.diff(value,z)
                 +sum(sp.diff(value,row[j])*row[j+1] for row in (n,v) for j in range(len(row)-1)))


def project(value):
    # Eliminate the actual new mass jets before the generic clock substitution.
    actual=generic.actual(value.xreplace(ZERO).subs(generic.D,3))
    mapping={generic.z:z}
    mapping.update({field:n[j] for j,field in enumerate(generic.n)})
    mapping.update({field:v[j] for j,field in enumerate(generic.v)})
    return clean(actual.subs(mapping,simultaneous=True))


def data(sector):
    return _data(kind(sector))


@cache
def _data(sector):
    physical=canonical.data(sector)
    d=H/2 if sector=="T" else H*(sp.Rational(1,2)+z)
    r,rg=project(physical["delta_log_frequency"]),project(physical["delta_log_g"])
    dd,dl=time(rg),time(r)
    return {"rate":d,"r":r,"rg":rg,"delta_d":dd,"delta_lambda":dl,
            "delta_U":time(dd)+2*d*dd,"c1":d+lam/2,"delta_c1":dd+dl/2}


def readout(sector,component):
    return _readout(kind(sector),output(component))


@cache
def _readout(sector,component):
    physical=canonical.data(sector)
    label,factor=("N",-sp.Integer(1)) if component=="energy" else ("Z",sp.Rational(1,3))
    A,B=[factor*project(value) for value in physical["weights"][label]]
    dA,dB=[factor*project(value) for value in physical["varied_weights"][label]]
    normalization=3*v[0] if component=="energy" else n[0]+3*v[0]
    return {"A":A,"B":B,"delta_A":clean(dA-normalization*A),
            "delta_B":clean(dB-normalization*B)}


def adiabatic(sector,component,order):
    sector,component,order=kind(sector),output(component),generic.order(order)
    return _adiabatic(sector,component,order)


@cache
def _adiabatic(sector,component,order):
    label,factor=("N",-sp.Integer(1)) if component=="energy" else ("Z",sp.Rational(1,3))
    generic_rows=variation.data(sector)
    varied=factor*project(generic_rows["coefficients"][label][order])
    baseline=factor*project(generic_rows["baseline"][label][order])
    normalization=3*v[0] if component=="energy" else n[0]+3*v[0]
    return clean(varied-normalization*baseline)


def baseline(sector,order):
    kind(sector)
    if type(order) is not int or not 0<=order<=4:
        raise ValueError("Require native reference coefficient order 0..4")
    return series.coefficient("transverse" if sector=="T" else "longitudinal",order)


def linear_bound(value):
    if not isinstance(value,sp.Expr) or value.has(sp.Float,sp.oo,-sp.oo,sp.zoo,sp.nan):
        raise ValueError("Require an exact finite linear two-source expression")
    fields=n+v
    polynomial=sp.Poly(sp.expand(value),*fields)
    total,reconstructions,max_order=sp.Integer(0),[],0
    for powers,coefficient in polynomial.terms():
        if coefficient==0:
            continue
        if sum(powers)!=1:
            raise ValueError("Require a linear prepared metric-source coefficient")
        derivative=powers.index(1)%len(n)
        if derivative>10:
            raise ValueError("The C10 source norm does not control this derivative")
        if coefficient.free_symbols-{u,z}:
            raise ValueError("Only fixed-clock time and momentum fraction may enter source coefficients")
        bound=wkb.box_bound(sp.factor(coefficient))
        total+=bound["absolute_upper"]
        reconstructions.append(bound["reconstruction"])
        max_order=max(max_order,derivative)
    return {"absolute_upper":total,"reconstructions":reconstructions,"highest_source_derivative":max_order}


@cache
def checks():
    from p8_proca_stress import subtraction
    N,a,q,m=sp.symbols("new_physical_lapse new_physical_scale new_physical_q new_mass",positive=True)
    W2=sp.Symbol("new_mode_frequency_squared",positive=True)
    out={}
    scale=lambda value:a*sp.diff(value,a)-2*q*sp.diff(value,q)
    change=lambda value:n[0]*sp.diff(value,N)+v[0]*scale(value)
    on=lambda value:sp.factor(value.subs(N,1).subs(q,m*m*z/(1-z)))
    Dt=lambda value:time(value)+2*lam*W2*sp.diff(value,W2)
    for sector in ("T","L"):
        g2=a/N if sector=="T" else a**3*m*m*q/(N*(q+m*m))
        frequency2=N*N*(q+m*m)
        d=data(sector)
        out[sector+"_literal_minimal_frequency_variation"]=clean(on(change(frequency2)/(2*frequency2))-d["r"])
        out[sector+"_literal_minimal_canonical_map_variation"]=clean(on(change(g2)/(2*g2))-d["rg"])
        rows,delta={},{}
        for component in ("energy","pressure"):
            derivative=(lambda value:sp.diff(value,N)) if component=="energy" else scale
            factor=-sp.Integer(1) if component=="energy" else sp.Rational(1,3)
            raw=(derivative(g2)/g2,-derivative(g2*frequency2)/(g2*frequency2))
            normalization=3*v[0] if component=="energy" else n[0]+3*v[0]
            target=readout(sector,component)
            for key,value in zip(("A","B"),raw):
                out[sector+"_"+component+"_literal_weight_"+key]=clean(factor*on(value)-target[key])
                out[sector+"_"+component+"_literal_contact_"+key]=clean(
                    factor*on(change(value))-normalization*target[key]-target["delta_"+key])
            rows[component]=sp.Matrix([[target["B"]*W2,0,target["A"]]])
            delta[component]=sp.Matrix([[(target["delta_B"]+2*d["r"]*target["B"])*W2,0,target["delta_A"]]])
        rate,dd,rf=d["rate"],d["delta_d"],d["r"]
        M=sp.Matrix([[2*rate,2,0],[-W2,0,1],[0,-2*W2,-2*rate]])
        dM=sp.Matrix([[2*dd,0,0],[-2*rf*W2,0,0],[0,-4*rf*W2,-2*dd]])
        E,P=rows["energy"],rows["pressure"]
        dE,dP=delta["energy"],delta["pressure"]
        base=E.applyfunc(Dt)-3*H*E+E*M+3*H*(E+P)
        varied=dE.applyfunc(Dt)-3*H*dE+dE*M+E*dM+3*H*(dE+dP)+3*v[1]*(E+P)
        out[sector+"_ordinary_mode_covariance_Ward"]=sp.ImmutableMatrix(base.applyfunc(clean))
        out[sector+"_actual_varied_mode_covariance_Ward"]=sp.ImmutableMatrix(varied.applyfunc(clean))
        for order in range(3):
            baseline_new=-project(variation.data(sector)["baseline"]["N"][order])
            prior=subtraction.adiabatic_rows("transverse" if sector=="T" else "longitudinal",0)[order]
            out[sector+"_actual_new_energy_adiabatic_baseline_"+str(order)]=clean(baseline_new-prior)
            pressure_baseline=project(variation.data(sector)["baseline"]["Z"][order])/3
            varied_energy=adiabatic(sector,"energy",order)
            varied_pressure=adiabatic(sector,"pressure",order)
            out[sector+"_actual_varied_adiabatic_Ward_"+str(order)]=clean(
                time(varied_energy)+(1-2*order)*lam*varied_energy+3*H*varied_pressure
                +3*v[1]*(baseline_new+pressure_baseline))
    return out
