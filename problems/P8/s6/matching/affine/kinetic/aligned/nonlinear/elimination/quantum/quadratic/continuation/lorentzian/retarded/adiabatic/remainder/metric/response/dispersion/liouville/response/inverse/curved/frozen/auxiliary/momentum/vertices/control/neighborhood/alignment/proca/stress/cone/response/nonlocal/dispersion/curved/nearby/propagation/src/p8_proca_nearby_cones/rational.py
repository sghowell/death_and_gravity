"""Primitive-free rational evolution and the actual local clock characteristic."""
import hashlib
import json
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model
from p8_proca_scalar_cone import background
from p8_proca_scalar_cone import principal as central
from sympy.polys.domains import QQ
from sympy.polys.fields import field

u,N,s=model.u,model.N,model.s
z=sp.Symbol("rescaled_trace",real=True)
FIELD,uf,Nf,zf=field((u,N,z),QQ)


@cache
def coefficients():
    c=model.coefficients()
    h=(1+u*u)**3
    D=1+N*N*(h-1)
    e=c["e_omega"]
    Ln=1/(2*N*D)
    Lt=3*u*(1-N*N)/(2*(1+u*u)*D)
    a0=-3*N/4
    m0=D/(2*N*h)
    d=sp.factor(c["b"]+model.I)
    V=sp.factor(N*c["f0"]+model.Iphi)
    IN=sp.factor(-c["boundary_primitive_s"].subs(s,1/N)/N**2)
    d0,V0,I0=(sp.factor(value*e) for value in (d,V,IN))
    b0=sp.factor(sp.diff(d0,N)-Ln*d0-I0)
    A=sp.factor(sp.diff(a0,N)+Ln*a0)
    M=sp.factor(sp.diff(m0,N)+Ln*m0)
    C=sp.factor(-sp.diff(V0,N)+Ln*V0+sp.diff(I0,u)-Lt*I0)
    B=-2*a0*b0
    return {"h":h,"D":D,"e":e,"Ln":Ln,"Lt":Lt,"a0":a0,"m0":m0,
            "d0":d0,"V0":V0,"I0":I0,"b0":b0,"A":A,"M":M,"B":B,"C":C}


@cache
def system():
    c=coefficients()
    names=("Ln","Lt","a0","m0","d0","V0","I0","b0","A","M","B","C","h","D")
    data={name:FIELD.from_expr(c[name]) for name in names}
    Ln,Lt,a0,m0,d0,V0,_I0,b0,A,M,B,C,h,D=(data[name] for name in names)
    u,N,z=uf,Nf,zf
    w2=-(A*z*z+B*z+C)/M
    Kz=Ln*z-b0
    Z0=-a0*z*z+V0-d0.diff(u)+Lt*d0+m0*w2+Lt*z
    pivot=A.diff(N)*z*z+B.diff(N)*z+C.diff(N)+M.diff(N)*w2+(2*A*z+B)*Kz+2*M*Ln*w2
    forcing=A.diff(u)*z*z+B.diff(u)*z+C.diff(u)+M.diff(u)*w2+(2*A*z+B)*Z0+2*M*(Lt-2*a0*z)*w2
    Bbar=-(2*A*z+B)/(3*M)
    H=-N*z/2
    G0=Bbar.diff(u)+Bbar.diff(z)*Z0-(Lt+H)*Bbar-m0*Bbar*Bbar-N*w2
    GN=Bbar.diff(N)+Bbar.diff(z)*Kz-Ln*Bbar
    speed=-4*h*M*M*(G0*pivot-GN*forcing)/(D*pivot*pivot)
    return {**data,"constraint_matter_square":w2,"trace_lapse_chain_coefficient":Kz,
            "trace_flow_at_fixed_lapse":Z0,"on_constraint_rescaled_lapse_Hessian":pivot,
            "constraint_time_forcing":forcing,"rescaled_principal_B":Bbar,
            "hat_Hubble":H,"gradient_fixed_lapse_part":G0,"gradient_lapse_velocity_part":GN,
            "clock_physical_speed_squared":speed,"clock_physical_speed_squared_excess":speed-1}


@cache
def identities():
    d=system()
    return {name:value.as_expr() for name,value in d.items()}


def polynomial_signature(poly):
    rows=[[list(powers),int(value.numerator),int(value.denominator)]
          for powers,value in sorted(poly.items())]
    raw=json.dumps(rows,separators=(",",":"),ensure_ascii=True).encode("ascii")
    return {"monomials":len(rows),"total_degree":max((sum(row[0]) for row in rows),default=0),
            "exact_sorted_monomial_rational_coefficient_sha256":hashlib.sha256(raw).hexdigest()}


@cache
def signatures():
    return {name:{"numerator":polynomial_signature(value.numer),
                  "denominator":polynomial_signature(value.denom)}
            for name,value in system().items()}


@cache
def checks():
    c,d=coefficients(),system()
    old=model.coefficients()
    e=c["e"]
    out={
        "literal_conformal_factor_fourth_power":sp.factor(e**4-N*N*c["h"]/c["D"]),
        "literal_conformal_fixed_lapse_log_derivative":sp.factor(sp.diff(e,N)/e-c["Ln"]),
        "literal_conformal_fixed_time_log_derivative":sp.factor(sp.diff(e,u)/e-c["Lt"]),
        "literal_log_derivative_mixed_partials":sp.factor(sp.diff(c["Ln"],u)-sp.diff(c["Lt"],N)),
        "literal_trace_square_rationalization":sp.factor(N/(4*old["a"])-e*c["a0"]),
        "literal_matter_square_rationalization":sp.factor(N/(2*old["U"])-e*c["m0"]),
        "literal_spatial_curvature_is_minus_matter_square":sp.factor(-N*e*old["B4"]+e*c["m0"]),
        "literal_volume_is_conformal_cube":sp.factor(old["U"]-e**3),
        "actual_shift_removes_only_original_primitive":
            sp.factor(old["b"]+model.I-c["d0"]/e),
        "actual_potential_retains_original_primitive_time_derivative":
            sp.factor(N*old["f0"]+model.Iphi-c["V0"]/e),
        "actual_primitive_endpoint_lapse_derivative":
            sp.factor(-old["boundary_primitive_s"].subs(s,1/N)/N**2-c["I0"]/e),
        "actual_trace_shift_fixed_phase_lapse_derivative":
            sp.factor(sp.diff(c["d0"]/e,N)-c["I0"]/e-c["b0"]/e)}
    for name in ("a0","m0","V0","A","M","C"):
        out["actual_even_time_coefficient_"+name]=sp.factor(c[name].subs(u,-u)-c[name])
    for name in ("d0","I0","b0","B"):
        out["actual_odd_time_coefficient_"+name]=sp.factor(c[name].subs(u,-u)+c[name])
    center=lambda value:value.evaluate([(uf,0),(zf,0)]).as_expr()
    P2=background.data()["matter_momentum_squared"]
    out["actual_central_constraint_matter_square"]=sp.cancel(center(d["constraint_matter_square"])-N*P2)
    out["actual_central_fixed_phase_lapse_Hessian"]=sp.cancel(
        center(d["on_constraint_rescaled_lapse_Hessian"])-sp.sqrt(N)*background.data()["fixed_phase_lapse_Hessian"])
    out["actual_central_constraint_time_forcing"]=center(d["constraint_time_forcing"])
    out["actual_central_clock_speed_reproduces_full_Euler_result"]=sp.cancel(
        center(d["clock_physical_speed_squared"])-central.formula()["clock_physical_speed_squared"])
    speed=d["clock_physical_speed_squared"]
    out["full_time_and_trace_parity_of_clock_speed_numerator"]=sp.Integer(
        sum(1 for powers in speed.numer if (powers[0]+powers[2])%2))
    out["full_time_and_trace_parity_of_clock_speed_denominator"]=sp.Integer(
        sum(1 for powers in speed.denom if (powers[0]+powers[2])%2))
    return out


@cache
def chain_checks():
    # Independent formal algebra, before matter-square elimination.
    ee,zz,ww,ln,lt,aa,mm,bb,ip,du,V,np=sp.symbols(
        "e z w Ln Lt a0 m0 b0 I_phi d_u potential lapse_velocity",nonzero=True)
    y,ell=zz/ee,ww/ee
    a,m=ee*aa,ee*mm
    pd=-a*y*y+V-ip+m*ell*ell
    bd=du-ip+bb*np/ee
    zd=(lt+ln*np)*zz+ee*(pd-bd)
    Z0=-aa*zz*zz+ee*V-ee*du+mm*ww*ww+lt*zz
    wd=(lt+ln*np)*ww-2*a*y*ww
    f,fu,fn,fz,fw,kz,kw,z0,w0,pivot,forcing=sp.symbols(
        "f f_u f_N f_z f_w K_z K_w Z_0 W_0 pivot forcing")
    Ndot=-forcing/pivot
    return {
        "primitive_phi_cancels_between_actual_trace_flow_and_trace_shift":
            sp.factor(zd-Z0-(ln*zz-bb)*np),
        "actual_rescaled_matter_flow_retains_lapse_velocity":
            sp.factor(wd-(lt-2*aa*zz)*ww-ln*ww*np),
        "rescaled_fixed_phase_lapse_Hessian_keeps_off_constraint_term":
            sp.factor((fn+fz*kz+fw*kw-ln*f)/ee-(pivot-ln*f).subs(pivot,fn+fz*kz+fw*kw)/ee),
        "actual_constraint_derivative_vanishes_on_solved_lapse_flow":
            sp.factor((fu+fz*z0+fw*w0+(fn+fz*kz+fw*kw)*Ndot).subs(
                {forcing:fu+fz*z0+fw*w0,pivot:fn+fz*kz+fw*kw},simultaneous=True)),
        "rescaled_hat_Hubble_is_actual_Hamiltonian_trace_derivative":
            sp.factor(2*a*y/3-2*aa*zz/3)}
