"""Actual homogeneous central time jets before fixed-phase lapse reduction."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import center, model
from p8_offclock_scalar import datum

u,N,P=model.u,model.N,datum.P
pd=sp.Symbol("fixed_canonical_trace_time_derivative",real=True)


@cache
def data():
    c=model.coefficients()
    Iphi=center.data()["actual_boundary_primitive_phi_at_bounce"].subs(center.r,1/N)
    bu=sp.factor(sp.diff(c["b"],u).subs(u,0)-Iphi)
    a=-1/(3*sp.sqrt(N))
    U=N**sp.Rational(3,2)
    f0=sp.factor(c["f0"].subs({u:0,model.I:0,model.Iphi:Iphi,model.Q:0}))
    P2=datum.data()["matter_momentum_squared"]
    pdot=sp.factor(N*(f0+P**2/(2*U)))
    on=sp.factor(pdot.subs(P**2,P2))
    hdot=sp.factor(N*(on-bu)/(6*a))
    derivative_before_family=sp.diff(N/(2*a),N)*pd-sp.diff(N*bu/(2*a),N)
    Fdot=sp.factor(derivative_before_family.subs(pd,on))
    return {"trace_a":a,"volume_U":U,"scalar_potential":f0,
            "actual_boundary_primitive_phi":Iphi,"actual_trace_linear_time_derivative":bu,
            "matter_momentum_squared":P2,"fixed_phase_lapse_Hessian":datum.data()["fixed_phase_lapse_Hessian"],
            "trace_momentum_time_derivative_before_family":pdot,
            "trace_momentum_time_derivative":on,
            "hat_Hubble_time_derivative":hdot,
            "lapse_time_derivative":sp.Integer(0),
            "matter_density_time_derivative":sp.Integer(0),
            "trace_force_time_derivative_before_family":derivative_before_family,
            "trace_force_time_derivative":Fdot}


@cache
def checks():
    d=data()
    c=model.coefficients()
    R,pr,l,N0,a,b,f,U=sp.symbols("hat_scale trace_momentum matter_density lapse a b f U",nonzero=True)
    canonical=sp.Symbol("canonical_scale_momentum",real=True)
    H=R**3*N0*((canonical/(3*R**2)-b)**2/(4*a)-f+l*l/(2*U*R**6))
    Rdot=sp.diff(H,canonical).subs(canonical,3*R**2*pr)
    pdot=-sp.diff(H,R).subs(canonical,3*R**2*pr)/(3*R**2)-2*Rdot*pr/R
    out={"actual_homogeneous_scale_Hamilton_equation":sp.factor(Rdot/R-N0*(pr-b)/(6*a)),
         "actual_homogeneous_trace_Hamilton_equation":sp.factor(
             pdot-N0*(-(pr-b)**2/(4*a)+f+l*l/(2*U*R**6))),
         "actual_center_trace_coefficient":sp.factor(c["a"].subs(u,0)-d["trace_a"]),
         "actual_center_volume_coefficient":sp.factor(c["U"].subs(u,0)-d["volume_U"]),
         "zero_center_trace_linear_coefficient":sp.factor(c["b"].subs({u:0,model.I:0})),
         "actual_bounce_trace_momentum_derivative":d["trace_momentum_time_derivative"].subs(N,1)+8,
         "actual_bounce_hat_Hubble_derivative":d["hat_Hubble_time_derivative"].subs(N,1)-4,
         "actual_bounce_trace_force_derivative":d["trace_force_time_derivative"].subs(N,1)-9}
    # Every even coefficient has zero first clock derivative at the slice.
    for name in ("U","a","e_omega","B4"):
        out["zero_center_first_time_derivative_"+name]=sp.factor(sp.diff(c[name],u).subs(u,0))
        out["actual_even_background_coefficient_"+name]=sp.factor(c[name].subs(u,-u)-c[name])
    out["actual_odd_boundary_primitive_integrand"]=sp.factor(
        c["boundary_primitive_s"].subs(u,-u)+c["boundary_primitive_s"])
    out["actual_odd_trace_linear_coefficient"]=sp.factor(
        c["b"].subs({u:-u,model.I:-model.I},simultaneous=True)+c["b"])
    out["actual_even_scalar_potential_with_even_Iphi"]=sp.factor(c["f0"].subs(u,-u)-c["f0"])
    target=model.adm.coefficients()
    X=model.adm.X
    out["global_exceptional_luminal_matter_relation"]=sp.factor(
        target["f"].diff(X)*2/X-model.adm.dictionary.target()["alpha3"].subs(model.adm.dictionary.x,-X))
    return out
