"""Exact central lapse acceleration and physical bounce acceleration, before family substitution."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model
from p8_proca_scalar_cone import background, domain

u,N,s=model.u,model.N,model.s
ell=sp.Symbol("fixed_nearby_matter_density",real=True)


@cache
def data():
    c,bg=model.coefficients(),background.data()
    first=sp.factor(sp.diff(c["boundary_primitive_s"],u).subs(u,0))
    third=sp.factor(sp.diff(c["boundary_primitive_s"],u,3).subs(u,0))
    I1=sp.factor(sp.integrate(first,(s,1,1/N)))
    I3=sp.factor(sp.integrate(third,(s,1,1/N)))
    a,U=bg["trace_a"],bg["volume_U"]
    bu=bg["actual_trace_linear_time_derivative"]
    fuu=sp.factor(sp.diff(c["f0"].subs(model.Iphi,0),u,2).subs(u,0)-I3/N)
    invUuu=sp.factor(sp.diff(1/c["U"],u,2).subs(u,0))
    Huu=N*(bu**2/(2*a)-fuu+ell*ell*invUuu/2)
    pd,hd=bg["trace_momentum_time_derivative"],bg["hat_Hubble_time_derivative"]
    ldd=-3*ell*hd
    numerator=sp.diff(Huu,N)+2*sp.diff(-N*bu/(2*a),N)*pd+sp.diff(N/(2*a),N)*pd**2+sp.diff(N*ell/U,N)*ldd
    Ndd=sp.factor(-numerator.subs(ell*ell,bg["matter_momentum_squared"])/bg["fixed_phase_lapse_Hessian"])
    Hphysical=sp.factor((hd+sp.Rational(3,2)*(1-N*N)+Ndd/(2*N))/N**2)
    return {"primitive_first_time_integrand":first,"primitive_third_time_integrand":third,
            "actual_boundary_primitive_phi":I1,"actual_boundary_primitive_phi_phi_phi":I3,
            "actual_scalar_potential_second_time_derivative":fuu,
            "inverse_volume_second_time_derivative":invUuu,
            "Hamiltonian_second_time_derivative_before_family":Huu,
            "fixed_phase_constraint_second_path_derivative_numerator":numerator,
            "actual_trace_momentum_first_time_derivative":pd,
            "actual_matter_density_second_time_derivative":ldd,
            "actual_hat_Hubble_first_time_derivative":hd,
            "actual_lapse_second_time_derivative":Ndd,
            "missing_third_primitive_jet_lapse_acceleration_error":sp.factor(sp.diff(I3,N)/bg["fixed_phase_lapse_Hessian"]),
            "actual_physical_Hubble_proper_time_derivative":Hphysical,
            "actual_physical_logscale_second_clock_time_derivative":sp.factor(N*N*Hphysical),
            "physical_Hubble_at_center":sp.Integer(0),
            "physical_scale_at_center":sp.sqrt(N),
            "fixed_phase_lapse_Hessian":bg["fixed_phase_lapse_Hessian"],
            "positive_matter_density_squared":bg["matter_momentum_squared"],
            "quantum_profile_and_response_included":False}


@cache
def checks():
    d=data()
    first=d["actual_boundary_primitive_phi"].subs(N,1/s)
    third=d["actual_boundary_primitive_phi_phi_phi"].subs(N,1/s)
    out={"old_actual_boundary_first_time_derivative_replayed":sp.factor(d["actual_boundary_primitive_phi"]-background.data()["actual_boundary_primitive_phi"]),
         "new_third_boundary_integrand_antiderivative":sp.factor(sp.diff(third,s)-d["primitive_third_time_integrand"]),
         "new_third_boundary_fixed_basepoint":third.subs(s,1),
         "old_first_boundary_integrand_antiderivative":sp.factor(sp.diff(first,s)-d["primitive_first_time_integrand"]),
         "original_clock_lapse_acceleration_zero":d["actual_lapse_second_time_derivative"].subs(N,1),
         "original_clock_physical_Hubble_acceleration":d["actual_physical_Hubble_proper_time_derivative"].subs(N,1)-4}
    a,U=background.data()["trace_a"],background.data()["volume_U"]
    bu=background.data()["actual_trace_linear_time_derivative"]
    pd=d["actual_trace_momentum_first_time_derivative"]
    hd=d["actual_hat_Hubble_first_time_derivative"]
    p,t=sp.symbols("independent_fixed_phase_trace independent_clock_path_parameter",real=True)
    # Expand a literal fixed-phase Hamiltonian through path degree two before imposing its path.
    H=N*((p-bu*t)**2/(4*a)-background.data()["scalar_potential"]
         -d["actual_scalar_potential_second_time_derivative"]*t*t/2
         +ell*ell/(2*U)+ell*ell*d["inverse_volume_second_time_derivative"]*t*t/4)
    F=sp.diff(H,N)
    literal=sp.expand(F.subs({p:pd*t,ell:ell-sp.Rational(3,2)*ell*hd*t*t},simultaneous=True)).coeff(t,2)
    out["independent_fixed_phase_constraint_path_Taylor"]=sp.factor(
        (2*literal-d["fixed_phase_constraint_second_path_derivative_numerator"]).subs(ell*ell,d["positive_matter_density_squared"]))
    out["actual_lapse_acceleration_satisfies_twice_differentiated_constraint"]=sp.factor(
        d["fixed_phase_lapse_Hessian"]*d["actual_lapse_second_time_derivative"]
        +d["fixed_phase_constraint_second_path_derivative_numerator"].subs(ell*ell,d["positive_matter_density_squared"]))
    h=(1+u*u)**3
    omega=-sp.log((h-1+N**-2)/h)/4
    out["actual_physical_logscale_fixed_lapse_second_time_jet"]=sp.factor(
        sp.diff(omega,u,2).subs(u,0)-sp.Rational(3,2)*(1-N*N))
    out["actual_physical_logscale_lapse_jet_at_center"]=sp.factor(sp.diff(omega,N).subs(u,0)-1/(2*N))
    out["physical_proper_Hubble_acceleration_keeps_lapse_acceleration"]=sp.factor(
        N*N*d["actual_physical_Hubble_proper_time_derivative"]-hd-sp.Rational(3,2)*(1-N*N)-d["actual_lapse_second_time_derivative"]/(2*N))
    out["physical_acceleration_is_not_just_hat_acceleration"]=sp.factor(
        d["actual_physical_Hubble_proper_time_derivative"]-(hd+sp.Rational(3,2)*(1-N*N)+d["actual_lapse_second_time_derivative"]/(2*N))/N**2)
    wrong_force_numerator=d["fixed_phase_constraint_second_path_derivative_numerator"]-sp.diff(d["actual_boundary_primitive_phi_phi_phi"],N)
    wrong_acceleration=sp.factor(-wrong_force_numerator.subs(ell*ell,d["positive_matter_density_squared"])/d["fixed_phase_lapse_Hessian"])
    out["missing_third_primitive_jet_negative_control"]=sp.factor(wrong_acceleration-d["actual_lapse_second_time_derivative"]-d["missing_third_primitive_jet_lapse_acceleration_error"])
    return out


@cache
def boxes():
    d=data()
    values={"matter_density_squared":d["positive_matter_density_squared"],
            "lapse_Hessian_times_N_seven_halves":sp.factor(d["fixed_phase_lapse_Hessian"]*N**sp.Rational(7,2)),
            "lapse_acceleration":d["actual_lapse_second_time_derivative"],
            "physical_Hubble_proper_acceleration":d["actual_physical_Hubble_proper_time_derivative"]}
    return {key:domain.enclosure(value,sp.Rational(1,10**6)) for key,value in values.items()}
