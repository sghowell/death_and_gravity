"""Literal covariant margin increment and freshly derived rational constraint flow."""
from functools import cache

import sympy as sp
from p8_aligned_margin import action as original_margin
from p8_proca_nearby_cones import rational as old

OLD_MARGIN=sp.Rational(1,10**6)
NEW_MARGIN=sp.Rational(1,200)
INCREMENT=NEW_MARGIN-OLD_MARGIN
u,N,z=old.u,old.N,old.z
FIELD,uf,Nf,zf=old.FIELD,old.uf,old.Nf,old.zf


@cache
def action():
    c=old.coefficients()
    X=sp.Symbol("positive_physical_clock_X",positive=True)
    lower=INCREMENT*(X-1)**2/c["h"]**2
    density=N*c["e"]**3*lower.subs(X,N**-2)
    deltaV=sp.factor(c["e"]*density)
    return {"physical_clock_X":X,"old_total_margin":OLD_MARGIN,
        "new_total_margin":NEW_MARGIN,"literal_margin_increment":INCREMENT,
        "added_physical_lower_scalar":lower,"added_hat_volume_Lagrangian_density":density,
        "added_hat_volume_Hamiltonian_density":-density,
        "added_rationalized_potential":deltaV,
        "clock_fixed_phase_lapse_Hessian_change":-8*INCREMENT/c["h"]**2,
        "separately_named_light_action_not_an_edit_of_any_frozen_parent":True,
        "no_old_quantum_profile_state_or_UV_matching_transferred":True}


@cache
def coefficients():
    c=dict(old.coefficients())
    c["V0"]=sp.factor(c["V0"]+action()["added_rationalized_potential"])
    c["C"]=sp.factor(-sp.diff(c["V0"],N)+c["Ln"]*c["V0"]
                    +sp.diff(c["I0"],u)-c["Lt"]*c["I0"])
    return c


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
    Ndot=-forcing/pivot
    zdot=Z0+Kz*Ndot
    Bbar=-(2*A*z+B)/(3*M)
    H=-N*z/2
    G0=Bbar.diff(u)+Bbar.diff(z)*Z0-(Lt+H)*Bbar-m0*Bbar*Bbar-N*w2
    GN=Bbar.diff(N)+Bbar.diff(z)*Kz-Ln*Bbar
    speed=-4*h*M*M*(G0*pivot-GN*forcing)/(D*pivot*pivot)
    return {**data,"constraint_matter_square":w2,"trace_lapse_chain_coefficient":Kz,
        "trace_flow_at_fixed_lapse":Z0,"on_constraint_rescaled_lapse_Hessian":pivot,
        "constraint_time_forcing":forcing,"lapse_flow":Ndot,"rescaled_trace_flow":zdot,
        "rescaled_principal_B":Bbar,"hat_Hubble":H,
        "conformal_log_flow":Lt+Ln*Ndot,
        "gradient_fixed_lapse_part":G0,"gradient_lapse_velocity_part":GN,
        "clock_physical_speed_squared":speed,"clock_physical_subluminal_margin":1-speed}


@cache
def checks():
    a,c,d=action(),coefficients(),system()
    X=a["physical_clock_X"]
    lower=a["added_physical_lower_scalar"]
    density=a["added_hat_volume_Lagrangian_density"]
    before=old.coefficients()
    out={"exact_new_total_margin_is_old_plus_increment":NEW_MARGIN-OLD_MARGIN-INCREMENT,
        "physical_X_and_older_negative_x_conventions_agree":sp.factor(lower-
            original_margin.deformation()["physical_lower_scalar"].subs(
                {original_margin.x:-X,original_margin.epsilon:INCREMENT},simultaneous=True)),
        "added_scalar_value_vanishes_on_entire_clock_trajectory":lower.subs(X,1),
        "added_scalar_clock_first_derivative_vanishes":sp.diff(lower,X).subs(X,1),
        "added_scalar_explicit_time_first_derivative_vanishes":sp.diff(lower,u).subs(X,1),
        "added_scalar_clock_first_derivative_total_time_derivative_vanishes":
            sp.diff(sp.diff(lower,X).subs(X,1),u),
        "literal_density_value_vanishes_on_clock":density.subs(N,1),
        "literal_density_lapse_first_vanishes_on_clock":sp.diff(density,N).subs(N,1),
        "literal_Hamiltonian_clock_lapse_Hessian_change_retained":
            sp.factor(-sp.diff(density,N,2).subs(N,1)-a["clock_fixed_phase_lapse_Hessian_change"]),
        "conformal_rationalization_retains_all_four_powers_of_e":sp.factor(
            a["added_rationalized_potential"]-INCREMENT*(N*N-1)**2/(N*c["D"]*c["h"])),
        "literal_new_constraint_constant_derived_before_matter_elimination":sp.factor(
            c["C"]-before["C"]+sp.diff(a["added_rationalized_potential"],N)
            -c["Ln"]*a["added_rationalized_potential"]),
        "new_lapse_flow_exactly_preserves_full_constraint":(
            d["constraint_time_forcing"]+d["on_constraint_rescaled_lapse_Hessian"]*d["lapse_flow"]).as_expr(),
        "new_trace_flow_keeps_the_actual_lapse_velocity":(
            d["rescaled_trace_flow"]-d["trace_flow_at_fixed_lapse"]
            -d["trace_lapse_chain_coefficient"]*d["lapse_flow"]).as_expr()}
    for name in ("e","Ln","Lt","a0","m0","d0","I0","b0","A","M","B","h","D"):
        out["unchanged_literal_derivative_coefficient_"+name]=sp.factor(c[name]-before[name])
    return {name:sp.factor(value) for name,value in out.items()}
