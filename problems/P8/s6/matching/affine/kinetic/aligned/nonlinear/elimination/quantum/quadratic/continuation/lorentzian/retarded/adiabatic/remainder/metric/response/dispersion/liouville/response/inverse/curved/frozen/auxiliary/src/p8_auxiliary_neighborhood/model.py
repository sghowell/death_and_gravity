"""Full post-temporal Hamiltonian in nine geometric canonical-jet invariants."""
from functools import cache

import sympy as sp
from p8_affine import lower
from p8_affine_nonlinear import adm, trace
from p8_aligned_margin import action as margin

u,N,s=adm.u,adm.N,adm.s
MARGIN=sp.Rational(1,10**6)
I,Iphi,Q=sp.symbols("boundary_primitive boundary_primitive_phi original_Q_lower",real=True)
dp,dc,j,shear,electric,magnetic,vector,chi_gradient,curvature=sp.symbols(
    "trace_momentum_deviation matter_momentum_deviation vector_momentum_divergence "
    "squared_shear_momentum normalized_squared_electric_momentum "
    "normalized_squared_magnetic_field squared_vector_potential "
    "squared_matter_gradient spatial_curvature",real=True)
VARIABLES=(dp,dc,j,shear,electric,magnetic,vector,chi_gradient,curvature)


@cache
def gamma():
    p=sp.Symbol("affine_p",positive=True)
    return {"p":p,"temporal":9*(2*p**3-1)/(22*p**3+3*p-11),
            "spatial":9*(8*p+5)/(-72*p**2+88*p+55)}


@cache
def coefficients():
    old=adm.transform()
    bg=trace.alignment.parent.old.background()
    h=bg["h"]
    ratio=(h-1+N**-2)/h
    eo=ratio**(-sp.Rational(1,4))
    U=ratio**(-sp.Rational(3,4))
    B=old["B"].subs(s,1/N)
    linear=U*old["Blinear"].subs(s,1/N)-I
    margin_scalar=margin.deformation()["physical_lower_scalar"].subs(
        {margin.x:-N**-2,margin.epsilon:MARGIN},simultaneous=True)
    potential=U*(old["Fnew"].subs(s,1/N)+margin_scalar)-Iphi/N
    source_delta=(N**-2-1)/h
    p=sp.sqrt(ratio)/2
    mass=gamma()
    return {"background":bg,"ratio":ratio,"U":U,"B":B,"a":sp.Rational(2,3)*U*B,
            "b":linear,"f0":potential,"delta":source_delta,
            "c":-3*bg["H"]*source_delta/N+sp.Rational(3,2)*Q/N,
            "gamma_t":mass["temporal"].subs(mass["p"],p),
            "gamma_s":mass["spatial"].subs(mass["p"],p),
            "e_omega":eo,"B4":-B,"boundary_primitive_s":old["primitive_s"],
            "Q_ODE_coefficient":lower.ode()["coefficient"],
            "Q_ODE_forcing":lower.ode()["forcing"]}


@cache
def generic():
    a,U,eo,gt,gs,B4=sp.symbols("trace_a hat_volume e_omega gamma_t gamma_s B4",nonzero=True)
    b,f,d,c,p0,l0=sp.symbols("trace_b scalar_f source_delta source_c background_trace_momentum background_matter_momentum",real=True)
    UB=sp.Rational(3,2)*a
    p,l=p0+dp,l0+dc
    trace_part=(p-b-d*j)**2/(4*a)-f+gt*j*j/(2*U)-c*j
    pieces={"trace_potential_and_temporal_Gauss":trace_part,
            "five_shear":-shear/UB,
            "electric":electric/(2*eo),
            "magnetic":magnetic/(4*eo),
            "spatial_vector_mass":eo*vector/(2*gs),
            "matter_momentum":l*l/(2*U),
            "matter_gradient":eo*chi_gradient/2,
            "spatial_curvature":-eo*B4*curvature}
    density=N*sum(pieces.values())
    K=(p-b-d*j)/(2*a)
    T=d*K+c-gt*j/U
    return {"a":a,"U":U,"eo":eo,"gt":gt,"gs":gs,"B4":B4,"b":b,"f":f,
            "d":d,"c":c,"p0":p0,"l0":l0,"pieces":pieces,
            "Hamiltonian_over_hat_volume":density,
            "trace_reconstruction":K,"temporal_reconstruction":T}


@cache
def specialization():
    g,c=generic(),coefficients()
    return {g["a"]:c["a"],g["U"]:c["U"],g["eo"]:c["e_omega"],
            g["gt"]:c["gamma_t"],g["gs"]:c["gamma_s"],g["B4"]:c["B4"],
            g["b"]:c["b"],g["f"]:c["f0"],g["d"]:c["delta"],g["c"]:c["c"],
            g["p0"]:-2*c["background"]["H"],g["l0"]:c["background"]["ell"]}


@cache
def checks():
    g=generic()
    original=trace.legendre()
    mapping={original["a"]:g["a"],original["volume"]:g["U"],original["gamma"]:g["gt"],
             original["b"]:g["b"],original["f"]:g["f"],original["delta"]:g["d"],
             original["c"]:g["c"],original["p"]:g["p0"]+dp,original["j"]:j}
    out={"full_joint_trace_temporal_Hamiltonian_bridge":sp.factor(
        original["H_after_temporal"].subs(mapping,simultaneous=True)-g["pieces"]["trace_potential_and_temporal_Gauss"]),
         "full_joint_temporal_reconstruction_bridge":sp.factor(
        original["T_joint"].subs(mapping,simultaneous=True)-g["temporal_reconstruction"]),
         "full_joint_trace_reconstruction_bridge":sp.factor(
        original["K_joint"].subs(mapping,simultaneous=True)-g["trace_reconstruction"])}
    zeta=sp.Symbol("positive_Maxwell_coefficient",positive=True)
    E,Bmag=sp.symbols("unscaled_squared_electric_momentum unscaled_squared_magnetic_field",real=True)
    out["electric_invariant_normalization_keeps_inverse_zeta"]=sp.factor(
        g["pieces"]["electric"].subs(electric,E/zeta)-E/(2*zeta*g["eo"]))
    out["magnetic_invariant_normalization_keeps_zeta"]=sp.factor(
        g["pieces"]["magnetic"].subs(magnetic,zeta*Bmag)-zeta*Bmag/(4*g["eo"]))
    p=gamma()["p"]
    out["literal_retuned_temporal_mass_rational_function"]=sp.factor(
        gamma()["temporal"]-1/(p/(3*(2*p**3-1))+sp.Rational(11,9)))
    out["literal_retuned_spatial_mass_rational_function"]=sp.factor(
        gamma()["spatial"]-1/(sp.Rational(11,9)-8*p*p/(8*p+5)))
    c=coefficients()
    bg=c["background"]
    old=adm.clock_lapse()
    out["old_scalar_background_equation_replayed"]=old["actual_background_lapse_equation"]
    out["old_scalar_lapse_Hessian_replayed"]=old["actual_secondary_lapse_Jacobian"]
    deform=margin.deformation()
    out["margin_lapse_Hessian_shift"]=sp.factor(
        -sp.diff(deform["hat_density"],margin.N,2).subs({margin.N:1,margin.epsilon:MARGIN})
        +8*MARGIN/bg["h"]**2)
    zero={N:1,I:0,Iphi:0,Q:0}
    out["clock_linear_trace_coefficient_zero"]=sp.factor(c["b"].subs(zero))
    out["clock_temporal_source_zero"]=sp.factor(c["c"].subs(zero))
    out["clock_trace_reconstruction"]=sp.factor(
        g["trace_reconstruction"].subs(specialization(),simultaneous=True).subs(
            dict.fromkeys(VARIABLES,0)).subs(zero)-3*bg["H"])
    out["clock_temporal_reconstruction"]=sp.factor(
        g["temporal_reconstruction"].subs(specialization(),simultaneous=True).subs(
            dict.fromkeys(VARIABLES,0)).subs(zero))
    return out
