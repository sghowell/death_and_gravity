"""Exact same-action homogeneous data at the central clock slice."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import center, model
from p8_coupled_energy.bounds import exact_expression

N=model.N
P=sp.Symbol("positive_homogeneous_matter_momentum",positive=True)
RADIUS=sp.Rational(1,10**28)


def enclosure(value):
    value=exact_expression(value,(N,))
    x=sp.Symbol("lapse_displacement",real=True)
    numerator,denominator=sp.fraction(sp.cancel(value))
    num=sp.Poly(numerator.subs(N,1+x),x,domain=sp.QQ)
    den=sp.Poly(denominator.subs(N,1+x),x,domain=sp.QQ)
    if den.nth(0)<0:
        num,den=-num,-den
    tail=lambda poly:sum(abs(c)*RADIUS**degree[0] for degree,c in poly.terms() if degree[0])
    nt,dt=tail(num),tail(den)
    lower,upper=den.nth(0)-dt,den.nth(0)+dt
    if lower<=0:
        raise ValueError("No nonzero continuous lapse-box denominator")
    endpoints=((num.nth(0)-nt)/lower,(num.nth(0)-nt)/upper,
               (num.nth(0)+nt)/lower,(num.nth(0)+nt)/upper)
    return {"lower":min(endpoints),"upper":max(endpoints),
            "numerator_center":num.nth(0),"numerator_tail_upper":nt,
            "denominator_center":den.nth(0),"denominator_tail_upper":dt,
            "radius":RADIUS}


@cache
def data():
    old=center.data()
    A=old["background"]
    H=A+(P**2-sp.Rational(1,100))/(2*sp.sqrt(N))
    P2=sp.factor(sp.Rational(1,100)+4*N**sp.Rational(3,2)*sp.diff(A,N))
    hessian=sp.factor(sp.diff(H,N,2).subs(P**2,P2))
    gt,gs=old["actual_gamma_t_at_bounce"],old["actual_gamma_s_at_bounce"]
    delta=N**-2-1
    kappa=sp.factor(gt-sp.Rational(3,2)*N**2*delta**2)
    return {"homogeneous_Hamiltonian":H,"matter_momentum_squared":P2,
            "fixed_phase_lapse_Hessian":hessian,
            "lapse_Hessian_times_N_seven_halves":sp.factor(N**sp.Rational(7,2)*hessian),
            "gamma_t":gt,"gamma_s":gs,"source_delta":delta,
            "longitudinal_kinetic_times_sqrt_N":kappa,
            "joint_temporal_auxiliary_pivot":-N**sp.Rational(5,2)/kappa,
            "canonical_gamma_vector_mixing":sp.factor(N**sp.Rational(3,2)*delta),
            "physical_spatial_metric_over_hat_metric":N,
            "physical_lapse":N,
            "datum_description":"phi=0, flat hat spatial metric I, zero homogeneous gravitational momentum and spatial vector data, positive matter momentum P(N). The lapse constraint is imposed at fixed canonical phase data, not differentiated along the constraint family."}


@cache
def checks():
    d=data()
    H,P2,h=d["homogeneous_Hamiltonian"],d["matter_momentum_squared"],d["fixed_phase_lapse_Hessian"]
    eps=model.MARGIN
    expectedP=(sp.Rational(69,4)+6*eps)/N-(sp.Rational(2849,100)-4*eps)*N+(sp.Rational(45,4)-10*eps)*N**3
    a1,a2,a3=sp.Rational(23,8)+eps,sp.Rational(2849,200)-2*eps,sp.Rational(9,8)-eps
    expectedh=(-3*a1-a2*N**2+15*a3*N**4)/(2*N**sp.Rational(7,2))
    T,trace_p=sp.symbols("temporal_vector trace_momentum",real=True)
    a=-1/(3*sp.sqrt(N))
    U=N**sp.Rational(3,2)
    delta,gt=d["source_delta"],d["gamma_t"]
    K=(trace_p+U*delta*T/gt)/(2*a+U*delta**2/gt)
    before=N*(trace_p*K-a*K*K-U*(T-delta*K)**2/(2*gt))
    pivot=sp.factor(sp.diff(before,T,2))
    return {"actual_homogeneous_lapse_constraint_family":sp.factor(sp.diff(H,N).subs(P**2,P2)),
            "actual_positive_matter_momentum_squared_formula":sp.factor(P2-expectedP),
            "fixed_phase_not_constraint_tangent_lapse_Hessian":sp.factor(h-expectedh),
            "original_clock_matter_datum":P2.subs(N,1)-sp.Rational(1,100),
            "original_clock_lapse_Hessian":h.subs(N,1)+sp.Rational(1199,400)+8*eps,
            "original_clock_vector_kinetic_pivot":d["longitudinal_kinetic_times_sqrt_N"].subs(N,1)-1,
            "full_joint_temporal_pivot_not_isolated_mass":sp.factor(pivot-d["joint_temporal_auxiliary_pivot"]),
            "joint_mixed_lapse_temporal_pivot_zero_at_actual_datum":sp.factor(
                sp.diff(before,T,N).subs({T:0,trace_p:0})),
            "offclock_gamma_vector_coupling":sp.factor(d["canonical_gamma_vector_mixing"]-(1-N*N)/sp.sqrt(N)),
            "nonzero_first_normal_derivative_of_coupling":sp.diff(d["canonical_gamma_vector_mixing"],N).subs(N,1)+2}


@cache
def domain():
    d=data()
    quantities={"matter_squared":d["matter_momentum_squared"],
                "matter_squared_secant":sp.cancel((d["matter_momentum_squared"]-sp.Rational(1,100))/(N-1)),
                "normalized_lapse_Hessian":d["lapse_Hessian_times_N_seven_halves"],
                "gamma_t":d["gamma_t"],"gamma_s":d["gamma_s"],
                "normalized_longitudinal_kinetic":d["longitudinal_kinetic_times_sqrt_N"]}
    proofs={name:enclosure(value) for name,value in quantities.items()}
    secant=max(abs(proofs["matter_squared_secant"][name]) for name in ("lower","upper"))
    gates={"positive_matter_root":proofs["matter_squared"]["lower"]>sp.Rational(1,200),
           "fixed_phase_lapse_Hessian_negative":proofs["normalized_lapse_Hessian"]["upper"]<-2,
           "finite_temporal_mass_pivot":proofs["gamma_t"]["lower"]>sp.Rational(9,10),
           "finite_spatial_mass_pivot":proofs["gamma_s"]["lower"]>sp.Rational(9,10),
           "positive_longitudinal_kinetic_limit":proofs["normalized_longitudinal_kinetic"]["lower"]>sp.Rational(9,10),
           "inside_original_auxiliary_lapse_disc":RADIUS<sp.Rational(1,10**14),
           "inside_original_nine_invariant_polydisc":10*secant*RADIUS<sp.Rational(1,10**24)}
    if not all(bool(value) for value in gates.values()):
        raise ValueError("The same-action off-clock datum left its certified auxiliary domain")
    return {"lapse_half_width":RADIUS,"positive_matter_momentum_deviation_upper":10*secant*RADIUS,
            "continuous_rational_proofs":proofs,"gates":{name:bool(value) for name,value in gates.items()}}
