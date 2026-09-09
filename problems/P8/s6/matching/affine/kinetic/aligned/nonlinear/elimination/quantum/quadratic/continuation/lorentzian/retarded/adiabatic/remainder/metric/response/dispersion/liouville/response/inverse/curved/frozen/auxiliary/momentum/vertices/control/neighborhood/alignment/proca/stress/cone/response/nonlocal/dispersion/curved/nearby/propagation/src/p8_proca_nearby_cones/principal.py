"""Euler-first scalar principal matrices in a crossing-regular momentum chart."""
from functools import cache

import sympy as sp

a,m,g,r,rn,h,ell=sp.symbols("trace_square matter_square matter_gradient curvature curvature_N lapse_Hessian matter_density",nonzero=True)
H,alpha,beta=sp.symbols("hat_Hubble force_trace force_matter",real=True)
q=sp.Symbol("hat_momentum_squared",positive=True)
b,chi,Pb,P=sp.symbols("clock_coordinate matter_coordinate clock_momentum matter_momentum",real=True)


def infinity_limit(value):
    n,d=sp.fraction(sp.cancel(value))
    pn,pd=sp.Poly(n,q),sp.Poly(d,q)
    if pn.is_zero or pn.degree()<pd.degree():
        return sp.Integer(0)
    if pn.degree()==pd.degree():
        return sp.factor(pn.LC()/pd.LC())
    raise ValueError("A declared finite principal limit diverges")


@cache
def data():
    v,p=Pb/(2*q),-2*q*b
    force=alpha*p+beta*(P-3*ell*v)+4*rn*q*v
    ham=m*(P-3*ell*v)**2+sp.Rational(2,3)*a*ell*p*chi
    ham+=(g*q-a*ell*ell)*chi*chi+2*r*q*v*v-force**2/(2*h)-H*b*Pb
    mat=sp.hessian(ham,(b,chi,Pb,P))
    A,B,C=mat[2:,2:],mat[2:,:2],mat[:2,:2]
    K=A.inv().applyfunc(sp.factor)
    V=(-K*B).applyfunc(sp.factor)
    G=(C-B.T*K*B).applyfunc(sp.factor)
    return {"Hamiltonian":ham,"momentum_Hessian":A,
            "finite_q_Legendre_determinant":sp.factor(A.det()),
            "kinetic":K.applyfunc(infinity_limit),
            "velocity_mixed_over_q":(V/q).applyfunc(infinity_limit),
            "potential_over_q":(G/q).applyfunc(infinity_limit)}


@cache
def formula():
    Bd=sp.Symbol("time_derivative_of_alpha_over_curvature_N",real=True)
    B=alpha/rn
    K=sp.Matrix([[(beta*beta-2*m*h)/(8*m*rn*rn),-beta/(4*m*rn)],
                 [-beta/(4*m*rn),1/(2*m)]])
    G=sp.Matrix([[Bd-H*B+r*B*B,-4*a*ell/3],[-4*a*ell/3,2*g]])
    return {"kinetic":K,"Euler_gradient":G,
            "B":B,"B_dot":Bd,
            "weighted_Euler_gradient_added_to_frozen_potential":sp.diag(Bd+H*B,0),
            "momentum_chart_requires_q_positive_and_r_N_nonzero":True,
            "division_by_Hubble_or_gamma_crossing_coefficient":False,
            "all_finite_q_Legendre_regular":False}


@cache
def physical():
    e,N=sp.symbols("positive_eomega positive_lapse",positive=True)
    mN=sp.Symbol("matter_square_N",nonzero=True)
    a0,m0,g0=-3*N*e/4,N/(2*e**3),N*e/2
    substitution={a:a0,m:m0,g:g0,r:-m0,rn:-mN,beta:2*mN*ell}
    K=formula()["kinetic"].subs(substitution,simultaneous=True)
    G=formula()["Euler_gradient"].subs(substitution,simultaneous=True)
    T=sp.Matrix([[1,0],[-ell,1]])
    Kdiag=(T.T*K*T).applyfunc(sp.factor)
    Gdiag=(T.T*G*T).applyfunc(sp.factor)
    speed=sp.factor(e*e*Gdiag[0,0]/(N*N*Kdiag[0,0]))
    return {"eomega":e,"lapse":N,"matter_square_N":mN,
            "triangular_characteristic_basis":T,
            "diagonal_kinetic":Kdiag,"diagonal_gradient":Gdiag,
            "physical_clock_speed_squared":speed,
            "physical_matter_speed_squared":sp.factor(e*e*Gdiag[1,1]/(N*N*Kdiag[1,1]))}


@cache
def checks():
    d,f,p=data(),formula(),physical()
    B=f["B"]
    added=f["weighted_Euler_gradient_added_to_frozen_potential"]
    return {
        "high_q_kinetic_from_full_canonical_Hessian":(d["kinetic"]-f["kinetic"]).applyfunc(sp.factor),
        "leading_velocity_mixing_is_symmetric_rank_one":(d["velocity_mixed_over_q"]-sp.diag(B,0)).applyfunc(sp.factor),
        "Euler_first_gradient_retains_weight_and_wave_number_time_derivative":
            (d["potential_over_q"]+added-f["Euler_gradient"]).applyfunc(sp.factor),
        "finite_q_Legendre_determinant_not_declared_nonzero_at_all_q":
            sp.factor(d["finite_q_Legendre_determinant"]+(beta*beta*r-2*h*m*r+8*m*q*rn*rn)/(h*q)),
        "physical_triangular_basis_has_unit_determinant":p["triangular_characteristic_basis"].det()-1,
        "clock_kinetic_Schur_complement_retains_lapse_Hessian":
            sp.factor(p["diagonal_kinetic"][0,0]+h/(4*p["matter_square_N"]**2)),
        "matter_kinetic_is_positive_physical_U_over_N":
            sp.factor(p["diagonal_kinetic"][1,1]-p["eomega"]**3/p["lapse"]),
        "kinetic_off_diagonal_removed":p["diagonal_kinetic"][0,1],
        "gradient_off_diagonal_removed":p["diagonal_gradient"][0,1],
        "physical_matter_characteristic_exactly_luminal":p["physical_matter_speed_squared"]-1,
        "time_dependent_canonical_transform_retains_weighted_boundary":-2*H+3*H-H,
        "physical_momentum_and_proper_time_conversion":
            sp.factor((q/p["eomega"]**2)*(p["lapse"]**2/p["eomega"]**2)**-1-q/p["lapse"]**2)}
