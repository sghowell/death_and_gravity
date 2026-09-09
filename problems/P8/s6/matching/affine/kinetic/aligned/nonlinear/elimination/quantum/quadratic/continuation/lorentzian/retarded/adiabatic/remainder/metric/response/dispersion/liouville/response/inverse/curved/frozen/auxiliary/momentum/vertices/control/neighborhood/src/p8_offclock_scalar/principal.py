"""Exact central-slice gamma phase blocks and their high-frequency orders."""
from functools import cache

import sympy as sp

from . import datum, spatial

q=sp.Symbol("positive_hat_squared_momentum",positive=True)


def clean(value):
    return spatial.root_relation(sp.factor(value.subs(spatial.k,sp.sqrt(q))))


def infinity(value):
    numerator,denominator=sp.fraction(sp.cancel(value))
    a,b=sp.Poly(numerator,q),sp.Poly(denominator,q)
    if a.is_zero or a.degree()<b.degree():
        return sp.Integer(0)
    if a.degree()>b.degree():
        raise ValueError("Unbounded high-frequency limit")
    return spatial.root_relation(a.LC()/b.LC())


@cache
def blocks():
    old=spatial.phase_matrix()
    T=sp.zeros(6)
    T[0,3]=1/(2*q)
    T[1,0]=-2*q
    T[2,1]=T[3,4]=T[4,2]=T[5,5]=1
    transformed=(T.T*old*T).applyfunc(clean)
    A=sp.ImmutableMatrix(transformed[3:,3:])
    B=sp.ImmutableMatrix(transformed[3:,:3])
    C=sp.ImmutableMatrix(transformed[:3,:3])
    return {"old_phase_matrix":old,"canonical_gamma_map":T,
            "momentum_Hessian":A,"momentum_coordinate_block":B,"coordinate_Hessian":C,
            "momentum_Hessian_limit":sp.ImmutableMatrix(A.applyfunc(infinity)),
            "mixing_over_q_limit":sp.ImmutableMatrix((B/q).applyfunc(infinity))}


@cache
def checks():
    d,b=datum.data(),blocks()
    n=spatial.t**2
    h=spatial.at_lapse(d["fixed_phase_lapse_Hessian"])
    r=-1/(4*n**3*h)
    expected=sp.Matrix([[r,-r*spatial.P,0],[-r*spatial.P,1/spatial.t+r*spatial.P**2,0],
                        [0,0,spatial.at_lapse(d["longitudinal_kinetic_times_sqrt_N"])/spatial.t]])
    expectedB=sp.zeros(3)
    expectedB[2,0]=spatial.at_lapse(d["canonical_gamma_vector_mixing"])
    oldJ=sp.diag(*([sp.Matrix([[0,1],[-1,0]])]*3))
    newJ=sp.zeros(6)
    newJ[:3,3:]=sp.eye(3)
    newJ[3:,:3]=-sp.eye(3)
    A,B,C=b["momentum_Hessian"],b["momentum_coordinate_block"],b["coordinate_Hessian"]
    out={"full_spatial_gamma_canonical_map":(b["canonical_gamma_map"].T*oldJ*b["canonical_gamma_map"]-newJ).applyfunc(sp.factor),
         "actual_three_scalar_high_q_kinetic_Hessian":(b["momentum_Hessian_limit"]-expected).applyfunc(clean),
         "actual_scalar_longitudinal_q_mixing":(B-q*expectedB).applyfunc(clean),
         "spatial_constraint_cancels_bare_scalar_momentum_square":spatial.pair("scalar_p","scalar_p")["kernel"],
         "no_linear_lapse_vector_force_at_actual_datum":spatial.pair("longitudinal_p","longitudinal_p")["lapse_force_left"],
         "no_bare_gamma_b_squared_term":C[0,0],
         "no_bare_gamma_b_longitudinal_coordinate_term":C[0,2],
         "no_bare_matter_longitudinal_coordinate_term":C[1,2],
         "actual_matter_gradient_leading_coefficient":infinity(C[1,1]/q)-spatial.t**3,
         "actual_longitudinal_gradient_coefficient":clean(C[2,2]/q-spatial.t**3/spatial.at_lapse(d["gamma_s"])),
         "actual_longitudinal_finite_q_kinetic_correction":clean(
             A[2,2]-expected[2,2]-10**6*spatial.t/q),
         "leading_light_kinetic_determinant":clean(expected[:2,:2].det()-r/spatial.t)}
    return out
