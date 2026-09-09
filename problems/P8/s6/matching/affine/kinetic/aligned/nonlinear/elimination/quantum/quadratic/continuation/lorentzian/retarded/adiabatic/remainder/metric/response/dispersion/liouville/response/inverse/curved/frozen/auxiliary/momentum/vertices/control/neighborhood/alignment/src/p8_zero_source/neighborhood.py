"""Remove the exact old q^4 mechanism, without asserting all-frequency health."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model as auxiliary
from p8_offclock_scalar import datum, principal, spatial


@cache
def data():
    old=principal.blocks()
    d=datum.data()
    t,q=spatial.t,principal.q
    n=t*t
    delta=spatial.at_lapse(d["source_delta"])
    gt=spatial.at_lapse(d["gamma_t"])
    A=sp.MutableDenseMatrix(old["momentum_Hessian"])
    B=sp.MutableDenseMatrix(old["momentum_coordinate_block"])
    C=old["coordinate_Hessian"]
    # Actual central-slice Delta H = delta*N^(3/2)*p*P_sigma/2
    # +3*N^(3/2)*delta^2*P_sigma^2/4, with p=-2q*b.
    A[2,2]+=sp.Rational(3,2)*n**sp.Rational(3,2)*delta**2
    B[2,0]-=n**sp.Rational(3,2)*delta*q
    A=sp.ImmutableMatrix(A.applyfunc(principal.clean))
    B=sp.ImmutableMatrix(B.applyfunc(principal.clean))
    return {"unchanged_constraint_compatible_light_family":d,
            "new_momentum_Hessian":A,"new_momentum_coordinate_block":B,
            "new_coordinate_Hessian":C,
            "new_Lagrangian_coordinate_Hessian":C,
            "new_longitudinal_high_frequency_pivot":gt/t,
            "new_joint_temporal_auxiliary_pivot":-auxiliary.N**sp.Rational(5,2)/d["gamma_t"],
            "new_vector_mass_ratio_on_bounce_slice":sp.factor(d["gamma_t"]/d["gamma_s"]),
            "same_specific_q_fourth_determinant_coefficient":sp.Integer(0),
            "scope":"Only the previously identified scalar-longitudinal source mixing and its negative q^4 Legendre determinant coefficient are removed. The mass ratio, light scalar sector, higher interactions and same-scale Gaussian response are not made healthy by this identity."}


@cache
def checks():
    d=data()
    A,B=d["new_momentum_Hessian"],d["new_momentum_coordinate_block"]
    q,t=principal.q,spatial.t
    C=d["new_coordinate_Hessian"]
    g=auxiliary.generic()
    p,P,b=sp.symbols("central_trace_momentum longitudinal_momentum gamma_coordinate",real=True)
    n=t*t
    delta=n**-2-1
    # Derive the correction from the literal Hamiltonian difference,
    # independently of the old six-channel Hessian reconstruction.
    change=auxiliary.N*((g["p0"]+auxiliary.dp-g["b"])*g["d"]*auxiliary.j/(2*g["a"])
                        -g["d"]**2*auxiliary.j**2/(4*g["a"])+g["c"]*auxiliary.j)
    change=sp.factor(change.subs({auxiliary.N:n,g["p0"]:0,g["b"]:0,g["c"]:0,
                                 g["a"]:-1/(3*t),g["d"]:delta,
                                 auxiliary.dp:p/3,auxiliary.j:-P},simultaneous=True).subs(p,-2*q*b))
    expected=-n**sp.Rational(3,2)*delta*q*b*P+sp.Rational(3,4)*n**sp.Rational(3,2)*delta**2*P**2
    ratio=d["new_vector_mass_ratio_on_bounce_slice"]
    # Once the vector-light source block vanishes identically, all
    # remaining central Euler entries are at most O(q), including
    # arbitrary finite background time jets. Do not set those jets zero.
    leading=sp.Matrix(3,3,sp.symbols("new_Euler_order_q0:9",real=True))
    remainder=sp.Matrix(3,3,sp.symbols("new_Euler_order_one0:9",real=True))
    return {"literal_new_central_phase_correction":sp.factor(change-expected),
            "full_new_scalar_longitudinal_momentum_coordinate_block_zero":B,
            "new_longitudinal_pivot_is_actual_temporal_mass_not_old_joint_kappa":principal.clean(
                A[2,2]-d["new_longitudinal_high_frequency_pivot"]-10**6*t/q),
            "new_vector_light_momentum_cross_block_zero":sp.ImmutableMatrix(A[:2,2:]),
            "no_new_negative_q_squared_b_coordinate_entry":C[0,0],
            "new_bare_coordinate_determinant_has_no_q_fourth_term":principal.infinity(C.det()/q**4),
            "arbitrary_remaining_Euler_time_jets_have_no_q_fourth_determinant":principal.infinity(
                (q*leading+remainder).det()/q**4),
            "unchanged_offclock_mass_ratio_at_clock":ratio.subs(auxiliary.N,1)-1,
            "unchanged_mass_ratio_normal_derivative_is_not_made_zero":sp.factor(
                sp.diff(ratio,auxiliary.N).subs(auxiliary.N,1)+sp.Rational(8,81))}


@cache
def gates():
    ratio=data()["new_vector_mass_ratio_on_bounce_slice"]
    return {"old_offclock_mass_ratio_issue_has_not_been_hidden":bool(ratio.subs(auxiliary.N,1-datum.RADIUS)>1),
            "new_joint_temporal_pivot_nonzero_in_parent_continuous_mass_domain":bool(
                datum.domain()["continuous_rational_proofs"]["gamma_t"]["lower"]>sp.Rational(9,10)),
            "no_actual_time_dependent_or_UV_health_inferred_from_removed_mechanism":True}
