"""Euler-first high-frequency determinant, with arbitrary finite time jets."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model
from p8_coupled_energy import scalars

from . import datum, principal, spatial


@cache
def actual_orders():
    d,b=datum.data(),principal.blocks()
    q=principal.q
    A,B,C=b["momentum_Hessian"],b["momentum_coordinate_block"],b["coordinate_Hessian"]
    coupling=spatial.at_lapse(d["canonical_gamma_vector_mixing"])
    gamma=sp.MutableDenseMatrix(C)
    gamma[0,0]-=coupling**2*q*q/A[2,2]
    return {"Lagrangian_coordinate_matrix":sp.ImmutableMatrix(gamma),
            "minus_q_squared_coefficient":spatial.at_lapse(
                d["canonical_gamma_vector_mixing"]**2*sp.sqrt(model.N)/d["longitudinal_kinetic_times_sqrt_N"]),
            "matter_gradient_coefficient":spatial.t**3,
            "longitudinal_gradient_coefficient":spatial.t**3/spatial.at_lapse(d["gamma_s"]),
            "Hamiltonian_mixing_block":B,
            "time_jet_order_statement":"A and A_dot are O(1), with positive invertible A_infinity. B and B_dot have first column O(q) and the other columns O(1). Therefore beta_dot and 3*H_hat*beta have only O(q) in their first column and O(1) elsewhere. All actual homogeneous time jets are finite; they are not set to zero."}


@cache
def generic_determinant():
    q=principal.q
    d,av,ev,cchi,cv=sp.symbols("nonzero_mixing positive_vector_pivot electric_remainder positive_matter_gradient positive_vector_gradient",real=True)
    a00,a01,a10,a20=sp.symbols("first_order00 first_order01 first_order10 first_order20",real=True)
    rest=sp.Matrix(3,3,sp.symbols("bounded_remainder0:9",real=True))
    C=rest+sp.Matrix([[-d*d*q*q/(av+ev/q)+a00*q,a01*q,0],
                      [a10*q,cchi*q,0],[a20*q,0,cv*q]])
    determinant=sp.factor(C.det())
    leading=principal.infinity(determinant/q**4)
    return {"matrix_with_all_allowed_time_jet_orders":C,
            "leading_determinant_over_q_fourth":leading,
            "expected_leading_coefficient":-d*d*cchi*cv/av,
            "residual":sp.factor(leading+d*d*cchi*cv/av)}


@cache
def checks():
    d,b=datum.data(),principal.blocks()
    q=principal.q
    gamma=actual_orders()["Lagrangian_coordinate_matrix"]
    A,C=b["momentum_Hessian"],b["coordinate_Hessian"]
    out={"arbitrary_finite_Euler_time_jets_do_not_remove_q_fourth_determinant":generic_determinant()["residual"],
         "actual_Lagrangian_q_squared_negative_entry":principal.clean(
             principal.infinity(gamma[0,0]/q**2)+actual_orders()["minus_q_squared_coefficient"]),
         "actual_matter_gradient_has_no_higher_order_remainder":principal.infinity(
             (C[1,1]-spatial.t**3*q)/q),
         "actual_b_matter_coordinate_coupling_is_at_most_q":principal.infinity(C[0,1]/q**2),
         "actual_vector_kinetic_block_decouples_at_slice":sp.ImmutableMatrix(A[:2,2:]).applyfunc(principal.clean)}
    coefficient=sp.factor(-actual_orders()["minus_q_squared_coefficient"]*
                          actual_orders()["matter_gradient_coefficient"]*
                          actual_orders()["longitudinal_gradient_coefficient"])
    target=spatial.at_lapse(-(1-model.N**2)**2*model.N**sp.Rational(5,2)/(
        d["gamma_s"]*d["longitudinal_kinetic_times_sqrt_N"]))
    out["actual_nonzero_offclock_negative_determinant_coefficient"]=principal.clean(coefficient-target)
    clock=gamma.subs({spatial.t:1,spatial.P:sp.Rational(1,10)})
    out["freezing_Hamiltonian_before_Euler_would_give_false_clock_sign"]=sp.factor(
        principal.infinity(clock.det()/q**3)+sp.Rational(1,100))
    original=scalars.data("gamma")
    actual=original["potential"]+original["antisymmetric_mixing"].applyfunc(scalars.derivative)+3*scalars.H*original["antisymmetric_mixing"]
    actual=actual.subs(scalars.actual_background_jets(),simultaneous=True).subs(
        {model.u:0,scalars.z:1/q}).applyfunc(sp.factor)
    # The original selected clock scalar operator, varied before freezing,
    # has positive q^2 determinant; the longitudinal block adds q.
    out["actual_Euler_first_clock_sign_is_positive_not_false_frozen_H_sign"]=sp.factor(
        principal.infinity(actual.det()/q**2)-sp.Rational(1199,100))
    return out


@cache
def result():
    d=datum.data()
    return {"same_action_lapse_domain":datum.domain(),
            "Euler_constant_determinant_leading_coefficient":sp.factor(
                -(1-model.N**2)**2*model.N**sp.Rational(5,2)/(
                    d["gamma_s"]*d["longitudinal_kinetic_times_sqrt_N"])),
            "frequency_variable":"q=|k_hat|^2 at the central hat-metric slice; physical squared momentum is q/N. The rescaling is positive and finite.",
            "scope":"For every fixed real N in the punctured certified interval, the Euler-first frozen three-scalar constant determinant is negative for sufficiently large q, while its positive-real-frequency leading coefficient is positive. Continuity gives at least one real growing root of this frozen local symbol. This is not an actual time-dependent instability or a UV exclusion; no numerical interacting cutoff has been inferred.",
            "clock_negative_control_without_Euler_time_jets":-sp.Rational(1,100),
            "actual_clock_Euler_leading_determinant":sp.Rational(1199,100)}
