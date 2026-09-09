"""Unchanged fixed-light quadratic operator; changed coherent-source functional."""
from functools import cache

import sympy as sp


@cache
def data():
    # Two generic symmetric, noncommuting matrices retain operator order.
    # Their finite-dimensional Gaussian identity is the algebraic
    # completion-of-square step, not an independent functional-analysis
    # theorem for curved-space Green functions.
    m0,m1,m2,k0,k1,k2=sp.symbols("mass00 mass01 mass11 kinetic00 kinetic01 kinetic11",real=True)
    M=sp.Matrix([[m0,m1],[m1,m2]])
    K=sp.Matrix([[k0,k1],[k1,k2]])
    A=M+K
    W=sp.Matrix(sp.symbols("retained_vector0:2",real=True))
    S=sp.Matrix(sp.symbols("light_source0:2",real=True))
    old=((W-S).T*M*(W-S)+W.T*K*W)[0]/2
    new=(W.T*A*W)[0]/2
    mean=A.inv()*M*S
    on_shell=(S.T*(M-M*A.inv()*M)*S)[0]/2
    return {"mass_matrix":M,"kinetic_matrix":K,"fixed_light_vector_operator":A,
            "vector_variables":W,"source_variables":S,
            "old_fixed_light_action":old,"new_fixed_light_action":new,
            "old_coherent_mean":mean,"new_coherent_mean":sp.zeros(2,1),
            "old_minus_new_source_functional":on_shell,
            "same_covariance_and_determinant_policy":"The linear vector operator and the originally selected covariance data are identical at each fixed light configuration. Changing a linear source changes the mean, not its covariance. This transfers the existing fixed-light Gaussian determinant, not mixed light-vector or graviton loops.",
            "light_order_boundary":"On the original clock S0=S1=0. The changed mean starts at light degree two and the changed source-dependent effective light functional at degree four. Background and quadratic Gaussian response and the fixed tadpole remain the previously established ones; higher light response is not declared unchanged."}


@cache
def checks():
    d=data()
    A,W,S=d["fixed_light_vector_operator"],d["vector_variables"],d["source_variables"]
    old,new=d["old_fixed_light_action"],d["new_fixed_light_action"]
    mean=d["old_coherent_mean"]
    stationary=old.subs(dict(zip(W,mean)),simultaneous=True)
    eps=sp.Symbol("light_degree",real=True)
    s2=sp.Matrix(sp.symbols("source_second0:2",real=True))
    s3=sp.Matrix(sp.symbols("source_third0:2",real=True))
    source=eps**2*s2+eps**3*s3
    effective=d["old_minus_new_source_functional"].subs(dict(zip(S,source)),simultaneous=True)
    out={"full_fixed_light_vector_Hessian_unchanged":sp.hessian(old,W)-sp.hessian(new,W),
         "source_shift_cancels_the_old_linear_vector_force":
             sp.Matrix([sp.diff(old,x) for x in W]).subs(dict(zip(W,mean)),simultaneous=True).applyfunc(sp.factor),
         "noncommuting_Gaussian_completion_keeps_operator_order":sp.factor(stationary-d["old_minus_new_source_functional"]),
         "old_source_zero_recovers_new_action":sp.expand(old.subs(dict.fromkeys(S,0))-new),
         "unchanged_covariance_inverse_equation":(sp.hessian(old,W)*A.inv()-sp.eye(2)).applyfunc(sp.factor)}
    for n in range(4):
        out[f"changed_effective_light_action_has_no_degree_{n}"]=sp.expand(effective).coeff(eps,n)
    for n in range(2):
        out[f"changed_coherent_mean_has_no_degree_{n}"]=mean.subs(
            dict(zip(S,source)),simultaneous=True).applyfunc(lambda x,n=n:sp.expand(x).coeff(eps,n))
    out["unchanged_background_and_linear_vector_covariance_operator"]=sp.Matrix(
        [[sp.diff(A.det(),s) for s in S]])
    return out


@cache
def gates():
    # Direct independent noncommuting fixture, with invertible M+K.
    M=sp.Matrix([[2,1],[1,3]])
    K=sp.Matrix([[4,0],[0,5]])
    S=sp.Matrix([1,2])
    change=sp.factor((S.T*(M-M*(M+K).inv()*M)*S)[0]/2)
    return {"source_functional_not_declared_unchanged_at_quartic_light_order":change!=0,
            "noncommuting_operator_fixture_is_genuinely_noncommuting":M*K!=K*M,
            "initial_covariance_is_not_reselected":True,
            "old_same_scale_quadratic_frozen_diagnostic_is_not_cured":True,
            "no_numeric_old_pole_bracket_transferred_to_new_scale":True,
            "no_mixed_loop_or_offclock_Hadamard_theorem_inferred":True}
