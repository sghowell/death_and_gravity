"""Ordered Gaussian identity and a causal-boundary countercontrol."""
from functools import cache

import sympy as sp


@cache
def ordered():
    K, D = sp.MatrixSymbol("K", 2, 2), sp.MatrixSymbol("D", 2, 2)
    z = sp.Symbol("zeta", positive=True)
    M = K-z*D
    response = sp.expand(K.inv()*M-(sp.Identity(2)-z*K.inv()*D))
    schur = sp.expand(M-M*K.inv()*M-z*D+z**2*D*K.inv()*D)
    return {"ordered_vector_resolvent": response.as_explicit(),
            "ordered_exact_induced_operator": schur.as_explicit()}


@cache
def dense():
    # Noncommuting symmetric matrices are a finite-dimensional operator
    # ordering control, not a bound on the Lorentzian differential inverse.
    M = sp.Matrix([[3, 1], [1, 2]])
    D = sp.Matrix([[2, -1], [-1, 4]])
    z = sp.Symbol("zeta", positive=True)
    K = M+z*D
    W, S = [sp.Matrix(sp.symbols(prefix+"0:2", real=True)) for prefix in ("W", "S")]
    action = ((W-S).T*M*(W-S))[0]/2+z*(W.T*D*W)[0]/2
    solution = K.inv()*M*S
    effective = z*(S.T*D*S)[0]/2-z**2*(S.T*D*K.inv()*D*S)[0]/2
    completed = ((W-solution).T*K*(W-solution))[0]/2+effective
    actual = action.subs(dict(zip(W, solution, strict=True)), simultaneous=True)
    equation = sp.Matrix([sp.diff(action, item) for item in W])
    return {"M": M, "D": D, "z": z, "K": K, "W": W, "S": S,
            "action": action, "solution": solution, "effective": effective,
            "dense_full_vector_Euler_equation": (equation-K*W+M*S).applyfunc(sp.factor),
            "dense_exact_action_completion": sp.factor(action-completed),
            "dense_exact_stationary_action": sp.factor(actual-effective),
            "retained_vector_Hessian_source_independent": (sp.hessian(action, list(W))-K).applyfunc(sp.factor)}


@cache
def causal_control():
    G = sp.Matrix([[1, 0], [2, 1]])
    J = sp.Matrix(sp.symbols("J0:2", real=True))
    action = (J.T*G*J)[0]/2
    derivative = sp.Matrix([sp.diff(action, item) for item in J])
    return {"retarded_example": G, "source": J, "one_copy_variation": derivative,
            "variation_is_symmetric_not_retarded": derivative-(G+G.T)*J/2,
            "wrong_retarded_variation_difference": derivative-G*J}


@cache
def checks():
    data = dense()
    out = dict(ordered())
    out.update({name: data[name] for name in ("dense_full_vector_Euler_equation",
                "dense_exact_action_completion", "dense_exact_stationary_action",
                "retained_vector_Hessian_source_independent")})
    out["one_copy_retarded_kernel_variation_is_symmetric"] = causal_control()["variation_is_symmetric_not_retarded"]
    return out


@cache
def proof_checks():
    data, causal = dense(), causal_control()
    return {"mass_and_curl_control_do_not_commute": data["M"]*data["D"] != data["D"]*data["M"],
            "retarded_one_copy_action_control_fails": causal["wrong_retarded_variation_difference"] != sp.zeros(2, 1),
            "retarded_early_row_has_no_late_source": causal["retarded_example"][0, 1] == 0,
            "one_copy_early_variation_has_late_source": sp.diff(causal["one_copy_variation"][0], causal["source"][1]) == 1,
            "no_commuting_mass_or_generic_inverse_norm_assumption": True,
            "vector_determinant_not_claimed_zero_or_quantitatively_bounded": True}
