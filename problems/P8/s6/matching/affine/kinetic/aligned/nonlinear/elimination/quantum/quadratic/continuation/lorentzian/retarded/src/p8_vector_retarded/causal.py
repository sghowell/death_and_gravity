"""Kubo/covariance equivalence and background-cone scope controls."""
from functools import cache

import sympy as sp


@cache
def covariance():
    F = sp.Matrix(2, 2, sp.symbols("F0:4", real=True))
    s0, s1, s2 = sp.symbols("covariance_00 covariance_01 covariance_11", real=True)
    Sigma = sp.Matrix([[s0, s1], [s1, s2]])
    a0, a1, a2, b0, b1, b2 = sp.symbols("output0:3 input0:3", real=True)
    A, B = sp.Matrix([[a0, a1], [a1, a2]]), sp.Matrix([[b0, b1], [b1, b2]])
    J = sp.Matrix([[0, 1], [-1, 0]])
    W = F*(Sigma+sp.I*J/2)
    connected = sp.trace(A*W*B*W.T)/2
    kubo = sp.I*(connected-sp.conjugate(connected))
    deltaM = -J*B
    deltaSigma = F*(deltaM*Sigma+Sigma*deltaM.T)*F.T
    causal_response = sp.trace(A*deltaSigma)/2
    return {"connected_Wick_covariance": sp.expand(connected),
            "source_generator": deltaM,
            "retarded_covariance_integrand": sp.expand(causal_response),
            "Kubo_equals_exact_canonical_covariance": sp.expand(kubo-causal_response),
            "real_response": sp.expand(sp.im(kubo))}


@cache
def checks():
    data = covariance()
    out = {name: data[name] for name in ("Kubo_equals_exact_canonical_covariance", "real_response")}
    H1, H2, E1, E2 = sp.symbols("symmetric1 symmetric2 commutator1 commutator2", real=True)
    difference = (H1+sp.I*E1)*(H2+sp.I*E2)-(H1-sp.I*E1)*(H2-sp.I*E2)
    out["commutator_factorization"] = sp.expand(difference-2*sp.I*(H1*E2+E1*H2))
    out["spacelike_current_commutator_zero"] = sp.expand(difference.subs({E1: 0, E2: 0}))
    # Every finite Taylor coefficient of a moving delta has support at the old front.
    e, time, position = sp.symbols("perturbation time position", real=True)
    shifted = sp.DiracDelta(position-(1+e)*time)
    for derivative in range(3):
        target = (-time)**derivative*sp.DiracDelta(position-time, derivative)
        out["moving_front_Taylor_control_"+str(derivative)] = (
            sp.diff(shifted, e, derivative).subs(e, 0)-target)
    return out
