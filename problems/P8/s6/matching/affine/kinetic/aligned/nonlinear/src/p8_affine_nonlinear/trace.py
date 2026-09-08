"""Exact coupled trace/temporal-vector Legendre transformation."""
from functools import cache

import sympy as sp
from p8_affine_aligned import alignment
from p8_affine_retuned import geometry

from . import adm


@cache
def legendre():
    a = sp.Symbol("a_trace", nonzero=True)
    volume, gamma = sp.symbols("volume gamma_t", positive=True)
    b, f, T, delta, c, K, p, j = sp.symbols("b f T delta c K p j", real=True)
    L = a*K**2+b*K+f+volume*(T-delta*K-c)**2/(2*gamma)
    a_total = a+volume*delta**2/(2*gamma)
    Ksolution = (p-b+volume*delta*(T-c)/gamma)/(2*a_total)
    H = sp.factor((p*K-L).subs(K, Ksolution)-T*j)
    Kjoint = (p-b-delta*j)/(2*a)
    Tjoint = delta*Kjoint+c-gamma*j/volume
    after = (p-b-delta*j)**2/(4*a)-f+gamma*j**2/(2*volume)-c*j
    gamma_eff = gamma+volume*delta**2/(2*a)
    return {"a": a, "volume": volume, "gamma": gamma, "b": b, "f": f,
            "T": T, "delta": delta, "c": c, "K": K, "p": p, "j": j,
            "L": L, "K_solution": Ksolution, "H_before_temporal": H,
            "K_joint": Kjoint, "T_joint": Tjoint, "H_after_temporal": after,
            "gamma_effective": gamma_eff,
            "metric_trace_momentum_Euler": sp.factor(sp.diff(L, K).subs(K, Ksolution)-p),
            "temporal_secondary_equation": sp.factor(sp.diff(H, T).subs(T, Tjoint)),
            "joint_trace_equation": sp.factor(sp.diff(L, K).subs({K: Kjoint, T: Tjoint}, simultaneous=True)-p),
            "full_temporal_elimination": sp.factor(H.subs(T, Tjoint)-after),
            "temporal_secondary_Jacobian": sp.factor(sp.diff(H, T, 2)+volume/gamma_eff)}


@cache
def clock_auxiliary():
    data = legendre()
    old = adm.clock_lapse()
    bg = alignment.parent.old.background()
    h, H = bg["h"], bg["H"]
    N, s = adm.N, adm.s
    transformed = adm.transform()
    volume = transformed["volume"].subs(s, 1/N)
    B = transformed["B"].subs(s, 1/N)
    I = old["I_N"]*(N-1)+old["I_NN"]*(N-1)**2/2
    linear = volume*transformed["Blinear"].subs(s, 1/N)-I
    potential = volume*transformed["Fnew"].subs(s, 1/N)-sp.diff(I, adm.u)/N
    p_affine = sp.sqrt((h-1+N**-2)/(4*h))
    gamma = geometry.update()["gamma_t"].subs(geometry.P, p_affine)
    delta = (N**-2-1)/h
    # Q=Q_N=0 on the clock; Q_NN=-3h'/h³ follows from the original ODE.
    Q = -3*sp.diff(h, adm.u)*(N-1)**2/(2*h**3)
    c = -3*H*delta/N+sp.Rational(3, 2)*Q/N
    substitutions = {data["a"]: sp.Rational(2, 3)*volume*B,
                     data["b"]: linear, data["f"]: potential,
                     data["gamma"]: gamma, data["volume"]: volume,
                     data["delta"]: delta, data["c"]: c, data["p"]: -2*H, data["j"]: 0}
    # Work only to the needed Taylor order before the Hessian; no higher
    # coefficient jet is inferred from this background rank calculation.
    substitution_jets = {key: sp.diff(value, N, 0).subs(N, 1)
                         +sp.diff(value, N).subs(N, 1)*(N-1)
                         +sp.diff(value, N, 2).subs(N, 1)*(N-1)**2/2
                         for key, value in substitutions.items()}
    Hbefore = N*(data["H_before_temporal"].subs(substitution_jets, simultaneous=True)+bg["ell"]**2/(2*volume))
    T = data["T"]
    jacobian = sp.hessian(Hbefore, (N, T)).subs({N: 1, T: 0}).applyfunc(sp.factor)
    expected = sp.diag(-2*bg["J"], -1)
    return {"jacobian": jacobian, "expected": expected,
            "actual_joint_auxiliary_Jacobian": (jacobian-expected).applyfunc(sp.factor)}


@cache
def checks():
    data = legendre()
    names = ("metric_trace_momentum_Euler", "temporal_secondary_equation", "joint_trace_equation",
             "full_temporal_elimination", "temporal_secondary_Jacobian")
    out = {name: data[name] for name in names}
    out["actual_joint_auxiliary_Jacobian"] = clock_auxiliary()["actual_joint_auxiliary_Jacobian"]
    return out
