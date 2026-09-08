"""Full spatial Maxwell Legendre map and one-form momentum constraint."""
from functools import cache

import sympy as sp


@cache
def legendre():
    N, kappa = sp.symbols("N kappa", positive=True)
    # G is the inverse positive spatial hat metric; the identities hold
    # for any invertible symmetric G, while positivity uses that geometry.
    g = sp.symbols("g0:6", real=True)
    G = sp.Matrix([[g[0], g[3], g[4]], [g[3], g[1], g[5]], [g[4], g[5], g[2]]])
    E, P, DW0 = [sp.Matrix(sp.symbols(prefix+"0:3", real=True)) for prefix in ("E", "Pi", "D_W0_")]
    L = kappa*(E.T*G*E)[0]/(2*N)
    solution = N*G.inv()*P/kappa
    momentum = sp.Matrix([sp.diff(L, value) for value in E])
    actual = ((P.T*(E+DW0))[0]-L).subs(dict(zip(E, solution, strict=True)), simultaneous=True)
    target = N*(P.T*G.inv()*P)[0]/(2*kappa)+(P.T*DW0)[0]
    return {"N": N, "kappa": kappa, "G": G, "E": E, "P": P, "DW0": DW0,
            "L": L, "solution": solution, "H": target,
            "full_vector_momentum_equations": (momentum.subs(dict(zip(E, solution, strict=True)), simultaneous=True)-P).applyfunc(sp.factor),
            "full_vector_Hamiltonian": sp.factor(actual-target)}


@cache
def spatial_generator():
    xi, W, P = [sp.Matrix(sp.symbols(prefix+"0:3", real=True)) for prefix in ("xi", "W", "Pi")]
    dW = sp.Matrix(3, 3, sp.symbols("dW0:9", real=True))
    dxi = sp.Matrix(3, 3, sp.symbols("dxi0:9", real=True))
    divP = sp.Symbol("div_Pi", real=True)
    actual = sum(P[j]*(xi[i]*dW[i, j]+W[i]*dxi[j, i]) for i in range(3) for j in range(3))
    target = sum(xi[i]*P[j]*(dW[i, j]-dW[j, i]) for i in range(3) for j in range(3))-(xi.T*W)[0]*divP
    boundary = (sum(dxi[j, i]*W[i]*P[j]+xi[i]*dW[j, i]*P[j] for i in range(3) for j in range(3))
                +(xi.T*W)[0]*divP)
    return {"actual": actual, "target": target,
            "vector_spatial_diffeomorphism_generator": sp.expand(actual-target-boundary)}


@cache
def checks():
    data = legendre()
    out = {"full_vector_momentum_equations": data["full_vector_momentum_equations"],
           "full_vector_Hamiltonian": data["full_vector_Hamiltonian"],
           "vector_spatial_diffeomorphism_generator": spatial_generator()["vector_spatial_diffeomorphism_generator"]}
    # After spatial IBP, Pi^i D_i W0 becomes -W0 D_i Pi^i. Pi is an
    # independent canonical momentum here, not its velocity expression.
    W0, divP, divergence = sp.symbols("W0 div_Pi divergence_PiW0", real=True)
    before = (data["P"].T*data["DW0"])[0]
    definition = divergence-before-W0*divP
    out["temporal_Gauss_spatial_boundary"] = sp.expand(before+W0*divP-divergence+definition)
    return out
