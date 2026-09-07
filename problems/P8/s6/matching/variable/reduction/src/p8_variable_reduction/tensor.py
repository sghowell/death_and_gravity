"""Literal physical-metric tensor test of the formal four-derivative action.

No stationary-reduction or coefficient-dictionary module is imported.
The homogeneous anisotropy is an off-shell test, not a FLRW solution.
"""
from functools import cache

import sympy as sp

t = sp.Symbol("t", real=True)
q, phi = (sp.Function(name)(t) for name in ("q", "phi"))


@cache
def geometry():
    metric = sp.diag(1, -sp.exp(q), -sp.exp(-q), -1)
    inverse = metric.inv()

    def partial(value, coordinate):
        return sp.diff(value, t) if coordinate == 0 else sp.S.Zero

    gamma = [[[sp.simplify(sum(
        inverse[a, d]*(partial(metric[d, c], b)+partial(metric[d, b], c)
                       -partial(metric[b, c], d))/2 for d in range(4)))
        for c in range(4)] for b in range(4)] for a in range(4)]
    ricci = sp.Matrix(4, 4, lambda b, c: sp.simplify(sum(
        partial(gamma[a][c][b], a)-partial(gamma[a][a][b], c)
        +sum(gamma[a][a][d]*gamma[d][c][b]-gamma[a][c][d]*gamma[d][a][b]
             for d in range(4)) for a in range(4))))
    scalar = sp.simplify(sp.trace(inverse*ricci))
    ricci_square = sp.simplify(sp.trace(inverse*ricci*inverse*ricci))
    gradient = sp.Matrix([sp.diff(phi, t), 0, 0, 0])
    hessian = sp.Matrix(4, 4, lambda a, b: sp.simplify(
        partial(gradient[b], a)-sum(gamma[c][a][b]*gradient[c] for c in range(4))))
    up = inverse*gradient
    box = sp.simplify(sp.trace(inverse*hessian))
    sandwich = (up.T*hessian*up)[0]
    invariants = (
        sp.trace(inverse*hessian*inverse*hessian), box**2, sandwich*box,
        (up.T*hessian*inverse*hessian*up)[0], sandwich**2)
    return {"metric": metric, "inverse": inverse, "ricci": ricci,
            "scalar": scalar, "ricci_square": ricci_square,
            "curvature_combination": sp.simplify(ricci_square-scalar**2/3),
            "hessian": hessian, "X": (gradient.T*inverse*gradient)[0],
            "box": box, "Li": tuple(map(sp.simplify, invariants))}


def euler(density):
    return sp.expand(sp.diff(density, q)-sp.diff(sp.diff(density, sp.diff(q, t)), t)
                     +sp.diff(sp.diff(density, sp.diff(q, t, 2)), t, 2))


@cache
def derive():
    g = geometry()
    F, K, F2, *Ai = (sp.Function(name)(t) for name in
                     ("F", "K", "F2", "A1", "A2", "A3", "A4", "A5"))
    kappa = sp.Function("kappa")(t)
    density = F+K*g["box"]+F2*g["scalar"]+sum(a*l for a, l in zip(Ai, g["Li"]))
    anisotropy = sp.expand(density-density.subs({sp.diff(q, t): 0, sp.diff(q, t, 2): 0}))
    eps = sp.Symbol("epsilon", real=True)
    curvature_quadratic = sp.expand(g["curvature_combination"].subs(q, eps*q).doit()).coeff(eps, 2)
    full_quadratic = anisotropy+kappa*curvature_quadratic
    equation = euler(full_quadratic)
    return {"F2": F2, "A1": Ai[0], "kappa": kappa, "density": density,
            "DHOST_tensor_density": anisotropy, "GT": -2*F2+2*g["X"]*Ai[0],
            "curvature_quadratic": curvature_quadratic,
            "full_quadratic": full_quadratic, "equation": equation,
            "fourth_symbol": equation.coeff(sp.diff(q, t, 4))}


@cache
def checks():
    g, d = geometry(), derive()
    qt, qtt = sp.diff(q, t), sp.diff(q, t, 2)
    pt, ptt = sp.diff(phi, t), sp.diff(phi, t, 2)
    expected_ricci = sp.diag(-qt**2/2, sp.exp(q)*qtt/2, -sp.exp(-q)*qtt/2, 0)
    expected_Li = (ptt**2+pt**2*qt**2/2, ptt**2, pt**2*ptt**2,
                   pt**2*ptt**2, pt**4*ptt**2)
    residuals = {f"literal_Ricci_{i}{j}": sp.simplify(g["ricci"][i, j]-expected_ricci[i, j])
                 for i in range(4) for j in range(4)}
    residuals.update({f"literal_clock_L{i+1}": sp.simplify(actual-expected)
                      for i, (actual, expected) in enumerate(zip(g["Li"], expected_Li))})
    residuals.update({
        "unit_volume": sp.simplify(g["metric"].det()+1),
        "literal_R_B": sp.simplify(g["scalar"]+qt**2/2),
        "literal_curvature_square": sp.simplify(g["curvature_combination"]-qtt**2/2-qt**4/6),
        "DHOST_tensor_kinetic": sp.simplify(d["DHOST_tensor_density"]-d["GT"]*qt**2/4),
        "curvature_tensor_second_variation": sp.simplify(d["curvature_quadratic"]-qtt**2/2),
        "DHOST_no_fourth_symbol": euler(d["DHOST_tensor_density"]).coeff(sp.diff(q, t, 4)),
        "nonzero_retained_fourth_symbol": sp.simplify(d["fourth_symbol"]-d["kappa"]),
        "moving_kappa_third_symbol": sp.simplify(d["equation"].coeff(sp.diff(q, t, 3))
                                                 -2*sp.diff(d["kappa"], t)),
    })
    if any(value != 0 for value in residuals.values()):
        raise ValueError("A literal physical-metric tensor identity failed")
    return residuals


def calibration():
    M, tau, c = sp.symbols("M tau c", positive=True)
    beta1_center = 32*M**2/(tau**2*c*(c-2))
    kappa_center = sp.cancel(M**4*2**3/(4*beta1_center))
    return {"pure_physical_TT_fourth_symbol": "kappa",
            "kappa_center": kappa_center,
            "strict_domain": "M,tau>0; 2<c<=4; no claim at c=2",
            "symbol_scope": "formal retained action; not parent ghost or physical cutoff",
            "DHOST_result": "second order for every F,K,F2,A1,...,A5 on a homogeneous timelike clock"}
