"""Literal norm-two diagonal-exponential TT action, including the matter source."""

from functools import cache

import sympy as sp

from . import model as m


@cache
def derive():
    d, z = sp.symbols("d z", real=True)
    # d=gamma-h, D=diag(1,-1,0), tr(D^2)=2. Both metric
    # spatial determinants are exactly independent of their TT amplitudes.
    generating = (1+m.c*z)*(1+m.y*z)*(1+2*m.y*sp.cosh(d/2)*z+m.y**2*z**2)
    e = sp.Poly(generating, z)
    interaction = -m.M4*m.Ng*m.a**3*sum(m.betas[j]*e.nth(j) for j in range(5))
    spatial_volume = m.a**3*m.r*(m.alpha**2+m.beta**2*m.y**2
                                  +2*m.alpha*m.beta*m.y*sp.cosh(d/2))
    matter = m.Ng*m.s*m.p*spatial_volume
    literal = sp.diff(interaction+matter, d, 2).subs(d, 0)/2
    mu = m.y*(m.M4*(m.betas[1]+m.betas[2]*(m.y+m.c)+m.betas[3]*m.c*m.y)
              -m.alpha*m.beta*m.r*m.s*m.p)
    hg, hf, dh, dg = sp.symbols("Hg Hf hdot gammadot", real=True)
    kg = sp.diag(hg+dh/(2*m.Ng), hg-dh/(2*m.Ng), hg)
    kf = sp.diag(hf+dg/(2*m.Nf), hf-dg/(2*m.Nf), hf)
    lg = m.G*m.Ng*m.a**3*(sp.trace(kg*kg)-sp.trace(kg)**2)/2
    lf = m.F*m.Nf*m.b**3*(sp.trace(kf*kf)-sp.trace(kf)**2)/2
    return {"d": d, "generating": generating, "interaction": interaction,
            "matter": matter, "literal_quadratic": literal, "mu": mu,
            "literal_g_inertia": sp.diff(lg, dh, 2)/2,
            "literal_f_inertia": sp.diff(lf, dg, 2)/2,
            # Coefficients per TT tensor contraction, rather than norm-two amplitude.
            "K_g": m.G*m.a**3/(8*m.Ng), "K_f": m.F*m.b**3/(8*m.Nf),
            "U_g_gradient": m.G*m.Ng*m.a/8, "U_f_gradient": m.F*m.Nf*m.b/8,
            "U_relative": m.Ng*m.a**3*mu/8}


@cache
def spatial_curvature():
    """Literal positive spatial metric diag(exp(h(z)), exp(-h(z)), 1)."""
    z = sp.Symbol("z", real=True)
    h = sp.Function("h")(z)
    metric = sp.diag(sp.exp(h), sp.exp(-h), 1)
    inverse = metric.inv()

    def derivative(value, index):
        return sp.diff(value, z) if index == 2 else sp.S.Zero

    connection = [[[sum(inverse[i, ell]*(derivative(metric[ell, k], j)
                    +derivative(metric[ell, j], k)-derivative(metric[j, k], ell))
                    for ell in range(3))/2 for k in range(3)]
                   for j in range(3)] for i in range(3)]
    ricci = sp.zeros(3)
    for i in range(3):
        for j in range(3):
            ricci[i, j] = sum(derivative(connection[k][i][j], k)
                              -derivative(connection[k][i][k], j)
                              +sum(connection[k][k][ell]*connection[ell][i][j]
                                   -connection[k][j][ell]*connection[ell][i][k]
                                   for ell in range(3)) for k in range(3))
    scalar = sp.simplify(sp.trace(inverse*ricci))
    return scalar, sp.diff(h, z)


def checks():
    d = derive()
    scalar, hz = spatial_curvature()
    pressure_branch = {m.betas[0]: 0, m.betas[1]: 0, m.betas[2]: 1,
                       m.betas[3]: 0, m.betas[4]: 0, m.alpha: 1, m.beta: 1,
                       m.p: 2*m.M4*m.y/(1+m.y)**2}
    return {"literal_source_aware_TT_stiffness": m.cancel(d["literal_quadratic"]+2*d["U_relative"]),
            "g_extrinsic_curvature_kinetic": m.cancel(d["literal_g_inertia"]-2*d["K_g"]),
            "f_extrinsic_curvature_kinetic": m.cancel(d["literal_f_inertia"]-2*d["K_f"]),
            "literal_spatial_gradient": sp.simplify(scalar+hz**2/2),
            "pressure_branch_factorization": m.cancel(d["mu"].subs(pressure_branch, simultaneous=True)
                                                        -m.M4*m.y*(m.y-1)*(m.y-m.c)/(1+m.y))}


def negative_controls():
    d = derive()
    omitted = sp.diff(d["interaction"], d["d"], 2).subs(d["d"], 0)/2
    return {"omitted_matter_pressure_changes_TT": m.cancel(omitted-d["literal_quadratic"])}
