"""Literal -2 beta edge action, both source balances and reciprocal weights.

Signature +---, EH=-G R_B/2, R_B=-6(DH+2H²). The coefficients beta
have mass dimension four. Relative to S6.5's -m4*beta normalization,
replace m4*beta there by 2*beta here. S6.16 uses this -2 normalization.
"""

from functools import cache

import sympy as sp
from p8_star.exact import coefficients, positive


@cache
def derive():
    ni, nj, ai, aj, y, c = sp.symbols("N_i N_j a_i a_j y c", positive=True)
    hi, hj = sp.symbols("H_i H_j", real=True)
    beta = sp.symbols("beta0:5", real=True)
    lag = -2*(ni*ai**3*beta[0]+beta[1]*(nj*ai**3+3*ni*ai**2*aj)
              +3*beta[2]*(nj*ai**2*aj+ni*ai*aj**2)
              +beta[3]*(3*nj*ai*aj**2+ni*aj**3)+beta[4]*nj*aj**3)
    point = {aj: y*ai, nj: c*ni}
    ri = sp.cancel((-sp.diff(lag, ni)/ai**3).subs(point))
    rj = sp.cancel((-sp.diff(lag, nj)/aj**3).subs(point))
    pi = sp.cancel((sp.diff(lag, ai)/(3*ni*ai**2)).subs(point))
    pj = sp.cancel((sp.diff(lag, aj)/(3*nj*aj**2)).subs(point))
    di_y = y*(c*hj-hi)
    ci = sp.factor(sp.diff(ri, y)*di_y+3*hi*(ri+pi))
    cj = sp.factor(sp.diff(rj, y)*di_y/c+3*hj*(rj+pj))
    poly = 2*(beta[1]+2*beta[2]*y+beta[3]*y**2)
    return {"beta": beta, "N_i": ni, "N_j": nj, "a_i": ai, "a_j": aj,
            "y": y, "c": c, "H_i": hi, "H_j": hj, "L_edge": lag,
            "P": poly, "rho_i": ri, "pressure_i": pi, "rho_j": rj, "pressure_j": pj,
            "null_i": sp.factor(ri+pi), "null_j": sp.factor(rj+pj),
            "D_i_y": di_y, "C_i": ci, "C_j": cj,
            "Bianchi_weight_i": ni**2*ai**3, "Bianchi_weight_j": nj**2*aj**3,
            "null_weight_i": ni*ai**3, "null_weight_j": nj*aj**3,
            "flux": 3*ni*nj*ai**2*poly.subs(y, aj/ai)*(aj*hj-ai*hi)}


@cache
def checks():
    d = derive()
    y, c, p, hi, hj = (d[key] for key in ("y", "c", "P", "H_i", "H_j"))
    ni, nj, ai, aj = (d[key] for key in ("N_i", "N_j", "a_i", "a_j"))
    z = sp.Symbol("z")
    generating = sp.Poly((1+z*c)*(1+z*y)**3, z)
    point = {aj: y*ai, nj: c*ni}
    reverse = {d["beta"][i]: d["beta"][4-i] for i in range(5)}
    reflected = d["L_edge"].subs(reverse, simultaneous=True).subs(
        {ni: nj, nj: ni, ai: aj, aj: ai}, simultaneous=True)
    return {
        "literal_elementary_symmetric_action": sp.expand(d["L_edge"].subs(point)+2*ni*ai**3*sum(b*generating.nth(i) for i, b in enumerate(d["beta"]))),
        "reverse_edge_beta_dictionary": sp.expand(reflected-d["L_edge"]),
        "i_null_polynomial": sp.factor(d["null_i"]-(y-c)*p),
        "j_null_polynomial": sp.factor(d["null_j"]-(c-y)*p/(c*y**3)),
        "i_undivided_Bianchi": sp.factor(d["C_i"]-3*p*c*(y*hj-hi)),
        "j_undivided_Bianchi": sp.factor(d["C_j"]+3*p*(y*hj-hi)/(c*y**3)),
        "proper_Bianchi_reciprocity_two_lapses": sp.factor(d["C_i"]+c**2*y**3*d["C_j"]),
        "null_reciprocity_one_lapse": sp.factor(d["null_i"]+c*y**3*d["null_j"]),
        "literal_oriented_flux": sp.factor(d["flux"].subs(point)-ni**2*ai**3*d["C_i"]),
    }


@cache
def vertex_checks():
    """Lapse-first EH/canonical source and Noether identity also at G=0."""
    g = sp.Symbol("G", nonnegative=True)
    a, n = sp.symbols("a N", positive=True)
    ad, add, nd, phid, phidd, v, vp = sp.symbols("ad add nd phid phidd V Vphi", real=True)
    lag = -3*g*a*ad**2/n+a**3*phid**2/(2*n)-n*a**3*v
    momentum = sp.diff(lag, ad)
    e_scale = sp.diff(lag, a)-(sp.diff(momentum, a)*ad+sp.diff(momentum, ad)*add+sp.diff(momentum, n)*nd)
    h, dh = ad/(n*a), add/(n**2*a)-ad**2/(n**2*a**2)-ad*nd/(n**3*a)
    rho, pressure = phid**2/(2*n**2)+v, phid**2/(2*n**2)-v
    rho_dot = phid*phidd/n**2-phid**2*nd/n**3+vp*phid
    scalar_euler = -n*a**3*vp-(3*a**2*ad*phid/n+a**3*phidd/n-a**3*phid*nd/n**2)
    hh, hp, rm, pm, rv, pv, dr_m, dr_v = sp.symbols("H Hprime rho_m pressure_m rho_v pressure_v Drho_m Drho_v", real=True)
    e0 = 3*g*hh**2-rm-rv
    es = g*(2*hp+3*hh**2)+pm+pv
    de0 = 6*g*hh*hp-dr_m-dr_v
    return {
        "lapse_first_EH_and_canonical_source": sp.cancel(sp.diff(lag, n)/a**3-(3*g*h**2-rho)),
        "scale_first_EH_and_canonical_source": sp.cancel(e_scale/(3*n*a**2)-(g*(2*dh+3*h**2)+pressure)),
        "canonical_source_null_nonnegative": sp.expand(rho+pressure-phid**2/n**2),
        "separate_scalar_on_shell_conservation": sp.cancel(rho_dot/n+3*h*(rho+pressure)+phid*scalar_euler/(n**2*a**3)),
        "vertex_Bianchi_including_zero_EH": sp.expand(de0+3*hh*(e0-es)+dr_m+3*hh*(rm+pm)+dr_v+3*hh*(rv+pv)),
        "zero_EH_vertex_same_identity": sp.expand((de0+3*hh*(e0-es)+dr_m+3*hh*(rm+pm)+dr_v+3*hh*(rv+pv)).subs(g, 0)),
    }


def evaluate(beta, y, c):
    d = derive()
    point = dict(zip(d["beta"], coefficients(beta), strict=True))
    point.update({d["y"]: positive(y, "ratio y"), d["c"]: positive(c, "relative lapse c")})
    return {key: sp.cancel(d[key].subs(point)) for key in
            ("P", "rho_i", "pressure_i", "rho_j", "pressure_j", "null_i", "null_j")}


def controls():
    d = derive()
    return {"wrong_one_lapse_Bianchi_weight": sp.factor(d["C_i"]+d["c"]*d["y"]**3*d["C_j"]),
            "wrong_two_lapse_null_weight": sp.factor(d["null_i"]+d["c"]**2*d["y"]**3*d["null_j"]),
            "raw_J_normalization_omission": sp.factor(d["P"]-(d["beta"][1]+2*d["beta"][2]*d["y"]+d["beta"][3]*d["y"]**2))}
