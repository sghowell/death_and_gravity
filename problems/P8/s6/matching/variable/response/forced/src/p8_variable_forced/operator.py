"""Literal full TT operator and source; no K>=1 transfer API is imported."""

from functools import cache

import sympy as sp
from p8_variable_beta import canonical


@cache
def derive():
    u = sp.Symbol("u", real=True)
    delta, kk = sp.symbols("delta K", nonnegative=True)
    c, d = 2+delta, 1+u**2
    ww = c*d**12+8
    theta = 6*u*(c*d**12-8)/(d*ww)
    omega = -24*sp.sqrt(2*c)*u*d**5/ww
    w1, w2 = c*d**12/ww, 8/ww
    cf, sw = c**2*d**8/4, sp.sqrt(8*c)*d**6/ww
    cross = kk/d**4*sw*(cf-1)
    denominator = c*d**4-2
    numerator = 16*(1-u**2)*(8/(c*d**14)+1/d**2)
    aa = kk/d**4*(w1+cf*w2)-sp.diff(theta, u)-theta**2
    cc = cross-2*sp.diff(omega, u)-2*omega*theta
    ee = cross-2*omega*theta
    ba = kk/d**4*(w2+cf*w1)+sp.diff(theta, u)-theta**2-4*omega**2
    ag, bg = 2*sp.sqrt(c)*d**3/sp.sqrt(ww), -4*sp.sqrt(2)/(d**3*sp.sqrt(ww))
    j_l, j_h = sp.sqrt(c)*d**9/sp.sqrt(ww), -2*sp.sqrt(2)*d**3/sp.sqrt(ww)
    return {"u": u, "delta": delta, "K": kk, "c": c, "d": d, "a": d**2,
            "theta": theta, "omega": omega, "w1": w1, "w2": w2,
            "D": denominator, "N": numerator, "mass": numerator/denominator,
            "A": aa, "C": cc, "E": ee, "dc": -2*omega/u, "fc": 2*omega/u,
            "b_analytic": ba, "B": numerator/denominator+ba,
            "ag": ag, "bg": bg, "jL": j_l, "jH": j_h}


def coupling(field, coordinate, e_coeff, f_coeff):
    return e_coeff*field+coordinate*f_coeff*sp.diff(field, coordinate)


def adjoint_coupling(field, coordinate, c_coeff, d_coeff):
    return c_coeff*field+coordinate*d_coeff*sp.diff(field, coordinate)


@cache
def checks():
    d, old = derive(), canonical.derive()
    u, delta, kk = (d[name] for name in ("u", "delta", "K"))
    values = {old["u"]: u, old["c"]: 2+delta, old["kbar"]: sp.sqrt(kk)}
    result = {
        "actual_A_including_K_zero": d["A"]-old["V_LL"].subs(values),
        "actual_B_including_K_zero": d["B"]-old["V_HH"].subs(values),
        "actual_algebraic_mass": d["mass"]-old["mass_squared"].subs(values),
        "literal_g_source_L": d["jL"]-d["a"]**3*d["ag"]/2,
        "literal_g_source_H": d["jH"]-d["a"]**3*d["bg"]/2,
        "full_physical_retarded_diagonal": d["ag"]*d["jL"]+d["bg"]*d["jH"]-2,
        "adjoint_derivative": d["dc"]+d["fc"],
        "adjoint_multiplication": d["C"]-d["E"]+sp.diff(u*d["fc"], u),
    }
    ll, hh = sp.Function("l")(u), sp.Function("Q")(u)
    bl = coupling(ll, u, d["E"], d["fc"])
    bsh = adjoint_coupling(hh, u, d["C"], d["dc"])
    result["full_formal_adjoint_boundary"] = hh*bl-ll*bsh-sp.diff(u*d["fc"]*ll*hh, u)
    return {name: sp.simplify(sp.factor(value)) for name, value in result.items()}
