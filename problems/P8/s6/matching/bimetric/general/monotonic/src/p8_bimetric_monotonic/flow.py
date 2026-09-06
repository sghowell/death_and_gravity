"""Lapse-aware all-time identities and separate algebraic-root formulas."""

from functools import cache

import sympy as sp
from p8_bimetric_general import background as bg

t = bg.t


@cache
def derive():
    old = bg.equations()
    y, c, h, hf = (old[key] for key in ("y", "c", "H_g", "H_f"))
    dg = lambda value: sp.diff(value, t)/bg.Ng
    df = lambda value: sp.diff(value, t)/bg.Nf
    inertia = bg.MG2+bg.MF2*y**2
    z = h/sp.sqrt(inertia)
    null_g, null_f = old["rho_g"]+old["p_g"], old["rho_f"]+old["p_f"]
    branch = hf-h/y
    full = -2*bg.MG2*dg(h)-2*bg.MF2*c*y**3*df(hf)
    reduced = -2*inertia*dg(h)+2*bg.MF2*y*dg(y)*h
    euler = old["EL_a"]-old["EL_Ng"]+c*y**3*(old["EL_b"]-old["EL_Nf"])
    null = null_g+c*y**3*null_f
    numerator = bg.Ng*sp.diff(bg.b, t)-bg.Nf*sp.diff(bg.a, t)
    return {"A": inertia, "Z": z, "P": bg.interaction_polynomial(y), "y": y, "c": c,
            "H_g": h, "H_f": hf, "n_g": null_g, "n_f": null_f, "weighted_matter_null": null,
            "Bianchi_numerator": numerator, "dynamic_branch_error": branch,
            "full_weighted_null": full, "dynamic_weighted_null": reduced,
            "weighted_Euler_residual": euler, "Dg_Z": dg(z),
            "off_branch_Z_identity_residual": euler-2*bg.MF2*c*y**3*df(branch)-2*inertia**sp.Rational(3, 2)*dg(z)-null}


def identities():
    data = derive()
    df = lambda value: sp.diff(value, t)/bg.Nf
    dg = lambda value: sp.diff(value, t)/bg.Ng
    c, y = data["c"], data["y"]
    return {
        "full_parent_null_balance": sp.factor(data["weighted_Euler_residual"]+data["full_weighted_null"]-data["weighted_matter_null"]),
        "off_branch_differentiated_constraint": sp.factor(data["full_weighted_null"]-data["dynamic_weighted_null"]
                                                          +2*bg.MF2*c*y**3*df(data["dynamic_branch_error"])),
        "undivided_Bianchi_to_dynamic_error": sp.factor(data["dynamic_branch_error"]-data["Bianchi_numerator"]/(bg.Ng*bg.Nf*bg.b)),
        "exact_normalized_Hubble_derivative": sp.simplify(data["dynamic_weighted_null"]+2*data["A"]**sp.Rational(3, 2)*data["Dg_Z"]),
        "lapse_retaining_off_branch_Z_identity": sp.simplify(data["off_branch_Z_identity_residual"]),
        "physical_volume_ratio_weight": sp.factor(c*y**3-bg.Nf*bg.b**3/(bg.Ng*bg.a**3)),
        "dynamic_Hubble_ratio_derivative": sp.factor(df(data["H_g"]/y)-dg(data["H_g"])/(c*y)+data["H_g"]*dg(y)/(c*y**2)),
    }


def root_checks():
    h, y, lapse = (sp.Function(name)(t) for name in ("H", "y", "Ng"))
    y0 = sp.Symbol("positive_root_y0", positive=True)
    ng = sp.Symbol("n_g", nonnegative=True)
    inertia = bg.MG2+bg.MF2*y**2
    zprime = sp.diff(h/sp.sqrt(inertia), t)/lapse
    root_value = zprime.subs({y: y0, sp.diff(y, t): 0,
                             sp.diff(h, t): -lapse*ng/(2*bg.MG2)}, simultaneous=True)
    zeros = dict.fromkeys(bg.BETAS[1:4], 0)
    old = bg.equations()
    return {"constant_ratio_root_interval_formula": sp.simplify(root_value+ng/(2*bg.MG2*sp.sqrt(bg.MG2+bg.MF2*y0**2))),
            "zero_polynomial_g_interaction_null": sp.factor((old["rho_g_interaction"]+old["p_g_interaction"]).subs(zeros)),
            "zero_polynomial_f_interaction_null": sp.factor((old["rho_f_interaction"]+old["p_f_interaction"]).subs(zeros)),
            "zero_polynomial_g_Hubble_equation": sp.factor(
                (old["EL_a"]-old["EL_Ng"]).subs(zeros)
                -2*bg.MG2*sp.diff(old["H_g"], t)/bg.Ng-old["rho_g"]-old["p_g"])}


def endpoint_checks():
    tau = sp.Symbol("tau", positive=True)
    h = 4*t/(tau**2+t**2)
    error = sp.Rational(8, 5)/tau
    return {"left_CD_endpoint": sp.factor(h.subs(t, -tau/2)+error),
            "right_CD_endpoint": sp.factor(h.subs(t, tau/2)-error),
            "static_parent_saturates_left_error": sp.factor(-h.subs(t, -tau/2)-error),
            "static_parent_saturates_right_error": sp.factor(h.subs(t, tau/2)-error)}
