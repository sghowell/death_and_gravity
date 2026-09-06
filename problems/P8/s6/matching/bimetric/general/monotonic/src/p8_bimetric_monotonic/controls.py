"""Exact backgrounds separating Z monotonicity from unweighted H claims."""

from functools import cache

import sympy as sp
from p8_bimetric_general import background as bg


def decoupled_de_sitter_checks():
    t = bg.t
    point = {bg.MG2: 1, bg.MF2: 1, bg.M4: 1, bg.a: sp.exp(t), bg.b: sp.exp(-t),
             bg.Ng: 1, bg.Nf: 1, bg.KG: 0, bg.KF: 0, bg.VG: 0, bg.VF: 0,
             **dict(zip(bg.BETAS, (3, 0, 0, 0, 3), strict=True))}
    old = bg.equations()
    checks = {key: sp.simplify(old[key].subs(point).doit()) for key in ("EL_Ng", "EL_Nf", "EL_a", "EL_b")}
    h, y = sp.Integer(1), sp.exp(-2*t)
    z = h/sp.sqrt(1+y**2)
    checks["positive_Z_derivative_in_decoupled_case"] = sp.simplify(
        sp.diff(z, t)-2*sp.exp(-4*t)/(1+sp.exp(-4*t))**sp.Rational(3, 2))
    checks["constant_physical_Hubble"] = sp.diff(h, t)
    return checks


@cache
def increasing_H_example():
    y = sp.Symbol("example_y", positive=True)  # y=1+t, so d/dt=d/dy.
    h = sp.sqrt((7*y**2+5/y)/12)
    c = y+1/h
    rho = 3*h**2+2-sp.Rational(15, 4)*y
    ng = -sp.diff(rho, y)/(3*h)
    pressure = ng-rho
    potential = rho-ng/2
    return {"y": y, "H": h, "c": c, "rho_g": rho, "n_g": ng, "p_g": pressure,
            "V_g": potential, "n_f": sp.Integer(0), "Z": h/sp.sqrt(1+y**2),
            "beta": (-2, sp.Rational(5, 4), 0, 0, sp.Rational(7, 4))}


def increasing_H_checks():
    data = increasing_H_example()
    y, h, c, rho, ng = (data[key] for key in ("y", "H", "c", "rho_g", "n_g"))
    spring = sp.Rational(5, 4)
    checks = {
        "g_Friedmann": sp.factor(3*h**2-rho+2-3*spring*y),
        "f_Friedmann": sp.factor(3*(h/y)**2-sp.Rational(7, 4)-spring/y**3),
        "g_null": sp.simplify(-2*sp.diff(h, y)-ng-(y-c)*spring),
        "f_null": sp.simplify(-2*sp.diff(h/y, y)/c-(c-y)*spring/(c*y**3)),
        "g_matter_conservation": sp.simplify(sp.diff(rho, y)+3*h*ng),
        "canonical_null_stress": sp.simplify(ng-(15-14*y+5/y**2)/(12*h)),
        "positive_null_numerator_interval_lower_bound": sp.Rational(15)-14*sp.Rational(11, 10)+5/sp.Rational(11, 10)**2-sp.Rational(2258, 605),
        "Hdot_at_centre": sp.simplify(sp.diff(h, y).subs(y, 1)-sp.Rational(3, 8)),
        "n_g_at_centre": sp.simplify(ng.subs(y, 1)-sp.Rational(1, 2)),
        "Zdot_at_centre": sp.simplify(sp.diff(data["Z"], y).subs(y, 1)+sp.sqrt(2)/16)}
    # Literal bridge to the pinned four lapse/scale expressions. A'=H A,
    # b=y A, y'=1; the positive scale itself need not be explicitly integrated.
    scale = sp.Symbol("positive_example_scale", positive=True)
    t = bg.t
    point = {bg.MG2: 1, bg.MF2: 1, bg.M4: 1, bg.a: scale, bg.b: y*scale,
             bg.Ng: 1, bg.Nf: c, bg.KG: ng, bg.VG: data["V_g"], bg.KF: 0, bg.VF: 0,
             sp.diff(bg.a, t): scale*h, sp.diff(bg.b, t): scale*(1+y*h),
             sp.diff(bg.a, t, 2): scale*(h**2+sp.diff(h, y)),
             sp.diff(bg.b, t, 2): scale*(2*h+y*h**2+y*sp.diff(h, y)),
             sp.diff(bg.Ng, t): 0, sp.diff(bg.Nf, t): sp.diff(c, y),
             **dict(zip(bg.BETAS, data["beta"], strict=True))}
    for key in ("EL_Ng", "EL_Nf", "EL_a", "EL_b"):
        checks["literal_"+key] = sp.simplify(bg.equations()[key].subs(point, simultaneous=True))
    return checks


def static_sharpness_checks():
    old = bg.equations()
    point = {bg.a: 1, bg.b: 1, bg.Ng: 1, bg.Nf: 1, bg.KG: 0, bg.KF: 0, bg.VG: 0, bg.VF: 0,
             bg.M4: 1, **dict(zip(bg.BETAS, (-3, 1, 0, 0, -1), strict=True))}
    return {key: sp.simplify(old[key].subs(point).doit()) for key in ("EL_Ng", "EL_Nf", "EL_a", "EL_b")}


def checks():
    return {"decoupled_de_sitter": decoupled_de_sitter_checks(),
            "interacting_increasing_H": increasing_H_checks(),
            "static_sharp_endpoint_control": static_sharpness_checks()}


def controls():
    t = bg.t
    return {"decoupled_Z_increases_although_H_is_constant": 2*sp.exp(-4*t)/(1+sp.exp(-4*t))**sp.Rational(3, 2),
            "interacting_H_increases_at_centre": sp.Rational(3, 8),
            "interacting_Z_decreases_at_centre": -sp.sqrt(2)/16,
            "interacting_matter_NEC_at_centre": sp.Rational(1, 2),
            "NEC_numerator_positive_on_nine_tenths_to_eleven_tenths": sp.Rational(2258, 605),
            "static_parent_attains_endpoint_error_threshold": sp.Rational(8, 5)}
