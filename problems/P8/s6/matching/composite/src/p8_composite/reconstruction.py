"""Exact local CD via a reconstructed potential and a regular ODE box.

This is not a solution for a potential fixed beforehand. It constructs
one analytic potential on the field interval traversed by the solution.
"""

from fractions import Fraction
from functools import cache

import sympy as sp


@cache
def derive():
    y, rho = sp.symbols("reconstructed_y reconstructed_rho", positive=True)
    h = sp.Symbol("prescribed_Hbar", real=True)
    r, p = 1+y, 2*y/(1+y)**2
    x = sp.sqrt(y**2+r**3*rho/3)
    z = -sp.sqrt(y**-2+r**3*rho/(3*y**3))
    denominator = x-y*z
    ng, nf = (r*h-y*z)/denominator, (x-r*h)/denominator
    v = y*r*(x*z-h*(x+z))/denominator
    w = -3*h*(rho+p)
    return {"y": y, "rho": rho, "h": h, "r": r, "p": p, "null": rho+p,
            "X": x, "Y": z, "D": denominator, "Ng": ng, "Nf": nf,
            "yprime": v, "rhoprime": w, "varphi_prime": sp.sqrt(rho+p),
            "potential_value": (rho-p)/2}


def checks():
    data = derive()
    y, rho, h, r, p, x, z, ng, nf, v, w = (data[key] for key in
        ("y", "rho", "h", "r", "p", "X", "Y", "Ng", "Nf", "yprime", "rhoprime"))
    dx = sp.diff(x, y)*v+sp.diff(x, rho)*w
    dz = sp.diff(z, y)*v+sp.diff(z, rho)*w
    root = sp.sqrt(rho+p)
    values = {"g_Friedmann": 3*x**2-3*y**2-r**3*rho,
              "f_Friedmann": 3*z**2-3/y**2-r**3*rho/y**3,
              "g_kinematics": h-v/r-ng*x, "f_kinematics": h+v/(y*r)-nf*z,
              "physical_clock": ng+nf-1,
              "g_null": -2*dx-ng*r**3*(rho+p),
              "f_null": -2*dz-nf*r**3*(rho+p)/y**3,
              "matter_conservation": w+3*h*(rho+p),
              "pressure_branch": 2*y-r**2*p,
              "canonical_scalar_equation_with_field_clock":
                  (w+sp.diff(p, y)*v)/(2*root)+3*h*root+(w-sp.diff(p, y)*v)/(2*root),
              "canonical_potential_and_density": root**2/2+data["potential_value"]-rho,
              "canonical_potential_and_pressure": root**2/2-data["potential_value"]-p}
    initial = {y: 1, rho: sp.Rational(1, 2), h: 0}
    for key, expected in {"X": sp.sqrt(sp.Rational(7, 3)), "Y": -sp.sqrt(sp.Rational(7, 3)),
                          "Ng": sp.Rational(1, 2), "Nf": sp.Rational(1, 2),
                          "yprime": -sp.sqrt(sp.Rational(7, 3)), "rhoprime": 0,
                          "varphi_prime": 1, "potential_value": 0}.items():
        values["initial_"+key] = data[key].subs(initial)-expected
    u = sp.Symbol("u", real=True)
    scale = (1+u**2)**2
    cd_h = 4*u/(1+u**2)
    values["exact_CD_scale_Hubble"] = sp.diff(scale, u)/scale-cd_h
    values["exact_CD_bounce_acceleration"] = sp.diff(cd_h, u).subs(u, 0)-4
    return {name: sp.simplify(value) for name, value in values.items()}


def positivity_box(delta=Fraction(1, 64)):
    """Exact arithmetic for the explicit |u|<=1/64 continuation proof.

The box is y in [3/4,5/4], rho in [1/4,3/4]. The accompanying proof
justifies each monotone endpoint estimate; this is not a sampled sign test.
"""
    if isinstance(delta, bool) or not isinstance(delta, (int, Fraction)):
        raise TypeError("Window radius must be an exact int/Fraction")
    if delta <= 0:
        raise ValueError("Window radius must be positive")
    q = Fraction
    delta = q(delta)
    h_upper = 4*delta
    x2_low = q(3, 4)**2+q(7, 4)**3*q(1, 4)/3
    x2_high = q(5, 4)**2+q(9, 4)**3*q(3, 4)/3
    z2_high = q(4, 3)**2+q(7, 3)**3*q(3, 4)/3
    null_low = q(1, 4)+2*q(3, 4)/q(9, 4)**2
    d_high = q(9, 4)+q(5, 4)*q(5, 2)
    lapse_numerator_low = 1-q(9, 4)*h_upper
    v_high = q(5, 4)*q(9, 4)*(q(9, 4)*q(5, 2)+h_upper*(q(9, 4)+q(5, 2)))/2
    w_high = 3*h_upper*q(5, 4)
    margins = {"X_lower_above_one": x2_low-1,
               "X_upper_below_9_over_4": q(9, 4)**2-x2_high,
               "minusY_upper_below_5_over_2": q(5, 2)**2-z2_high,
               "null_lower_above_half": null_low-q(1, 2),
               "positive_lapse_numerator": lapse_numerator_low,
               "lapse_lower_above_one_seventh": lapse_numerator_low/d_high-q(1, 7),
               "ratio_speed_below_nine": 9-v_high,
               "density_speed_below_quarter": q(1, 4)-w_high,
               "ratio_no_exit": q(1, 4)-9*delta,
               "density_no_exit": q(1, 4)-q(1, 4)*delta}
    if any(value <= 0 for value in margins.values()):
        raise ValueError("A local-CD positivity/continuation margin failed")
    return {"delta_u": delta, "h_upper": h_upper, "y_min": q(3, 4), "y_max": q(5, 4),
            "rho_min": q(1, 4), "rho_max": q(3, 4),
            "X2_lower": x2_low, "X2_upper": x2_high, "Y2_upper": z2_high,
            "null_lower": null_low, "D_upper": d_high,
            "ratio_speed_upper": v_high, "density_speed_upper": w_high,
            "strict_positive_margins": margins}


def scale_checks():
    mass, tau = sp.symbols("M tau", positive=True)
    nbar, vbar, hbar, rbar = sp.symbols("nbar Vbar hbar rhobar", real=True)
    return {"physical_Hubble": sp.simplify((hbar/tau)**2*tau**2-hbar**2),
            "physical_scalar_kinetic": sp.simplify((mass*sp.sqrt(nbar)/tau)**2-mass**2*nbar/tau**2),
            "physical_energy_restore": sp.simplify((mass**2*nbar/(2*tau**2)+mass**2*vbar/tau**2)
                                                     /(mass**2/tau**2)-nbar/2-vbar),
            "Einstein_to_interaction_units": sp.simplify((3*mass**2*(hbar/tau)**2-mass**2*rbar/tau**2)
                                                          /(mass**2/tau**2)-3*hbar**2+rbar)}
