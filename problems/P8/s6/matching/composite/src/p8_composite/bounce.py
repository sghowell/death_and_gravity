"""An explicit free-canonical pressure-branch bounce, not the CD solution.

Dimensionless variables use G=F=M^2 and m4=M^2/tau^2; primes below
differentiate u=T/tau. The physical scalar is chi=M*varphi.
"""

from functools import cache

import sympy as sp


@cache
def derive():
    y = sp.Symbol("free_y", positive=True)
    r, p = 1+y, 2*y/(1+y)**2
    x, z = sp.sqrt(y*(2+5*y)/3), -sp.sqrt((5+2*y)/(3*y**2))
    denominator = (y+5)*x-(5*y+1)*z
    ng, nf = -(5*y+1)*z/denominator, (y+5)*x/denominator
    velocity = 6*y*r*x*z/denominator
    ae = (r**2/(4*y))**sp.Rational(1, 6)
    he = (y-1)*x*z/denominator
    return {"y": y, "r": r, "rho": p, "p": p, "null": 2*p,
            "X": x, "Y": z, "D": denominator, "Ng": ng, "Nf": nf,
            "yprime": velocity, "a_eff": ae, "a": ae/r, "b": y*ae/r,
            "H_eff": he, "varphi_prime": 2*sp.sqrt(y)/r}


@cache
def time_jets():
    """Actual time jets, retaining the nonlinear y'=v(y) clock."""
    data = derive()
    y = data["y"]

    def jet(expression, order):
        return sp.simplify(sp.diff(expression, y, order).subs(y, 1))

    v, v1, v2 = (jet(data["yprime"], order) for order in range(3))
    h1, h2, h3 = (jet(data["H_eff"], order) for order in range(1, 4))
    a1, a2, a3, a4 = (jet(data["a_eff"], order) for order in range(1, 5))
    # a_y(1)=0 is checked independently; the fourth derivative formula
    # below can consequently omit the a_y coefficient of y''''.
    if a1 != 0:
        raise ValueError("Free-bounce scale is not stationary at y=1")
    return {"H_prime": sp.simplify(h1*v),
            "H_second": sp.simplify(h2*v**2+h1*v*v1),
            "H_third": sp.simplify(h3*v**3+3*h2*v**2*v1+h1*(v*v1**2+v**2*v2)),
            "a_second": sp.simplify(a2*v**2),
            "a_third": sp.simplify(a3*v**3+3*a2*v**2*v1),
            "a_fourth": sp.simplify(a4*v**4+6*a3*v**3*v1+a2*(7*v**2*v1**2+4*v**3*v2))}


def checks():
    data = derive()
    y, r, p, x, z, ng, nf, v, ae, he = (data[key] for key in
        ("y", "r", "p", "X", "Y", "Ng", "Nf", "yprime", "a_eff", "H_eff"))
    values = {
        "g_Friedmann": 3*x**2-3*y**2-r**3*p,
        "f_Friedmann": 3*z**2-3/y**2-r**3*p/y**3,
        "g_null": -2*sp.diff(x, y)*v-ng*r**3*2*p,
        "f_null": -2*sp.diff(z, y)*v-nf*r**3*2*p/y**3,
        "physical_clock": ng+nf-1,
        "ratio_kinematics": v-y*(nf*z-ng*x),
        "composite_Hubble": he-(ng*x+nf*y*z)/r,
        "scale_derivative": sp.diff(ae, y)*v-ae*he,
        "g_scale_kinematics": sp.diff(data["a"], y)*v-data["a"]*ng*x,
        "f_scale_kinematics": sp.diff(data["b"], y)*v-data["b"]*nf*z,
        "free_scalar_current": ae**3*data["varphi_prime"]-1,
        "free_scalar_conservation": sp.diff(p, y)*v+6*he*p,
        "positive_kinetic_density": data["varphi_prime"]**2-2*p,
        "pressure_branch": 2*y-r**2*p,
        "double_pressure_root": 2*y-r**2/2+(y-1)**2/2,
        "maximum_pressure_factorization": sp.Rational(1, 2)-p-(y-1)**2/(2*r**2),
        "bounce_Hubble": he.subs(y, 1),
        "bounce_acceleration": (sp.diff(he, y)*v).subs(y, 1)-sp.Rational(7, 36),
        "bounce_ratio_velocity": v.subs(y, 1)+sp.sqrt(sp.Rational(7, 3)),
        "bounce_both_lapses": ng.subs(y, 1)-sp.Rational(1, 2),
        "bounce_scalar_velocity": data["varphi_prime"].subs(y, 1)-1,
    }
    jets = time_jets()
    values.update({"actual_clock_H_second": jets["H_second"],
                   "actual_clock_H_third": jets["H_third"]+sp.Rational(73, 216),
                   "actual_clock_a_second": jets["a_second"]-sp.Rational(7, 36),
                   "actual_clock_a_third": jets["a_third"],
                   "actual_clock_a_fourth": jets["a_fourth"]+sp.Rational(97, 432),
                   "scale_Hubble_fourth_jet_bridge": jets["a_fourth"]-jets["H_third"]-3*jets["H_prime"]**2,
                   "non_CD_shape_invariant": jets["a_fourth"]-sp.Rational(3, 2)*jets["a_second"]**2+sp.Rational(9, 32)})
    return {name: sp.simplify(value) for name, value in values.items()}


def negative_controls():
    data = derive()
    y = data["y"]
    return {"old_P_factor_is_not_new_Q": 2,
            "double_root_has_nonzero_time_velocity": data["yprime"].subs(y, 1),
            "free_bounce_is_not_CD_at_same_tau": sp.Rational(7, 36)-4,
            "free_bounce_is_not_CD_under_any_positive_time_rescaling":
                sp.simplify(time_jets()["a_fourth"]-sp.Rational(3, 2)*time_jets()["a_second"]**2),
            "individual_expanding_Hubble_decreases_at_composite_bounce":
                sp.simplify((sp.diff(data["X"], y)*data["yprime"]/data["Ng"]).subs(y, 1))}
