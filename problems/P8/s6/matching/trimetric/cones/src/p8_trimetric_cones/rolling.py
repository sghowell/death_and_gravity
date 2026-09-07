"""A nonempty family of actual local rolling solutions, not a CD witness.

The full beta4 extension has B=-6q,beta4g=beta4f=-q,G=F and a free scalar.
t is proper time of the common Einstein metric r, NOT physical h time.
For each A0>1 the displayed analytic ODE has a local solution with positive
lapses/scales and canonical matter. No interval/cutoff estimate is imported.
"""

from functools import cache

import sympy as sp


@cache
def derive():
    # A=1+x with x>0 lets the exact engine certify all branch signs.
    x, q, g, eps = sp.symbols("x q G epsilon", positive=True)
    a = 1+x
    n = a/(6*a-5)
    rho = 12*q*(a-1)/(eps*a)
    h2 = 2*q*x*(x**2+3*x+3)/(3*g)
    h = sp.sqrt(h2)
    aprime = -6*h*a*(a-1)/(6*a-5)
    hprime = -6*q*a**3*(a-1)/(g*(6*a-5))
    return {"x": x, "A": a, "N": n, "q": q, "G": g, "epsilon": eps,
            "rho": rho, "p": rho, "V": sp.Integer(0), "B": -6*q, "beta4": -q,
            "H_r_squared": h2, "H_r": h, "A_prime": aprime, "H_r_prime": hprime,
            "psi_prime": n*sp.sqrt(2*rho), "D_h_psi": sp.sqrt(2*rho),
            "H_h": h/a, "a_r_log_prime": h, "physical_speed": 6*a-5,
            "physical_null_stress": 2*eps*rho}


@cache
def checks():
    d = derive()
    x, a, n, q, g, eps, rho, h, aprime, hp = (d[key] for key in
        ("x", "A", "N", "q", "G", "epsilon", "rho", "H_r", "A_prime", "H_r_prime"))
    energy_rate = sp.diff(rho, x)*aprime+6*(h+aprime/a)*rho
    # The scalar current is a_h³ D_h psi. Divide its derivative by this
    # strictly positive current, so no scalar normalization is discarded.
    scalar_log_rate = 3*(h+aprime/a)+sp.diff(rho, x)*aprime/(2*rho)
    return {
        "actual_u_lapse": sp.cancel(-6*q+6*q/a+eps*rho/2),
        "actual_u_space": sp.cancel(-6*q+2*q/n+4*q/a-eps*rho/2),
        "both_Einstein_lapse": sp.cancel(3*g*d["H_r_squared"]-2*q*a**3+2*q),
        "both_Einstein_spatial": sp.cancel(g*(2*hp+3*d["H_r_squared"])-2*q*n*a**2+2*q),
        "Friedmann_constraint_is_preserved": sp.simplify(sp.diff(d["H_r_squared"], x)*aprime-2*h*hp),
        "actual_free_scalar_energy_equation": sp.simplify(energy_rate),
        "actual_free_scalar_current_equation": sp.simplify(scalar_log_rate),
        "actual_physical_Hubble_conversion": sp.simplify((h+aprime/a)/n-d["H_h"]),
        "actual_tensor_speed_gap": sp.cancel(a/n-1-a*d["physical_null_stress"]/(4*q)),
    }


def initial_data():
    d = derive()
    return {key: sp.simplify(d[key].subs(d["x"], 1)) for key in
            ("A", "N", "rho", "H_r_squared", "H_r_prime", "D_h_psi", "physical_speed")}
