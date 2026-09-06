"""Full ODE-derived CD-centre first/second jets at two specified matter states."""

from functools import cache

import sympy as sp
from p8_composite import reconstruction

A = sp.Symbol("A", nonnegative=True)
RADICAL = sp.sqrt(442)


@cache
def initial_jets(initial_y):
    """u=mT, G=F=M²,m4=M²m²,rho0=1/2, h0=h''0=0, h'0=A.

    Only the exact initial ratios 1 and 2 are certified by this API.
    In CD, A=4/(m*tau)². Setting h''0=0 is part of this jet contract.
    """
    if type(initial_y) is not int or initial_y not in (1, 2):
        raise ValueError("Only exact specified initial ratios 1 and 2")
    old = reconstruction.derive()
    y, rho, h = (old[key] for key in ("y", "rho", "h"))
    coordinates = (y, rho, h)
    initial = {y: sp.S(initial_y), rho: sp.Rational(1, 2), h: 0}
    c_expr = old["Nf"]/old["Ng"]
    c = sp.simplify(c_expr.subs(initial))
    yp = sp.simplify(old["yprime"].subs(initial))
    velocity = (yp, sp.S.Zero, A)
    ypp = sp.simplify(sum(sp.diff(old["yprime"], z).subs(initial)*v for z, v in zip(coordinates, velocity, strict=True)))
    rpp = -3*A*(sp.Rational(1, 2)+2*initial_y/sp.S(initial_y+1)**2)
    acceleration = (ypp, rpp, sp.S.Zero)
    cp = sp.simplify(sum(sp.diff(c_expr, z).subs(initial)*v for z, v in zip(coordinates, velocity, strict=True)))
    cpp = sum(sp.diff(c_expr, z).subs(initial)*v for z, v in zip(coordinates, acceleration, strict=True))
    cpp += sum(sp.diff(c_expr, z, zz).subs(initial)*v*vv
               for z, v in zip(coordinates, velocity, strict=True)
               for zz, vv in zip(coordinates, velocity, strict=True))
    cpp = sp.simplify(cpp)
    yy, cc = sp.symbols("ratio_y ratio_c", positive=True)
    at = {yy: sp.S(initial_y), cc: c}
    # ln(f²) excluding its irrelevant constant, but retaining 3 ln(Ae).
    logs = {"sum": -3*sp.log(1+yy)+sp.log(1+cc)+sp.log(cc+yy**3)-sp.log(cc),
            "relative": -3*sp.log(1+yy)+sp.log(1+cc)+3*sp.log(yy)-sp.log(cc+yy**3)}
    result = {"A": A, "initial_y": initial_y, "c": c, "y_prime": yp,
              "y_second": ypp, "rho_second": rpp, "c_prime": cp, "c_second": cpp}
    for label, logarithm in logs.items():
        grad = [sp.diff(logarithm, z).subs(at) for z in (yy, cc)]
        theta = sp.simplify(sum(a*b for a, b in zip(grad, (yp, cp), strict=True))/2)
        theta_prime = 3*A/2+sum(a*b for a, b in zip(grad, (ypp, cpp), strict=True))/2
        theta_prime += sum(sp.diff(logarithm, z, zz).subs(at)*v*vv/2
                           for z, v in zip((yy, cc), (yp, cp), strict=True)
                           for zz, vv in zip((yy, cc), (yp, cp), strict=True))
        result["theta_"+label] = theta
        result["N_"+label] = sp.simplify(theta_prime+theta**2)
    w1 = cc/(cc+yy**3)
    sc = sp.sqrt(w1*(1-w1))
    difference = (3*yp/initial_y-cp/c)/2
    omega = sp.simplify(sc.subs(at)*difference)
    scprime = sp.diff(sc, yy).subs(at)*yp+sp.diff(sc, cc).subs(at)*cp
    difference_prime = (3*(ypp/initial_y-yp**2/initial_y**2)-cpp/c+cp**2/c**2)/2
    result["omega"] = omega
    result["omega_prime"] = sp.simplify(scprime*difference+sc.subs(at)*difference_prime)
    mu = sp.S(initial_y)*(initial_y-1)*(initial_y-c)/(initial_y+1)
    mass = sp.simplify(mu/(1+c)**2*(1+c/initial_y**3))
    cg2, cf2 = ((1+yy)/(1+cc))**2, (cc*(1+yy)/(yy*(1+cc)))**2
    result.update({"mu": mu, "mass_squared": mass,
                   "c_light_squared": sp.simplify((w1*cg2+(1-w1)*cf2).subs(at)),
                   "c_heavy_squared": sp.simplify(((1-w1)*cg2+w1*cf2).subs(at)),
                   "D": sp.simplify((sc*(cf2-cg2)).subs(at)),
                   "Routh_frequency_squared": sp.simplify(mass-result["N_relative"]),
                   "heavy_diagonal_at_zero_k": sp.simplify(mass-result["N_relative"]-4*omega**2)})
    return result


def checks():
    symmetric, asymmetric = initial_jets(1), initial_jets(2)
    root = sp.sqrt(sp.Rational(7, 3))
    n = (144*A**2-6*A+67)/84
    om = (3*A-4)/(3*root)
    d = RADICAL
    mass2 = (107-5*d)/42
    omega2 = ((376+16*d)*A-1037)**2/(16*(60996+3047*d))
    return {"symmetric_actual_ratio_clock": sp.simplify(symmetric["y_prime"]+root),
            "symmetric_actual_lapse_clock": sp.simplify(symmetric["c_prime"]+(5+12*A)/(3*root)),
            "symmetric_ratio_second_not_omitted": sp.simplify(symmetric["y_second"]-symmetric["y_prime"]**2),
            "symmetric_lapse_second_not_omitted": sp.simplify(symmetric["c_second"]-symmetric["c_prime"]**2),
            "symmetric_theta_sum_zero": symmetric["theta_sum"],
            "symmetric_theta_relative_zero": symmetric["theta_relative"],
            "symmetric_omega_prime_zero": symmetric["omega_prime"],
            "symmetric_omega": sp.simplify(symmetric["omega"]-om),
            "symmetric_N_sum": sp.simplify(symmetric["N_sum"]-n),
            "symmetric_N_relative": sp.simplify(symmetric["N_relative"]-n+4*om**2),
            "symmetric_zero_algebraic_mass": symmetric["mass_squared"],
            "symmetric_locked_luminal": symmetric["c_light_squared"]-1,
            "symmetric_heavy_diagonal_equals_light": sp.simplify(symmetric["heavy_diagonal_at_zero_k"]+n),
            "symmetric_positive_Routh_countercontrol": symmetric["Routh_frequency_squared"].subs(A, sp.Rational(4, 25))-sp.Rational(153, 100),
            "symmetric_Routh_not_parametrically_above_locked_normalization": sp.simplify(3*n-symmetric["Routh_frequency_squared"]-(36*A**2+30*A+1)/7),
            "symmetric_positive_N_completion": sp.simplify(n-sp.Rational(12, 7)*(A-sp.Rational(1, 48))**2-sp.Rational(51, 64)),
            "asymmetric_positive_mass_formula": sp.simplify(asymmetric["mass_squared"]-mass2),
            "asymmetric_moving_eigenvector_formula": sp.simplify(asymmetric["omega"]**2-omega2),
            "asymmetric_locked_cone_formula": sp.simplify(asymmetric["c_light_squared"]-(46631-1928*d)/6517)}
