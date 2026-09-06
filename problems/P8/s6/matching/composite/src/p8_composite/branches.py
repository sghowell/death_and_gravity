"""Undivided branch and nondegenerate-intersection identities.

The shared-matter pressure factor depends on time through p: it is not
the fixed polynomial root set in the separate-matter S6.5 theorem.
"""

import sympy as sp

from . import background as bg


def checks():
    data = bg.equations()
    y, c, r, s = (data[key] for key in ("y", "c", "r", "s"))
    error = data["H_f"]-data["H_g"]/y
    inertia = bg.G+bg.F*y**2
    z = data["H_g"]/sp.sqrt(inertia)
    full = -2*bg.G*sp.diff(data["H_g"], bg.t)/bg.Ng-2*bg.F*c*y**3*sp.diff(data["H_f"], bg.t)/bg.Nf
    answer = {
        "undivided_kinematic_error": sp.factor(error-data["B"]/(bg.Ng*bg.Nf*bg.b)),
        "effective_Hubble_off_branch": sp.factor(data["H_eff"]-data["H_g"]/r-bg.BETA*c*y*error/(s*r)),
        "weighted_normalized_Hubble_off_branch": sp.simplify(
            full+2*inertia**sp.Rational(3, 2)*sp.diff(z, bg.t)/bg.Ng
            +2*bg.F*c*y**3*sp.diff(error, bg.t)/bg.Nf),
    }
    n = sp.Symbol("nonnegative_scalar_null", nonnegative=True)
    dg_h = -s*r**3*n/(2*inertia)
    answer["dynamic_bounce_effective_acceleration"] = sp.factor(dg_h/(s*r)+r**2*n/(2*inertia))
    # At B=H_eff=0, both constituent velocities vanish. Retain both
    # lapse factors before imposing the Q=0 acceleration equations.
    acceleration = sp.diff(data["H_eff"], bg.t)/data["N_eff"]
    substitutions = {sp.diff(bg.a, bg.t): 0, sp.diff(bg.b, bg.t): 0,
                     sp.diff(bg.a, bg.t, 2): -bg.a*bg.Ng**2*bg.ALPHA*r**3*n/(2*bg.G),
                     sp.diff(bg.b, bg.t, 2): -bg.b*bg.Nf**2*bg.BETA*r**3*n/(2*bg.F*y**3)}
    intersection = -r**2*n*(bg.ALPHA**2/bg.G+bg.BETA**2*c**2/(bg.F*y**2))/(2*s**2)
    answer["Q_B_zero_intersection_acceleration"] = sp.factor(acceleration.subs(substitutions, simultaneous=True)-intersection)
    yy, pp = sp.symbols("y pressure", real=True)
    aa, bb = bg.ALPHA, bg.BETA
    polynomial = bg.M4*(bg.BETAS[1]+2*bg.BETAS[2]*yy+bg.BETAS[3]*yy**2)-aa*bb*(aa+bb*yy)**2*pp
    discriminant = 4*(bg.M4**2*(bg.BETAS[2]**2-bg.BETAS[1]*bg.BETAS[3])
                      +bg.M4*aa*bb*pp*(bb**2*bg.BETAS[1]+aa**2*bg.BETAS[3]-2*aa*bb*bg.BETAS[2]))
    answer["pressure_quadratic_discriminant_without_division"] = sp.expand(sp.discriminant(polynomial, yy)-discriminant)
    return answer


def positive_factors():
    y, c = sp.symbols("y c", positive=True)
    r, s = bg.ALPHA+bg.BETA*y, bg.ALPHA+bg.BETA*c
    return {"inertia": bg.G+bg.F*y**2, "effective_scale_ratio": r,
            "effective_lapse_ratio": s, "weighted_null": s*r**3,
            "intersection_acceleration_coefficient": r**2*(bg.ALPHA**2/bg.G+bg.BETA**2*c**2/(bg.F*y**2))/(2*s**2)}
