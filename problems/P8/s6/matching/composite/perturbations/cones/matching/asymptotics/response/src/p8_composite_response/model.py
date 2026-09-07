"""Exact regular background and physical-metric flux equations.

u=m*T, epsilon=1/Y, e=epsilon*z.  The metric variables here are
g=h_g/sigma and f=h_f/(epsilon*sigma), not canonical normal modes.
All dimensionless action coefficients remain regular at epsilon=0;
the endpoint itself is not a regular two-metric geometry.
"""

from functools import cache

import sympy as sp
from p8_composite_asymptotics import scaled


def expressions(epsilon, u, z, root, q):
    """Algebraic engine shared with outward-rounded interval automatic differentiation.

The square root is supplied through ``value.sqrt()`` for the interval jet
engine; SymPy expressions use ``sp.sqrt``. q is kbar squared.
"""
    e = epsilon*z
    cd = 1+11*epsilon*u*u/4
    ae = cd*cd
    j = 11*u/(cd*z)
    argument = 1+e*(root*root-1)
    x = argument.sqrt() if hasattr(argument, "sqrt") else sp.sqrt(argument)
    lam = (1+e)*(-x*(root+j)+j*e*e*root)/(x+e*root)
    numerator = x-j*e*(1+e)
    denominator = root+j*(1+e)
    kap = numerator/denominator
    zflow = -z*lam
    rflow = ((root*root-1)*lam*(2-e)/(1+e)-j*(2-e+3*e*root*root))/(2*root)
    ag = ae**3*z*z*(e+kap)/(1+e)**3
    af = ae**3*(e+kap)/(kap*(1+e)**3)
    vg = ae*z*z*q/((1+e)*(e+kap))
    vf = ae*kap*q/((1+e)*(e+kap))
    w = ae**3*z*z*(1-e)*(1-kap)/((1+e)**4*(e+kap))
    lock_c = af+epsilon*epsilon*ag
    lock_v = vf+epsilon*epsilon*vg
    output_g, output_f = z/(1+e), 1/(1+e)
    zero = epsilon*0
    matrix = [
        [zero, 1/ag, zero, zero, zero, zero],
        [-vg-w, zero, epsilon*w, zero, zero, zero],
        [zero, zero, zero, 1/af, zero, zero],
        [epsilon*w, zero, -vf-epsilon*epsilon*w, zero, zero, zero],
        [zero, zero, zero, zero, zero, 1/lock_c],
        [zero, zero, zero, zero, -lock_v, zero],
    ]
    forcing = [zero, output_g, zero, output_f, zero, 1+zero]
    return {
        "e": e, "CD_denominator": cd, "A_e": ae, "j": j, "X_bar": x,
        "lambda": lam, "kappa": kap, "lapse_numerator": numerator,
        "lapse_denominator": denominator, "root_denominator": x+e*root,
        "zprime": zflow, "Rprime": rflow,
        "N_g": e/(e+kap), "N_f": kap/(e+kap),
        "eta": 3*(root*root-1)/(1+e)**3,
        "null_over_e": (2-e+3*e*root*root)/(1+e)**3,
        "A_g": ag, "A_f": af, "V_g": vg, "V_f": vf, "W": w,
        "C_lock": lock_c, "V_lock": lock_v,
        "output_g": output_g, "output_f": output_f,
        "matrix": matrix, "forcing": forcing,
    }


@cache
def derive():
    epsilon = sp.Symbol("epsilon", nonnegative=True)
    u = sp.Symbol("u", real=True)
    z, root, q = sp.symbols("z R q", positive=True)
    return {"epsilon": epsilon, "u": u, "z": z, "R": root, "q": q,
            **expressions(epsilon, u, z, root, q)}


def derivative(value, data=None):
    d = derive() if data is None else data
    return (sp.diff(value, d["u"])+sp.diff(value, d["z"])*d["zprime"]
            +sp.diff(value, d["R"])*d["Rprime"])


def initial_root(epsilon):
    """The pinned eta(0)=9 family, not a finite-epsilon R(0)=2 family."""
    argument = 1+3*(1+epsilon)**3
    return argument.sqrt() if hasattr(argument, "sqrt") else sp.sqrt(argument)


@cache
def checks():
    d = derive()
    e, eps, root, z, q = (d[key] for key in ("e", "epsilon", "R", "z", "q"))
    old = scaled.core()
    change = {old["e"]: e, old["eta"]: d["eta"], old["j"]: d["j"]}
    # Root R>0 is prescribed, so sqrt(R²)=R. The old X is positive too.
    old_lam = old["lambda"].subs(change, simultaneous=True)
    old_eta = old["etaprime"].subs(change, simultaneous=True)
    eta_derivative = derivative(d["eta"], d)
    g, f, gp, fp, source = sp.symbols("g f gprime fprime p", real=True)
    action = (d["A_g"]*gp**2+d["A_f"]*fp**2-d["V_g"]*g**2-d["V_f"]*f**2
              -d["W"]*(g-eps*f)**2+2*source*(z*g+f)/(1+e))
    lock, lockp = sp.symbols("L Lprime", real=True)
    locked = action.subs({g: eps*lock, f: lock, gp: eps*lockp, fp: lockp}, simultaneous=True)
    v = sp.Symbol("v", positive=True)
    at0 = {eps: 0, root: v-d["j"].subs(eps, 0)}
    out = {
        "old_full_root_lambda": sp.simplify(old_lam-d["lambda"]),
        "old_full_density_flow": sp.simplify(eta_derivative-old_eta),
        "positive_canonical_null_source": sp.cancel(d["eta"]*e+2/(1+e)**2-d["null_over_e"]),
        "physical_clock_lapse_sum": sp.cancel(d["N_g"]+d["N_f"]-1),
        "literal_locked_action": sp.cancel(locked-(d["C_lock"]*lockp**2-d["V_lock"]*lock**2+2*source*lock)),
        "full_g_force": sp.cancel(sp.diff(action, g)/2-(-(d["V_g"]+d["W"])*g+eps*d["W"]*f+source*d["output_g"])),
        "full_f_force": sp.cancel(sp.diff(action, f)/2-(eps*d["W"]*g-(d["V_f"]+eps**2*d["W"])*f+source*d["output_f"])),
        "limit_Ag": sp.cancel(d["A_g"].subs(at0, simultaneous=True)-z*z/v),
        "limit_Af": sp.cancel(d["A_f"].subs(at0, simultaneous=True)-1),
        "limit_Vg": sp.cancel(d["V_g"].subs(at0, simultaneous=True)-z*z*v*q),
        "limit_Vf": sp.cancel(d["V_f"].subs(at0, simultaneous=True)-q),
        "limit_W": sp.cancel(d["W"].subs(at0, simultaneous=True)-z*z*(v-1)),
        "limit_lock_C": sp.cancel(d["C_lock"].subs(at0, simultaneous=True)-1),
        "limit_lock_V": sp.cancel(d["V_lock"].subs(at0, simultaneous=True)-q),
    }
    return out
