"""Analytic scaled full CD ODE; epsilon=0 is not a regular two-metric geometry."""

from functools import cache

import sympy as sp
from p8_composite import reconstruction


@cache
def core():
    e = sp.Symbol("e", nonnegative=True)
    eta, j, ac = sp.symbols("eta j a_current", real=True)
    xb = sp.sqrt(1+eta*e*(1+e)**3/3)
    rb = sp.sqrt(1+eta*(1+e)**3/3)
    lam = (1+e)*(-xb*(rb+j)+j*e**2*rb)/(xb+e*rb)
    kappa = (xb-j*e*(1+e))/(rb+j*(1+e))
    return {"e": e, "eta": eta, "j": j, "a_current": ac,
            "X_bar": xb, "R_bar": rb, "lambda": lam, "kappa": kappa,
            "eprime": -e*lam,
            "etaprime": -3*j*(eta*e+2/(1+e)**2)+2*eta*lam,
            "jprime": ac+j*lam,
            "null_over_e": eta*e+2/(1+e)**2}


@cache
def cd_family():
    """Exact epsilon>0 CD reconstruction in (z=e/epsilon,eta), physical u=mT."""
    old = core()
    epsilon = sp.Symbol("epsilon", nonnegative=True)
    u = sp.Symbol("u", real=True)
    z = sp.Symbol("z", positive=True)
    eta = sp.Symbol("eta", real=True)
    a0 = sp.Symbol("a0", positive=True)
    denominator = 1+a0*epsilon*u**2/4
    e = epsilon*z
    j = a0*u/(denominator*z)
    ac = a0*(1-a0*epsilon*u**2/4)/(denominator**2*z)
    at = {old["e"]: e, old["eta"]: eta, old["j"]: j, old["a_current"]: ac}
    lam, kap = (old[key].subs(at, simultaneous=True) for key in ("lambda", "kappa"))
    return {"epsilon": epsilon, "u": u, "z": z, "eta": eta, "a0": a0,
            "e": e, "j": j, "a_current": ac, "A_e": denominator**2,
            "h": a0*epsilon*u/denominator,
            "lambda": lam, "kappa": kap, "zprime": -z*lam,
            "etaprime": old["etaprime"].subs(at, simultaneous=True),
            "N_g": e/(e+kap), "N_f": kap/(e+kap),
            "null_over_e": old["null_over_e"].subs(at, simultaneous=True)}


def family_derivative(value, data=None):
    data = cd_family() if data is None else data
    return (sp.diff(value, data["u"])+sp.diff(value, data["z"])*data["zprime"]
            +sp.diff(value, data["eta"])*data["etaprime"])


@cache
def canonical_jets():
    """Regular coefficient templates. Rates are in u=mT, with M=m=1 units.

    f_relative²=e²*A_e³*shape_relative/4; its vanishing factor e² is
    removed before logarithmic differentiation. No inverse of e occurs.
    """
    e = sp.Symbol("e", nonnegative=True)
    k = sp.Symbol("kappa", positive=True)
    lam, ld, h, hd, kd, kdd = sp.symbols("lambda lambda_prime h hprime kappa_prime kappa_second", real=True)
    ed, edd = -e*lam, e*(lam**2-ld)
    shape_s = (e+k)*(1+k*e**2)/(k*(1+e)**3)
    shape_r = (e+k)/((1+e)**3*(1+k*e**2))
    coordinates, velocities, accelerations = (e, k), (ed, kd), (edd, kdd)

    def logarithmic_jets(shape):
        logarithm = sp.log(shape)
        first = sum(sp.diff(logarithm, x)*dx for x, dx in zip(coordinates, velocities, strict=True))
        second = sum(sp.diff(logarithm, x)*ddx for x, ddx in zip(coordinates, accelerations, strict=True))
        second += sum(sp.diff(logarithm, x, xx)*dx*dxx
                      for x, dx in zip(coordinates, velocities, strict=True)
                      for xx, dxx in zip(coordinates, velocities, strict=True))
        return first, second

    ls, lss = logarithmic_jets(shape_s)
    lr, lrr = logarithmic_jets(shape_r)
    ts, tr = (3*h+ls)/2, (3*h+lr)/2-lam
    ns, nr = (3*hd+lss)/2+ts**2, (3*hd+lrr)/2-ld+tr**2
    omega = e*sp.sqrt(k)/(1+k*e**2)*(lam-kd/(2*k))
    omega_prime = sum(sp.diff(omega, x)*dx for x, dx in
                      zip((e, k, lam, kd), (ed, kd, ld, kdd), strict=True))
    mass = (1-e)*(1-k)*(1+k*e**2)/((1+e)*(e+k)**2)
    mass_prime = sp.diff(mass, e)*ed+sp.diff(mass, k)*kd
    mass_second = (sp.diff(mass, e)*edd+sp.diff(mass, k)*kdd
                   +sp.diff(mass, e, 2)*ed**2+2*sp.diff(mass, e, k)*ed*kd
                   +sp.diff(mass, k, 2)*kd**2)
    source = (e*k-1)/((1+e)*sp.sqrt(k))
    return {"e": e, "kappa": k, "lambda": lam, "lambda_prime": ld,
            "h": h, "hprime": hd, "kappa_prime": kd, "kappa_second": kdd,
            "shape_sum": shape_s, "shape_relative": shape_r,
            "theta_sum": ts, "theta_relative": tr, "N_sum": ns, "N_relative": nr,
            "omega": omega, "omega_prime": omega_prime, "mass_squared": mass,
            "mass_squared_prime": mass_prime, "mass_squared_second": mass_second,
            "physical_source_relative": source}


def checks():
    d = core()
    old = reconstruction.derive()
    e, eta, j = (d[key] for key in ("e", "eta", "j"))
    change = {old["y"]: 1/e, old["rho"]: eta*e**2, old["h"]: j*e}
    # The positive-root scaling is used explicitly: X=Xbar/e, Y=-e*Rbar.
    literal_y = old["yprime"].xreplace({old["X"]: d["X_bar"]/e, -old["Y"]: e*d["R_bar"]})
    literal_y = literal_y.subs(change, simultaneous=True)
    # A second direct algebraic route avoids any need to simplify sqrt(e²)
    # at the endpoint: start from the positive e>0 root dictionary.
    direct_lambda = (1/e)*(1+1/e)*(-d["X_bar"]*d["R_bar"]
                      -j*e*(d["X_bar"]/e-e*d["R_bar"]))/(d["X_bar"]/e+d["R_bar"])*e
    direct_rho = -3*j*e*(eta*e**2+2*e/(1+e)**2)
    direct_eta = direct_rho/e**2+2*eta*d["lambda"]
    family = cd_family()
    u = family["u"]
    ae = family["A_e"]
    k = sp.Symbol("kappa", positive=True)
    w1 = k*e**2/(1+k*e**2)
    source = (w1-e/(1+e))/(e*sp.sqrt(k)/(1+k*e**2))
    return {
        "immutable_reconstruction_ratio_bridge": sp.cancel(e*literal_y-d["lambda"]),
        "positive_root_scaled_y_flow": sp.cancel(direct_lambda-d["lambda"]),
        "full_shared_matter_scaled_density_flow": sp.cancel(direct_eta-d["etaprime"]),
        "full_CD_physical_scale": sp.cancel(sp.diff(ae, u)/ae-family["h"]),
        "full_CD_scaled_hprime": sp.cancel(sp.diff(family["h"], u)/family["e"]-family["a_current"]),
        "full_CD_hsecond_scaled_limit": sp.cancel(sp.diff(family["h"], u, 2)/family["e"]).subs(family["epsilon"], 0),
        "full_CD_lapse_sum": sp.cancel(family["N_g"]+family["N_f"]-1),
        "physical_composite_source_projection": sp.cancel(source-(e*k-1)/((1+e)*sp.sqrt(k))),
        "null_source_has_positive_scaled_limit": d["null_over_e"].subs(e, 0)-2,
    }
