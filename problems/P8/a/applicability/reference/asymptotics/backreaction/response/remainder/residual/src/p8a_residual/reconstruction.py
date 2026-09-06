"""Exact conserved physical stress on the actual conformal metric.

FK signature +---, R=+6*(Hdot+2H^2). S incorporates the raw-scheme Box-R
anomaly, not a new finite Wick shift. J'=U'S and L'=h U^2, with the actual
common-past data fixing both integrals and the radiation integration constant.
"""

import sympy as sp


def potential(a, eta):
    return -sp.diff(a, eta, 2)/a


def effective_s(remainder, retarded_kernel, u):
    remainder, retarded_kernel, u = map(sp.sympify, (remainder, retarded_kernel, u))
    return remainder-retarded_kernel/2+u/10


def histories(a, s, eta, eta_start):
    x = sp.Dummy("s", real=True)
    u = potential(a, eta)
    h = sp.diff(a, eta)/a
    return {"J": sp.Integral((sp.diff(u, eta)*s).subs(eta, x), (x, eta_start, eta)),
            "L": sp.Integral((h*u**2).subs(eta, x), (x, eta_start, eta))}


def curvature_tensor(a, eta):
    """Physical components of FK I_R2; coefficients are not set by A.3."""
    proper = lambda value: sp.diff(value, eta)/a
    h = sp.diff(a, eta)/a**2
    r = 6*(proper(h)+2*h**2)
    rtt = 3*(proper(h)+h**2)
    rho = 2*r*rtt-r**2/2-6*h*proper(r)
    pressure = -2*r*(proper(h)+3*h**2)+r**2/2+2*proper(proper(r))+4*h*proper(r)
    return {"density": rho, "pressure": pressure, "EED": (rho+3*pressure)/2}


def trace(a, wick, eta, *, hbar=1, gamma=0):
    """Conserved GS raw-H_lambda trace in FK physical conventions."""
    h = sp.diff(a, eta)/a
    proper_h = h/a
    box = lambda value: (sp.diff(value, eta, 2)+2*h*sp.diff(value, eta))/a**2
    r = 6*sp.diff(a, eta, 2)/a**3
    v1 = r**2/288-proper_h**2*(sp.diff(proper_h, eta)/a+proper_h**2)/60-box(r)/120
    return -box(wick)/2-sp.sympify(hbar)*v1/(4*sp.pi**2)-6*sp.sympify(gamma)*box(r)


def components(a, s, j, ell, eta, *, hbar=1, gamma=0, radiation_constant=0):
    """Complete trace/conservation reconstruction with explicit finite freedom.

    s=Rmode-K/2+U/10, j'=U'*s, ell'=(a'/a)*U^2. The generic homogeneous
    radiation_constant has physical units a^4*rho. The specified prepared
    state fixes it to zero; conservation alone does not do so. The calibrated
    numerical theorem uses gamma=0, not an assumption of scheme independence.
    """
    h = sp.diff(a, eta)/a
    u = potential(a, eta)
    first, second = sp.diff(s, eta), sp.diff(s, eta, 2)
    rho = h**4/960+(-h*first+(h**2-u)*s+j)/8-ell/32
    pressure = (5*h**4+4*h**2*u)/2880+(second-3*h*first+(3*h**2+u)*s+j)/24+(u**2-ell)/96
    eed = (3*h**4+2*h**2*u)/960+(second-4*h*first+4*h**2*s+2*j)/16+(u**2-2*ell)/64
    pref = sp.sympify(hbar)/(sp.pi**2*a**4)
    finite = curvature_tensor(a, eta)
    c = sp.sympify(radiation_constant)/a**4
    return {"density": pref*rho+sp.sympify(gamma)*finite["density"]+c,
            "pressure": pref*pressure+sp.sympify(gamma)*finite["pressure"]+c/3,
            "EED": pref*eed+sp.sympify(gamma)*finite["EED"]+c}


def fix_radiation_constant(scale_at_start, prescribed_density, reconstructed_density):
    """Fix, rather than invent, the homogeneous solution from actual past data."""
    return sp.sympify(scale_at_start)**4*(prescribed_density-reconstructed_density)


def identities():
    eta = sp.Symbol("eta", positive=True)
    a, s, j, ell = (sp.Function(name)(eta) for name in ("a", "S", "J", "L"))
    h = sp.diff(a, eta)/a
    u = potential(a, eta)
    values = components(a, s, j, ell, eta)
    wick = (s-u/10)/(4*sp.pi**2*a**2)
    rules = {sp.diff(j, eta): sp.diff(u, eta)*s, sp.diff(ell, eta): h*u**2}
    finite = curvature_tensor(a, eta)
    r = 6*sp.diff(a, eta, 2)/a**3
    box_r = (sp.diff(r, eta, 2)+2*h*sp.diff(r, eta))/a**2
    return {"exact_raw_scheme_trace": sp.simplify(values["density"]-3*values["pressure"]-trace(a, wick, eta)),
            "exact_full_conservation": sp.simplify((sp.diff(values["density"], eta)
                + 3*h*(values["density"]+values["pressure"])).subs(rules)),
            "exact_EED_definition": sp.simplify(values["EED"]-(values["density"]+3*values["pressure"])/2),
            "finite_I_trace": sp.simplify(finite["density"]-3*finite["pressure"]+6*box_r),
            "finite_I_conservation": sp.simplify(sp.diff(finite["density"], eta)
                + 3*h*(finite["density"]+finite["pressure"])),
            "effective_anomaly_Wick_term": sp.simplify(wick-(s/(4*sp.pi**2*a**2)+r/(240*sp.pi**2)))}


def controls():
    eta = sp.Symbol("eta", positive=True)
    a, s, j, ell = (sp.Function(name)(eta) for name in ("a", "S", "J", "L"))
    h = sp.diff(a, eta)/a
    u = potential(a, eta)
    results = {}
    for name, jj, ll in (("discard_J", sp.Integer(0), ell), ("discard_L", j, sp.Integer(0))):
        value = components(a, s, jj, ll, eta)
        defect = sp.diff(value["density"], eta)+3*h*(value["density"]+value["pressure"])
        results[name] = sp.simplify(defect.subs({sp.diff(j, eta): sp.diff(u, eta)*s,
                                               sp.diff(ell, eta): h*u**2}))
    wick_without_anomaly = s/(4*sp.pi**2*a**2)
    full = components(a, s, j, ell, eta)
    results["omit_U_over_10_changes_trace"] = sp.simplify(
        full["density"]-3*full["pressure"]-trace(a, wick_without_anomaly, eta))
    results["wrong_past_constant_density_shift"] = 1/a**4
    results["wrong_past_constant_can_still_be_conserved"] = True
    return results
