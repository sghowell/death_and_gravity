"""Conformal-vacuum stress in the explicit FK dictionary.

All jet components below have the common factor hbar/(2880*pi**2)
removed. beta_m multiplies the conserved FK variation of integral R**2.
It is independent of every scalar prescription in earlier checkpoints.
"""

import sympy as sp

H0, H1, H2, H3 = sp.symbols("H Hdot Hddot Hthird", real=True)
BETA_M = sp.Symbol("beta_M", real=True)


def reference_jets(h=H0, hdot=H1, hddot=H2, hthird=H3, beta_m=BETA_M):
    h, hd, hdd, h3, beta = map(sp.sympify, (h, hdot, hddot, hthird, beta_m))
    i_rho = 18*hd**2-108*h**2*hd-36*h*hdd
    i_trace = -36*(h3+4*hd**2+7*h*hdd+12*h**2*hd)
    rho = 186*h**4+beta*i_rho
    trace = 744*h**2*(h**2+hd)+beta*i_trace
    return {"rho": sp.expand(rho), "pressure": sp.expand((rho-trace)/3),
            "trace": sp.expand(trace), "EED": sp.expand(rho-trace/2),
            "I_rho": sp.expand(i_rho), "I_trace": sp.expand(i_trace),
            "I_EED": sp.expand(i_rho-i_trace/2)}


def reference_stress(h, hdot, hddot, hthird, beta_m, *, hbar):
    """Formal physical components; hypotheses and the common scheme are external."""
    return {name: sp.sympify(hbar)*value/(2880*sp.pi**2)
            for name, value in reference_jets(h, hdot, hddot, hthird, beta_m).items()}


def jet_coefficients():
    result = {}
    for name, expression in reference_jets().items():
        pieces = {"universal": expression.subs(BETA_M, 0),
                  "beta_M": sp.diff(expression, BETA_M)}
        result[name] = {piece: {",".join(map(str, powers)): str(value)
                               for powers, value in sp.Poly(poly, H0, H1, H2, H3).terms()
                               if value != 0}
                        for piece, poly in pieces.items()}
    return result


def identities():
    t = sp.Symbol("t", real=True)
    h = sp.Function("H")(t)
    hd, hdd, h3 = (sp.diff(h, t, n) for n in (1, 2, 3))
    r, r00, rii = 6*(hd+2*h**2), 3*(hd+h**2), -(hd+3*h**2)
    ricci_squared = r00**2+3*rii**2
    h3_rho = r00**2-sp.Rational(2, 3)*r*r00-ricci_squared/2+r**2/4
    h3_trace = -ricci_squared+r**2/3
    box_r = sp.diff(r, t, 2)+3*h*sp.diff(r, t)
    i_rho = 2*r*r00-r**2/2-6*h*sp.diff(r, t)
    data = reference_jets(h, hd, hdd, h3)
    euler = 24*h**2*(h**2+hd)
    central_a = sp.Rational(31, 180)
    # HH has g00=-1, Ric00=-FK Ric00 and the same scalar R. Its H3_00
    # therefore equals -H3_FK, while physical covariant T00 is unchanged.
    hh_h3_rho = -h3_rho
    hh_rho_pi2 = -2*central_a*hh_h3_rho/16
    return {
        "FK_scalar_and_G00": sp.simplify(r00-r/2+3*h**2),
        "quadratic_tensor_density": sp.simplify(h3_rho-3*h**4),
        "quadratic_tensor_trace": sp.simplify(h3_trace-euler/2),
        "HH_to_FK_physical_density_sign": sp.simplify(hh_rho_pi2-31*h**4/480),
        "universal_Euler_anomaly": sp.simplify(data["trace"].subs(BETA_M, 0)-31*euler),
        "FK_R_squared_density": sp.simplify(data["I_rho"]-i_rho),
        "FK_R_squared_trace": sp.simplify(data["I_trace"]+6*box_r),
        "full_reference_conservation": sp.simplify(sp.diff(data["rho"], t)
                                                     +3*h*(data["rho"]+data["pressure"])),
        "EED_not_energy_when_trace_nonzero": sp.expand(data["EED"]-data["rho"]+data["trace"]/2),
        "radiation_finite_ambiguity_vanishes": sp.simplify(
            data["I_EED"].subs(h, 1/(2*t)).doit()),
    }
