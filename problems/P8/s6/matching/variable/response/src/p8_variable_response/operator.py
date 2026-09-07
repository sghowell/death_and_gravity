"""Full coupled operator and actual physical Cauchy maps; no locked ansatz.

The canonical variables are divided by M and use u=T/tau. Thus the physical
Cauchy vector below is (gamma_g, tau*gamma_g,T, gamma_f, tau*gamma_f,T).
The derivative vector is not silently identified with canonical momenta.
"""

from functools import cache

import sympy as sp
from p8_variable_beta import canonical

from .exact import MU, number, parameters, side


@cache
def derive():
    d = canonical.derive()
    u, c, k = (d[key] for key in ("u", "c", "kbar"))
    den = 1+u**2
    a = den**2
    om, theta, cross = d["omega"], d["theta_sum"], d["q_squared"]*d["D"]
    mass = sp.factor(d["mass_squared"])
    fields = {
        "A": d["V_LL"],
        "C": cross-2*sp.diff(om, u)-2*om*theta,
        "d_cross": -2*om/u,
        "E": cross-2*om*theta,
        "f_cross": 2*om/u,
        "b_analytic": d["V_HH"]-mass,
    }
    fields = {key: sp.factor(value) for key, value in fields.items()}
    fs = sp.sqrt(c*den**12+8)/(2*sp.sqrt(c)*den**3)
    fr = sp.sqrt(2)*den**3/sp.sqrt(c*den**12+8)
    return {**d, **fields, "d": den, "a": a, "b": 2/den**2,
            "u": u, "c": c, "kbar": k, "mass": mass,
            "B": mass+fields["b_analytic"], "f_sum": fs, "f_relative": fr}


def first_order():
    """u generator of (l,l_u,Q,Q_u), including BOTH moving-weight terms."""
    d = derive()
    return sp.Matrix([[0, 1, 0, 0],
                      [-d["A"], 0, -d["C"], -d["u"]*d["d_cross"]],
                      [0, 0, 0, 1],
                      [-d["E"], -d["u"]*d["f_cross"], -d["B"], 0]])


@cache
def physical_map():
    """Exact delta-dependent map, including derivatives of every coefficient."""
    d = derive()
    u = d["u"]
    ag, bg = 1/d["f_sum"], -d["w2"]/d["f_relative"]
    af, bf = ag, d["w1"]/d["f_relative"]
    return sp.Matrix([[ag, 0, bg, 0], [sp.diff(ag, u), ag, sp.diff(bg, u), bg],
                      [af, 0, bf, 0], [sp.diff(af, u), af, sp.diff(bf, u), bf]])


def physical_at(u_value, delta, kbar=1, *, tau=1):
    """Map to raw proper-time derivatives for a specified positive exact tau."""
    delta, _ = parameters(delta, kbar)
    u_value, tau = number(u_value, "u"), number(tau, "tau")
    if (sp.Rational(1, 10)-abs(u_value)).is_nonnegative is not True or tau.is_positive is not True:
        raise ValueError("Require |u|<=1/10 and tau>0")
    d = derive()
    return sp.diag(1, 1/tau, 1, 1/tau)*physical_map().subs({d["u"]: u_value, d["c"]: 2+delta})


def outer_map(r, sign):
    """Y=(l,l_u,Q/sqrt(r),(u Q_u-Q/2)/(mu sqrt(r))) to derivative data."""
    sign = side(sign)
    return sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0],
                      [0, 0, sp.sqrt(r), 0],
                      [0, 0, sign/(2*sp.sqrt(r)), sign*MU/sp.sqrt(r)]])


def rotation(log_radius):
    angle = MU*log_radius
    return sp.diag(sp.eye(2), sp.Matrix([[sp.cos(angle), sp.sin(angle)],
                                       [-sp.sin(angle), sp.cos(angle)]]))


def outer_matrix(r, sign, delta=0):
    """s=log(r) generator. delta=0 is only the punctured limiting operator."""
    sign = side(sign)
    d = derive()
    values = {d["u"]: sign*r, d["c"]: 2+delta}
    a, cc, dd, ee, ff, bb = (d[key].subs(values) for key in
                            ("A", "C", "d_cross", "E", "f_cross", "B"))
    return sp.Matrix([[0, sign*r, 0, 0],
                      [-sign*r*a, 0, -sign*r**sp.Rational(3, 2)*(cc+dd/2),
                       -sign*r**sp.Rational(3, 2)*dd*MU],
                      [0, 0, 0, MU],
                      [-r**sp.Rational(3, 2)*ee/MU,
                       -sign*r**sp.Rational(5, 2)*ff/MU,
                       -(r**2*bb-sp.Rational(1, 4))/MU, 0]])


def inner_generator(u, radius, a_coeff, c_coeff, d_coeff, e_coeff, f_coeff, b_bounded):
    """Exact z generator, r_e=radius, Q=sqrt(r_e)*psi, z'=1/r_e.

    b_bounded means B-10/r_e^2, NOT the analytic b_analytic=B-mass.
    The relation radius^2=u^2+delta/8 is part of this chart.
    """
    v = sp.Rational(3, 4)*(1-u**2/radius**2)
    return sp.Matrix([
        [0, radius, 0, 0],
        [-radius*a_coeff, 0,
         -radius**sp.Rational(3, 2)*c_coeff-u**2*d_coeff/(2*sp.sqrt(radius)),
         -u*d_coeff*sp.sqrt(radius)*MU],
        [0, 0, 0, MU],
        [-radius**sp.Rational(3, 2)*e_coeff/MU,
         -u*radius**sp.Rational(3, 2)*f_coeff/MU,
         -(MU**2+v+radius**2*b_bounded)/MU, 0],
    ])


def boundary_map(eta, sign):
    """Outer Y to inner (l,l_u,psi,psi_z/mu) at u=sign*r_match."""
    sign = side(sign)
    return sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0],
                      [0, 0, (1+eta)**(-sp.Rational(1, 4)), 0],
                      [0, 0, sign*eta/(2*MU*(1+eta)**sp.Rational(3, 4)),
                       sign*(1+eta)**sp.Rational(1, 4)]])


@cache
def checks():
    d = derive()
    u, c = d["u"], d["c"]
    r = sp.Symbol("r", positive=True)
    output = {
        "sum_normalization": sp.factor(d["f_sum"]**2-d["f_sum_squared"]),
        "relative_normalization": sp.factor(d["f_relative"]**2-d["f_relative_squared"]),
        "analytic_correction_not_mass_remainder": sp.factor(d["B"]-d["V_HH"]),
        "indicial_mu": MU**2+sp.Rational(1, 4)-10,
        "physical_map_determinant": sp.simplify(physical_map().det()-1/(d["f_sum"]**2*d["f_relative"]**2)),
    }
    for sign in (-1, 1):
        w = outer_map(r, sign)
        direct = w.inv()*(sign*r*first_order().subs(u, sign*r)*w-r*sp.diff(w, r))
        diff = direct-outer_matrix(r, sign, c-2)
        output[f"full_outer_chain_rule_side_{sign}"] = sum(sp.factor(v)**2 for v in diff)
    # Full inner chain rule with independent, generic operator coefficients.
    delta = sp.Symbol("delta", positive=True)
    re = sp.sqrt(u**2+delta/8)
    aa, cc, dd, ee, ff, bb = sp.symbols("A C d E f b", real=True)
    wi = sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, sp.sqrt(re), 0],
                    [0, 0, u/(2*re**sp.Rational(3, 2)), MU/sp.sqrt(re)]])
    generic = sp.Matrix([[0, 1, 0, 0], [-aa, 0, -cc, -u*dd],
                         [0, 0, 0, 1], [-ee, -u*ff, -10/re**2-bb, 0]])
    transformed = wi.inv()*(re*generic*wi-re*sp.diff(wi, u))
    output["full_inner_chain_rule"] = sum(sp.simplify(v)**2 for v in
                                          transformed-inner_generator(u, re, aa, cc, dd, ee, ff, bb))
    for sign in (-1, 1):
        eta = sp.Symbol("eta", positive=True)
        direct = (wi.inv()*outer_map(r, sign)).subs({u: sign*r, delta: 8*r**2*eta})
        # r,eta,1+eta are positive; factor before reducing fractional powers.
        output[f"matching_boundary_side_{sign}"] = sum(
            sp.simplify(sp.powdenest(sp.factor(v), force=True))**2
            for v in direct-boundary_map(eta, sign))
    return output
