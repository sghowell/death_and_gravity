"""Finite polynomial classification and genuine background countercontrols."""

from functools import cache

import sympy as sp

from . import potential
from .exact import coefficients


def classify(beta):
    """Classify a constant rational J polynomial, not a solution's branch history."""
    beta = coefficients(beta)
    b1, b2, b3 = beta[1:4]
    if b1 == b2 == b3 == 0:
        return {"kind": "endpoint_only", "positive_roots": (), "identically_zero": True}
    roots = []
    if b3 == 0:
        if b2 != 0:
            roots = [(-b1/(2*b2), 1)]
    else:
        discriminant = b2**2-b1*b3
        if discriminant == 0:
            roots = [(-b2/b3, 2)]
        elif discriminant > 0:
            roots = [((-b2+sign*sp.sqrt(discriminant))/b3, 1) for sign in (-1, 1)]
    positive_roots = tuple((sp.simplify(root), multiplicity)
                           for root, multiplicity in roots if root.is_positive is True)
    return {"kind": "genuine", "positive_roots": positive_roots, "identically_zero": False}


@cache
def algebraic_control():
    """Actual expanding central stiff scalar and algebraic de Sitter leaf.

    This is a background solution only. It is not a scalar/vector stability
    claim for the algebraic branch, and is not a counterexample to no-bounce.
    """
    t = sp.Symbol("T", positive=True)
    beta = (0, 1, sp.Rational(-1, 2), 0, sp.Rational(1, 2))
    a = t**sp.Rational(1, 3)
    nu, ni = sp.S.One, 1/(3*t)
    hu, hi = sp.diff(a, t)/a, sp.diff(a, t)/(ni*a)
    phi = sp.sqrt(sp.Rational(2, 3))*sp.log(t)
    rho = sp.diff(phi, t)**2/2
    nh = 2*rho
    d = potential.derive()
    point = dict(zip(d["beta"], beta, strict=True)) | {d["R"]: 1, d["N"]: nu/ni}
    fields = {key: sp.factor(d[key].subs(point)) for key in
              ("A", "B", "J", "rho_i", "pressure_i", "rho_u", "pressure_u", "null_i", "null_u")}
    residuals = {
        "leaf_00": sp.factor(3*hi**2-fields["rho_i"]),
        "leaf_space": sp.factor(-2*sp.diff(hi, t)/ni-3*hi**2-fields["pressure_i"]),
        "central_00": sp.factor(3*hu**2-rho-fields["rho_u"]),
        "central_space": sp.factor(-2*sp.diff(hu, t)-3*hu**2-rho-fields["pressure_u"]),
        "canonical_scalar": sp.simplify(sp.diff(phi, t, 2)+3*hu*sp.diff(phi, t)),
        "unfactored_bianchi": sp.factor(fields["J"]*(hu-hi)),
        "actual_K_dynamic": sp.factor(sp.diff(hu, t)+nh/2),
    }
    false_k = sp.Integer(2)
    false_defect = sp.simplify(sp.diff(hu/sp.sqrt(false_k), t)+nh/(2*false_k**sp.Rational(3, 2)))
    return {"T": t, "beta": beta, "G_u": 1, "G_i": 1,
            "a_u": a, "a_i": a, "n_u": nu, "n_i": ni, "R": 1,
            "N": nu/ni, "c": ni, "H_u": hu, "H_i": hi,
            "phi": phi, "rho": rho, "nh": nh, "fields": fields,
            "residuals": residuals, "false_K_all": false_k,
            "false_K_all_defect": false_defect}


@cache
def endpoint_control():
    """Disconnected G_u=0 exception: u has no geometric equation at all."""
    t = sp.Symbol("T", real=True)
    a = 1+t**2
    h = sp.diff(a, t)/a
    d = potential.derive()
    empty = dict.fromkeys(d["beta"], sp.S.Zero)
    # These are actual variations of the vanished potential, not prescribed
    # numerical residuals. The remaining leaf has a_i=n_i=1 and zero Riemann.
    potential_u = sp.diff(d["L"], d["n_u"]).subs(empty)
    potential_leaf = sp.diff(d["L"], d["n_i"]).subs(empty)
    leaf_h = sp.diff(sp.S.One, t)
    source_phi = sp.S.Zero
    central_euler = sp.factor(potential_u)
    leaf_euler = sp.factor(3*leaf_h**2+potential_leaf)
    scalar_euler = sp.diff(a**3*sp.diff(source_phi, t), t)
    return {"T": t, "G_u": 0, "G_i": 1, "beta": (0, 0, 0, 0, 0),
            "a_u": a, "a_i": 1, "n_u": 1, "n_i": 1,
            "phi": 0, "rho": 0, "nh": 0, "H_u": h,
            "H_u_prime_at_zero": sp.diff(h, t).subs(t, 0),
            "central_Euler": central_euler, "leaf_Euler": leaf_euler,
            "scalar_Euler": scalar_euler}


def original_auxiliary_mapping():
    """S6.13/14 has leaf beta3=p_i, beta0=b_i and central endpoint B."""
    p, endpoint, central_b = sp.symbols("p_i b_i B_central", real=True)
    d = potential.derive()
    point = dict(zip(d["beta"], (endpoint, 0, 0, p, 0), strict=True))
    return {"p_i": p, "b_i": endpoint, "B_central": central_b,
            "J": sp.factor(d["J"].subs(point)),
            "rho_i": sp.factor(d["rho_i"].subs(point)),
            "rho_u_link": sp.factor(d["rho_u"].subs(point)),
            "central_endpoint_density": 2*central_b,
            "central_endpoint_pressure": -2*central_b}


def checks():
    d = algebraic_control()
    original = original_auxiliary_mapping()
    r = potential.derive()["R"]
    return dict(d["residuals"]) | {
        "original_auxiliary_J": sp.expand(original["J"]-original["p_i"]*r**2),
        "original_auxiliary_u_density": sp.cancel(original["rho_u_link"]-6*original["p_i"]/r),
        "algebraic_false_formula_nonzero_control": sp.simplify(
            d["false_K_all_defect"]+1/(6*sp.sqrt(2)*d["T"]**2)),
    }
