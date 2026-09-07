"""Literal June-2026 coefficients and their constant-clock Laurent limits."""

from functools import cache

import sympy as sp


def exact(value):
    """Admit exact finite real numbers, never binary floats or booleans."""
    if isinstance(value, (bool, float)):
        raise TypeError("An exact finite real number is required")
    out = sp.sympify(value)
    if out.has(sp.Float) or out.is_number is not True or out.is_real is not True or out.is_finite is not True:
        raise ValueError("An exact finite real number is required")
    return out


def literal(x, a, g, *, printed=True):
    """Source HTML Eqs. (4), (19); a=a1(phi), g=g1(phi).

    The printed A4 parentheses differ from the earlier standard-Ia source.
    Both are exposed, not silently identified. This difference is regular at
    X=0 when F2(phi,0)!=0 and does not affect the theorem's leading term.
    """
    f = sp.Rational(1, 2)-g*(1+x)
    ax, fx = -a*(1+x), -g
    a3 = 2*(ax-2*fx)*(ax*x-2*f)/(x*(3*ax*x-4*f))
    disputed = 4*(3*f+16*x*fx*ax**2) if printed else 4*(3*f+16*x*fx)*ax**2
    a4 = (-16*x*ax**3+disputed-x**2*f*a3**2
          -(16*x**2*fx-12*x*f)*a3*ax-16*fx*(3*f+4*x*fx)*ax
          +8*f*(x*fx-f)*a3+48*f*fx**2)/(8*(f-x*ax)**2)
    a5 = ((4*fx-2*ax+x*a3)*(-2*ax**2-3*x*ax*a3+4*fx*ax+4*f*a3)
          /(8*(f-x*ax)**2))
    return {"F2": f, "A1": ax, "A2": -ax, "A3": a3, "A4": a4, "A5": a5}


def residues(a, g):
    """Exact theorem-domain data; r=0 is an exclusion, not a false verdict."""
    a, g = exact(a), exact(g)
    f, r = sp.Rational(1, 2)-g, 2*g-a
    if f == 0:
        raise ValueError("F2(phi0,0)=0 is outside this nondegenerate theorem")
    return {"F2_at_zero": f, "r": r, "A3": r, "A4": -r,
            "A5": -r**2/(2*f), "generic_obstruction": r != 0,
            "healthy_tensor_sign": f > 0}


def background(z):
    """Literal benchmark background, dimensionless z=phi/tau.

    a1 may require a removable limit at H=0. The analytic numerator below
    does not divide by H and is the object used in the open-set argument.
    """
    sigma = 1/(1+sp.exp(-z))
    scale = sigma*(1+z**2)**sp.Rational(1, 6)+(1-sigma)*(1+z**2)**sp.Rational(1, 20)
    h = sp.diff(scale, z)/scale
    n = z**2*sp.tanh(z+sp.Rational(1, 10))+sp.tanh(z)
    d = (1+z**2)*h
    g = sp.Rational(2, 3)/sp.cosh(z+sp.Rational(1, 10))**2
    return {"scale": scale, "h": h, "N": n, "D": d, "Z": 9*d-2*n,
            "g1": g, "a1_over_g1": (2*n/d-1)/4, "r_over_g1": (9*d-2*n)/(4*d)}


@cache
def checks():
    x, a, g, z = sp.symbols("X a g z", real=True)
    out = {}
    for printed in (True, False):
        data = literal(x, a, g, printed=printed)
        for key, expected in (("A3", 2*g-a), ("A4", a-2*g),
                              ("A5", -(2*g-a)**2/(1-2*g))):
            out[f"{'printed' if printed else 'Ia'}_{key}_residue"] = sp.factor(
                sp.cancel(x*data[key]).subs(x, 0)-expected)
    printed, ia = literal(x, a, g), literal(x, a, g, printed=False)
    difference = sp.factor(printed["A4"]-ia["A4"])
    out["printed_parenthesis_regular_difference"] = sp.factor(
        difference-3*printed["F2"]*(1-printed["A1"]**2)
        /(2*(printed["F2"]-x*printed["A1"])**2))
    out["difference_has_no_pole"] = sp.cancel(x*difference).subs(x, 0)
    b = background(z)
    out["scale_at_zero"] = b["scale"].subs(z, 0)-1
    out["h_at_zero"] = b["h"].subs(z, 0)
    out["h_slope"] = sp.diff(b["h"], z).subs(z, 0)-sp.Rational(13, 60)
    out["Z_at_zero"] = b["Z"].subs(z, 0)
    out["Z_nonzero_derivative"] = sp.diff(b["Z"], z).subs(z, 0)+sp.Rational(1, 20)
    # l'Hopital uses the nonzero D' just independently differentiated.
    out["removable_a1_over_g1"] = (2/sp.diff(b["D"], z).subs(z, 0)-1)/4-sp.Rational(107, 52)
    out["removable_r_over_g1"] = (-sp.Rational(1, 20)
                                       /(4*sp.diff(b["D"], z).subs(z, 0))+sp.Rational(3, 52))
    return {key: sp.simplify(value) for key, value in out.items()}


def stationary_conditions(f0, w2, chi, f0_prime, w2_prime):
    """Formal zero-derivative necessities for a hypothetical regular extension.

    These retain the displayed coefficient values. They are NOT actual
    Euler equations of the undefined literal X=0 action, and are not solved.
    """
    return {"zero_stress": f0-w2*chi**2, "phi": f0_prime-w2_prime*chi**2,
            "chi": w2*chi}
