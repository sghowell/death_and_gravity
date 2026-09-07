"""Exact past-square identity and endpoint-matched cubic samplers."""

import sympy as sp
from p8a_maxwell.domain import nonnegative


def cubic(value):
    value = sp.sympify(value)
    return 3*value**2-2*value**3


def moments():
    x = sp.Symbol("x", real=True)
    p = cubic(x)
    return {"zeroth": sp.integrate(p*p, (x, 0, 1)),
            "first": sp.integrate(sp.diff(p, x)**2, (x, 0, 1)),
            "second": sp.integrate(sp.diff(p, x, 2)**2, (x, 0, 1)),
            "pole_first": sp.integrate(sp.cancel((sp.diff(p, x)/x)**2), (x, 0, 1)),
            "pole_zero": sp.integrate(sp.cancel((p/x**2)**2), (x, 0, 1)),
            "cross": sp.integrate(p*sp.diff(p, x), (x, 0, 1))}


def lower(h, ratio):
    """tau times 3 integral_past[(u'-Hu)^2-u'^2] under H<=-h/tau."""
    h = nonnegative(h, positive=True)
    ratio = nonnegative(ratio, positive=True)
    return 3*h+sp.Rational(39, 35)*h*h*ratio


def identities():
    t = sp.Symbol("t", real=True)
    u, g, h = (sp.Function(name)(t) for name in ("u", "g", "H"))
    rho = 3*(sp.diff(h, t)+h*h)
    data = moments()
    return {"past_Ricci_square_with_endpoint_K": sp.expand(
                rho*u*u-3*(sp.diff(u, t)-h*u)**2+3*sp.diff(u, t)**2
                -sp.diff(3*h*u*u, t)),
            "future_index_square_with_initial_K": sp.expand(
                3*sp.diff(g, t)**2+rho*g*g-3*(sp.diff(g, t)-h*g)**2
                -sp.diff(3*h*g*g, t)),
            "cubic_zeroth_norm": data["zeroth"]-sp.Rational(13, 35),
            "cubic_first_norm": data["first"]-sp.Rational(6, 5),
            "cubic_second_norm": data["second"]-12,
            "relative_pole_first_norm": data["pole_first"]-12,
            "relative_pole_zero_norm": data["pole_zero"]-sp.Rational(13, 3),
            "history_cross_endpoint": data["cross"]-sp.Rational(1, 2)}
