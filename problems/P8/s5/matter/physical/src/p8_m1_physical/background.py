"""Rational compact time and the pinned nonlinear invariant coefficients."""

from functools import cache

import sympy as sp
from p8_m1 import nonlinear, series

rtime = sp.Symbol("rtime", real=True)


def point(value):
    if value is None:
        return rtime
    value = sp.Rational(value)
    if abs(value) > 1:
        raise ValueError("Rational compact time must lie in [-1,1]")
    return value


@cache
def functions(time_point=None):
    r = point(time_point)
    x = 2*r/(1+r**2)
    y = (1-r**2)/(1+r**2)
    J = series.compact_even(nonlinear.d*nonlinear.functions()["J"])
    J = sp.factor(J.subs(series.x, x))
    H, l = 4*x, y**11/10
    theta, lam = x*(4-y**6), 1-sp.Rational(3, 2)*y**6
    return {"x": x, "y": y, "H": H, "Hdot": 4*(1-2*x**2), "Adot": 4+40*x**2,
            "l": l, "theta": sp.factor(theta), "lam": sp.factor(lam),
            "w": sp.factor(l*(sp.Rational(3, 2)*y**6-1)), "J": J}


@cache
def lapse_jets(time_point=None):
    bg = functions(time_point)
    mapping = {series.x: bg["x"], series.y: bg["y"]}
    return {key: tuple(sp.factor(expr.subs(mapping, simultaneous=True)) for expr in row)
            for key, row in series.compact_jets().items()}


def compact_checks():
    r = rtime
    bg = functions()
    u = 2*r/(1-r**2)
    return {
        "unit_semicircle": sp.cancel(bg["x"]**2+bg["y"]**2-1),
        "inverse_time": sp.cancel(bg["x"]-u*bg["y"]),
        "lapse_J_bridge": sp.cancel(lapse_jets()["A"][2]+2*bg["J"]),
        "Theta_bridge": sp.cancel(lapse_jets()["B"][1]-2*bg["theta"]),
        "Lambda_bridge": sp.cancel(lapse_jets()["C"][1]+bg["lam"]/2),
        "matter_bridge": sp.cancel(lapse_jets()["L"][1]+bg["w"]),
    }
