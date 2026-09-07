"""Exact pole removal and the coefficient-space analytic inverse.

Q=D*q is an identity for the literal positive-delta parent. No locked ansatz
or truncation of the full TT equations is used to define the analytic branch.
"""

from functools import cache

import sympy as sp
from p8_variable_response import operator

from .exact import integer


@cache
def derive():
    previous = operator.derive()
    u = previous["u"]
    delta, momentum = sp.symbols("delta K", real=True)
    substitutions = {previous["c"]: 2+delta, previous["kbar"]**2: momentum}
    fields = {name: sp.factor(previous[name].subs(substitutions)) for name in
              ("A", "C", "d_cross", "E", "f_cross", "b_analytic", "omega",
               "f_sum", "f_relative", "w1", "w2")}
    d = 1+u**2
    pole = (2+delta)*d**4-2
    numerator = 16*(1-u**2)*(8/((2+delta)*d**14)+1/d**2)
    return {**fields, "u": u, "delta": delta, "K": momentum, "D": pole,
            "N": numerator, "mass": numerator/pole,
            "light_q": fields["C"]*pole+u*fields["d_cross"]*sp.diff(pole, u),
            "light_qprime": u*fields["d_cross"]*pole,
            "heavy_q": sp.diff(pole, u, 2)+numerator+pole*fields["b_analytic"]}


def diagonal(index):
    index = integer(index)
    return 8*(index**2+3*index+12)


def inverse_monomial(delta_degree, u_degree):
    """Exact polynomial L0 inverse of delta^j*u^n, including all delta shifts."""
    j, n = integer(delta_degree), integer(u_degree)
    d = derive()
    coefficient = sp.Rational(1, diagonal(n))
    output = coefficient*d["delta"]**j*d["u"]**n
    while n >= 2:
        coefficient *= -sp.Rational(n*(n-1), diagonal(n-2))
        j, n = j+1, n-2
        output += coefficient*d["delta"]**j*d["u"]**n
    return output


def leading_operator(expression):
    d = derive()
    u, delta = d["u"], d["delta"]
    return sp.expand((delta+8*u**2)*sp.diff(expression, u, 2)
                     +32*u*sp.diff(expression, u)+96*expression)


def low_jets():
    d = derive()
    u, delta, momentum = d["u"], d["delta"], d["K"]
    x = sp.Symbol("x", real=True)
    return {"q_odd": 3*u/40,
            "q_even": -momentum*(4*u**2/55+7*delta/2640),
            "l_even": 1-(momentum+sp.Rational(18, 5))*u**2/2,
            "H_odd": 3*x**3/5+3*x/40,
            "H_even": -momentum*(32*x**4/55+31*x**2/330+sp.Rational(7, 2640)),
            "x": x}


@cache
def checks():
    d, old, jets = derive(), operator.derive(), low_jets()
    u, delta, momentum = d["u"], d["delta"], d["K"]
    q = sp.Function("q")(u)
    out = {"literal_mass_numerator": sp.factor(
        old["mass"].subs(old["c"], 2+delta)*d["D"]-d["N"])}
    exact_heavy = sp.diff(d["D"]*q, u, 2)+(d["mass"]+d["b_analytic"])*d["D"]*q
    reduced = d["D"]*sp.diff(q, u, 2)+2*sp.diff(d["D"], u)*sp.diff(q, u)+d["heavy_q"]*q
    out["full_pole_removed_heavy_equation"] = sp.factor(exact_heavy-reduced)
    out["full_pole_removed_light_equation"] = sp.expand(
        d["C"]*d["D"]*q+u*d["d_cross"]*sp.diff(d["D"]*q, u)
        -d["light_q"]*q-d["light_qprime"]*sp.diff(q, u))
    origin = {u: 0, delta: 0}
    out["light_center"] = sp.simplify(d["A"].subs(origin)-momentum-sp.Rational(18, 5))
    out["moving_weight_center"] = sp.simplify(d["f_cross"].subs(origin)+sp.Rational(48, 5))
    out["E_constant"] = sp.simplify(d["E"].subs(origin))
    out["E_delta"] = sp.simplify(sp.diff(d["E"], delta).subs(origin)-2*momentum/5)
    out["E_u_squared"] = sp.simplify(
        sp.diff(d["E"], u, 2).subs(origin)/2-16*momentum/5+sp.Rational(864, 25))
    out["odd_reduced_jet"] = leading_operator(jets["q_odd"])-48*u/5
    out["even_reduced_jet"] = sp.expand(leading_operator(jets["q_even"])+momentum*(2*delta+64*u**2)/5)
    x = jets["x"]
    out["odd_inner_polynomial"] = sp.factor(
        sp.diff(jets["H_odd"], x, 2)+80*jets["H_odd"]/(1+8*x**2)-48*x/5)
    out["even_inner_polynomial"] = sp.factor(
        sp.diff(jets["H_even"], x, 2)+80*jets["H_even"]/(1+8*x**2)
        +momentum*(2+64*x**2)/5)
    for n in range(9):
        out[f"full_inverse_monomial_{n}"] = sp.expand(leading_operator(inverse_monomial(0, n))-u**n)
    # At K=0 the same constant physical field on both metrics is exact.
    common = d["f_sum"]/d["f_sum"].subs(u, 0)
    out["K_zero_common_light"] = sp.simplify((sp.diff(common, u, 2)+d["A"]*common).subs(momentum, 0))
    out["K_zero_common_heavy"] = sp.simplify((d["E"]*common+u*d["f_cross"]*sp.diff(common, u)).subs(momentum, 0))
    return out
