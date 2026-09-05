"""Quartic invariant reduction, with exact symbolic jets and uniform bounds."""

from functools import cache

import sympy as sp
from p8_s5.lapse_series import stationary_series

from . import nonlinear as model

x, y = sp.symbols("x y", real=True)


@cache
def symbolic_series():
    jets = {key: sp.symbols(f"{key}0:5", real=True) for key in ("A", "B", "C", "L", "D", "E", "M", "Z")}
    sigma, rho, eta, shear2, z = model.VARIABLES
    L = tuple(jets["B"][k]*sigma+jets["C"][k]*rho+jets["L"][k]*eta for k in range(4))
    Q = tuple(jets["D"][k]*sigma**2+jets["E"][k]*shear2+jets["M"][k]*eta**2+jets["Z"][k]*z for k in range(3))
    return jets, stationary_series(jets["A"], L, Q)


def coefficients(expression):
    return {",".join(map(str, powers)): sp.factor(value)
            for powers, value in sp.Poly(expression, *model.VARIABLES).terms() if value != 0}


def compact_even(expression):
    """Map an even rational function regular after compactification to x.

    This verifies the numerator has no poles except d=1+u^2, and is polynomial
    on x=u/sqrt(d). It is not an unevaluated limit or a sampled tail bound.
    """
    numerator, denominator = sp.fraction(sp.cancel(expression))
    constant, factors = sp.factor_list(denominator, model.u)
    n = 0
    for base, power in factors:
        if sp.expand(base-model.d) != 0:
            raise ValueError(f"Unexpected compact denominator {base}")
        n += power
    result = 0
    for (power,), value in sp.Poly(numerator, model.u).terms():
        if power % 2 or power > 2*n:
            raise ValueError("Not an even compact polynomial")
        result += value/constant*x**power*(1-x**2)**(n-power//2)
    return sp.expand(result)


def reduce_y(expression):
    """Use y^2=1-x^2 with y>=0, retaining at most a single y factor."""
    return sp.expand(sum(value*y**(power[0] % 2)*(1-x**2)**(power[0]//2)
                         for power, value in sp.Poly(sp.expand(expression), y).terms()))


@cache
def compact_jets():
    """Local fixed units ell=tau*sqrt(d); no time-dependent canonical claim."""
    f = model.functions()
    poly = sp.Poly(f["F"], model.X)
    generic = model.generic_lapse_jets()
    values = ((1-x**2)**3, -6*x*(1-x**2)**3, 6*(8*x**2-1)*(1-x**2)**3,
              4*x, y**11/10, *(compact_even(model.d*poly.nth(k)) for k in range(3)))
    mapping = dict(zip(generic["symbols"], values))
    return {key: tuple(reduce_y(expr.subs(mapping)) for expr in row)
            for key, row in generic["jets"].items()}


def polynomial_bound(expression):
    """l1 coefficient bound on |x|<=1, 0<=y<=1, including both tails."""
    return sum(abs(c) for c in sp.Poly(expression, x, y).coeffs())


@cache
def compact_checks():
    generic = model.generic_lapse_jets()
    R, Rd, Rdd, HH, l, f0, f1, f2 = generic["symbols"]
    c = compact_jets()
    u, d = model.u, model.d
    # Every generic jet has the required dimension under ell rescaling.
    ell = sp.Symbol("ell", positive=True)
    scale = {R: R, Rd: ell*Rd, Rdd: ell**2*Rdd, HH: ell*HH, l: ell*l,
             f0: ell**2*f0, f1: ell**2*f1, f2: ell**2*f2}
    degrees = {"A": 2, "B": 1, "C": 0, "L": 1, "D": 0, "E": 0, "M": 0, "Z": 0}
    checks = {f"jet_scaling_{key}_{k}": sp.expand(expr.subs(scale, simultaneous=True)-ell**degrees[key]*expr)
              for key, row in generic["jets"].items() for k, expr in enumerate(row)}
    checks.update({
        "compact_H": sp.simplify(4*u/sp.sqrt(d)-sp.sqrt(d)*model.H),
        "compact_Rd": sp.simplify(-6*u/d**sp.Rational(7, 2)-sp.sqrt(d)*sp.diff(1/d**3, u)),
        "compact_Rdd": sp.cancel(sp.diff(1/d**3, u, 2)-6*(7*u**2-1)/d**5),
        "compact_J": sp.expand(c["A"][2]+2*compact_even(d*model.functions()["J"])),
        "compact_background_lapse": c["A"][1],
    })
    poly = sp.Poly(model.functions()["F"], model.X)
    checks.update({f"compact_F_{k}": sp.cancel(compact_even(d*poly.nth(k)).subs(x, u/sp.sqrt(d))-d*poly.nth(k))
                   for k in range(3)})
    return checks


def coefficient_bound(expression, jets_bound):
    """Exact rational bound: the sole variable denominator is A2, |A2|>1/5."""
    jet_symbols, _ = symbolic_series()
    a2 = jet_symbols["A"][2]
    numerator, denominator = sp.fraction(sp.cancel(expression))
    den_poly = sp.Poly(denominator, a2)
    if len(den_poly.terms()) != 1 or den_poly.LC().free_symbols:
        raise ValueError("Unexpected stationary denominator")
    inverse = sp.Rational(5)**den_poly.degree()/abs(den_poly.LC())
    symbols = sorted(numerator.free_symbols, key=str)
    if not symbols:
        return abs(numerator)*inverse
    polynomial = sp.Poly(numerator, *symbols)
    return sp.factor(inverse*sum(abs(value)*sp.prod(jets_bound[s]**p for s, p in zip(symbols, powers))
                                for powers, value in polynomial.terms()))


@cache
def coefficient_report():
    symbolic, result = symbolic_series()
    compact = compact_jets()
    bounds = {symbolic[key][k]: polynomial_bound(compact[key][k]) for key in symbolic for k in range(5)}
    bounce = {symbolic[key][k]: compact[key][k].subs({x: 0, y: 1}) for key in symbolic for k in range(5)}
    report = {}
    for name, orders in result.items():
        report[name] = [{powers: {"jet_expression": str(value),
                                  "at_bounce": str(sp.factor(value.subs(bounce))),
                                  "uniform_absolute_upper": str(coefficient_bound(value, bounds))}
                        for powers, value in coefficients(expr).items()} for expr in orders]
    return report
