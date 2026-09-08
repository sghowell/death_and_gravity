"""Eighth-order reference variation from the exact linearized Riccati equation."""
from functools import cache

import sympy as sp
from p8_vector_metric_local import variation as prior
from p8_vector_regularity import spectral

from . import source

add, multiply, scale = spectral.add, spectral.multiply, spectral.scale


def truncated(polynomial, maximum):
    return {j: source.clean(value) for j, value in polynomial.items() if j <= maximum}


def time(polynomial):
    return {j: source.clean(source.time(value)-2*j*source.lam*value)
            for j, value in polynomial.items()}


def inverse(S, maximum):
    out = {0: sp.Integer(1)}
    for order in range(1, maximum+1):
        out[order] = source.clean(-sum(S.get(j, 0)*out[order-j] for j in range(1, order+1)))
    return out


def coefficient(sector, order):
    source.kind(sector)
    if type(order) is not int or not 1 <= order <= 4:
        raise ValueError("Require native varied reference coefficient order 1..4")
    return _coefficient(sector, order)


@cache
def _coefficient(sector, order):
    data = source.data(sector)
    S = {j: source.baseline(sector, j) for j in range(order+1)}
    dS = {j: coefficient(sector, j) for j in range(1, order)}
    inv = inverse(S, order-1)
    L = add({0: source.lam}, truncated(multiply(time(S), inv), order-1))
    dL = add({0: data["delta_lambda"]}, time(truncated(multiply(dS, inv), order-1)))
    rest = (-2*multiply(S, dS).get(order, 0)
            -2*data["r"]*multiply(S, S).get(order, 0)
            -(data["delta_U"] if order == 1 else 0)
            -time(dL).get(order-1, 0)/2
            +multiply(L, dL).get(order-1, 0)/2)
    return source.clean(rest/2)


def varied_P(sector, order):
    return source.clean(coefficient(sector, order)+2*order*source.data(sector)["r"]*source.baseline(sector, order))


@cache
def low_checks():
    out = {}
    for sector in ("T", "L"):
        for order, key in ((1, "delta_P2"), (2, "delta_P4")):
            out[sector+"_frozen_dimensional_variation_order_"+str(2*order)] = source.clean(
                varied_P(sector, order)-source.project(prior.data(sector)[key]))
    return out


@cache
def reference(sector):
    source.kind(sector)
    data = source.data(sector)
    S = {j: source.baseline(sector, j) for j in range(5)}
    dS = {j: coefficient(sector, j) for j in range(1, 5)}
    R, dR = time(S), time(dS)
    T = add(scale(time(R), -sp.Rational(1, 2)), scale(R, source.lam/2))
    dT = add(scale(time(dR), -sp.Rational(1, 2)),
             scale(R, data["delta_lambda"]/2), scale(dR, source.lam/2))
    S2, SdS = multiply(S, S), multiply(S, dS)
    Q = {j-1: -value for j, value in S2.items() if j >= 2}
    dQ = add(scale(SdS, -2, -1), scale(add({0: sp.Integer(1)}, scale(S2, -1)), 2*data["r"], -1),
             {0: 2*varied_P(sector, 1)})
    dQ = {j: source.clean(value) for j, value in dQ.items() if source.clean(value) != 0}
    N = add(multiply(Q, S2), multiply(T, S), scale(multiply(R, R), sp.Rational(3, 4)))
    dN = add(multiply(dQ, S2), scale(multiply(Q, SdS), 2), multiply(dT, S),
              multiply(T, dS), scale(multiply(R, dR), sp.Rational(3, 2)))
    residual_numerator = add(multiply(dN, S), scale(multiply(N, dS), -2))
    low = {j: source.clean(residual_numerator.get(j, 0)) for j in range(4)}
    highest = {name: source.clean(residual_numerator.get(4, 0)).coeff(row[10])\
               for name, row in (("N", source.n), ("Z", source.v))}
    return {"S": S, "dS": dS, "R": R, "dR": dR, "T": T, "dT": dT,
            "Q": Q, "dQ": dQ, "N": N, "dN": dN,
            "varied_residual_low_coefficients": low,
            "tenth_source_derivative_fixtures": {name: sp.factor(value.subs({source.u: 0, source.z: 0}))\
                                                 for name, value in highest.items()},
            "frequency_variation_coefficients": add(scale(S, data["r"]), dS)}
