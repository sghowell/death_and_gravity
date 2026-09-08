"""Continuous source-linear bounds for the eighth-order reference response."""
from functools import cache

import sympy as sp
from p8_vector_regularity import frequency
from p8_vector_state import wkb

from . import source, tangent

add, multiply, scale = tangent.add, tangent.multiply, tangent.scale


@cache
def reference(sector):
    source.kind(sector)
    data = tangent.reference(sector)
    upper, reconstructions = {}, []
    for key in ("S", "R", "T"):
        row = {}
        for j, coefficient in data[key].items():
            bound = wkb.box_bound(sp.factor(coefficient))
            row[j] = bound["absolute_upper"]
            reconstructions.append(bound["reconstruction"])
        upper[key] = row
    for key in ("dS", "dR", "dT"):
        row = {}
        for j, coefficient in data[key].items():
            bound = source.linear_bound(source.clean(coefficient))
            row[j] = bound["absolute_upper"]
            reconstructions.extend(bound["reconstructions"])
        upper[key] = row
    s = source.data(sector)
    rb = source.linear_bound(s["r"])["absolute_upper"]
    S, dS, R, dR, T, dT = [upper[key] for key in ("S", "dS", "R", "dR", "T", "dT")]
    S2, SdS = multiply(S, S), multiply(S, dS)
    Q = {j-1: coefficient for j, coefficient in S2.items() if j >= 2}
    dQ = add(scale(SdS, 2, -1), scale(S2, 2*rb, -1))
    dQ = {j: value for j, value in dQ.items() if j >= 1}
    N = add(multiply(Q, S2), multiply(T, S), scale(multiply(R, R), sp.Rational(3, 4)))
    dN = add(multiply(dQ, S2), scale(multiply(Q, SdS), 2), multiply(dT, S),
              multiply(T, dS), scale(multiply(R, dR), sp.Rational(3, 2)))
    numerator = add(multiply(dN, S), scale(multiply(N, dS), 2))
    mass = wkb.MASS_TIME_MIN
    residual = sp.ceiling(8*sum(value/mass**(2*(j-4)) for j, value in numerator.items() if j >= 4))
    K = data["frequency_variation_coefficients"]
    frequency_bound = sum(source.linear_bound(source.clean(value))["absolute_upper"]/mass**(2*j)
                          for j, value in K.items())
    slope_bound = sum(source.linear_bound(source.time(value)+(1-2*j)*source.lam*value)["absolute_upper"]/mass**(2*j)
                      for j, value in K.items())
    delta_d = source.linear_bound(s["delta_d"])["absolute_upper"]
    delta_L = 2*slope_bound+8*frequency_bound
    old = frequency.reference("transverse" if sector == "T" else "longitudinal")
    return {"polynomial_majorants": upper, "box_reconstructions": reconstructions,
            "varied_residual_numerator_majorants": numerator,
            "delta_residual_over_inverse_frequency_eighth_upper": residual,
            "delta_W_over_frequency_upper": frequency_bound,
            "delta_W_prime_over_frequency_upper": slope_bound,
            "delta_log_W_rate_upper": delta_L,
            "delta_c_upper": delta_d+delta_L/2,
            "base_residual_over_inverse_frequency_eighth_upper": old["residual_over_inverse_frequency_power_upper"],
            "base_reference_checks": old["proof_checks"],
            "source_A_variation_upper": source.linear_bound(s["delta_A"])["absolute_upper"],
            "source_B_variation_upper": source.linear_bound(s["delta_B"])["absolute_upper"],
            "source_log_frequency_variation_upper": rb}
