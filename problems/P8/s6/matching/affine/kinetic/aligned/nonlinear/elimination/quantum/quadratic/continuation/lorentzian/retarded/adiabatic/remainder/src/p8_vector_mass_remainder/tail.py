"""Integrable varied-reference readout tail after orders zero, two and four."""
from functools import cache

import sympy as sp
from p8_vector_mass_adiabatic import variation as prior
from p8_vector_state import wkb

from . import envelopes, source, tangent

add, multiply, scale = tangent.add, tangent.multiply, tangent.scale


def numerator(S, dS, R, dR, A, dA, B, dB, c, dc, r, adiabatic, absolute=False):
    sign = 1 if absolute else -1
    S2 = multiply(S, S)
    S3, S4 = multiply(S2, S), multiply(S2, S2)
    Z = add(scale(S2, 4*c**2), scale(multiply(R, S), 4*c), multiply(R, R))
    dZ = add(scale(S2, 8*c*dc), scale(multiply(S, dS), 8*c**2),
             scale(multiply(R, S), 4*dc), scale(multiply(dR, S), 4*c),
             scale(multiply(R, dS), 4*c), scale(multiply(R, dR), 2))
    N = add(scale(S4, 4*A), scale(S2, 4*B), scale(Z, A, 1))
    dN = add(scale(S4, 4*dA), scale(multiply(S3, dS), 16*A),
              scale(S2, 4*dB), scale(multiply(S, dS), 8*B),
              scale(Z, dA+sign*2*r*A, 1), scale(dZ, A, 1))
    return add(multiply(add(scale(N, r), dN), S), scale(multiply(N, dS), sign*3),
               scale(multiply(S4, adiabatic), sign*4))


@cache
def reference_tail(sector):
    source.kind(sector)
    data, s = tangent.reference(sector), source.data(sector)
    adiabatic = {j: source.project(value) for j, value in prior.data(sector)["readout_coefficients"].items()}
    exact = numerator(data["S"], data["dS"], data["R"], data["dR"],
                      s["A"], s["delta_A"], s["B"], s["delta_B"],
                      s["c1"], s["delta_c1"], s["r"], adiabatic)
    low = {j: source.clean(exact.get(j, 0)) for j in range(3)}
    reference = envelopes.reference(sector)
    upper = reference["polynomial_majorants"]

    def linear(value):
        result = source.linear_bound(source.clean(value))
        if any(value != 0 for value in result["reconstructions"]):
            raise ValueError("A source-coefficient reconstruction failed")
        return result["absolute_upper"]

    def base(value):
        result = wkb.box_bound(sp.factor(value))
        if result["reconstruction"] != 0:
            raise ValueError("A reference-coefficient reconstruction failed")
        return result["absolute_upper"]

    majorant = numerator(upper["S"], upper["dS"], upper["R"], upper["dR"],
                         base(s["A"]), linear(s["delta_A"]), base(s["B"]), linear(s["delta_B"]),
                         base(s["c1"]), linear(s["delta_c1"]), linear(s["r"]),
                         {j: linear(value) for j, value in adiabatic.items()}, absolute=True)
    mass = wkb.MASS_TIME_MIN
    # delta Q_ref-delta Q_ad = omega*numerator/(16*S^4), and S>=1/2.
    bound = sp.ceiling(sum(value/mass**(2*(j-3)) for j, value in majorant.items() if j >= 3))
    return {"low_tail_numerator_coefficients": low,
            "tail_numerator_majorants": majorant,
            "varied_reference_tail_over_inverse_frequency_fifth_upper": bound,
            "all_majorants_nonnegative": all(value >= 0 for value in majorant.values())}


@cache
def algebra_checks():
    S, R, t, A, B, c, r = sp.symbols("S R t A B c r", nonzero=True)
    dS, dR, dA, dB, dc = sp.symbols("dS dR dA dB dc")
    F = A*S+B/S+A*t*(c+R/(2*S))**2/S
    variation = r*F+sum(sp.diff(F, variable)*change for variable, change in
                       ((S, dS), (R, dR), (A, dA), (B, dB), (c, dc), (t, -2*r*t)))
    polynomial = numerator({0: S}, {0: dS}, {0: R}, {0: dR}, A, dA, B, dB, c, dc, r, {})
    represented = sum(value*t**j for j, value in polynomial.items())/(4*S**4)
    return {"exact_varied_readout_denominator": sp.factor(variation-represented)}
