"""Exact eighth-order reference residual and differentiated envelopes."""
from functools import lru_cache

import sympy as sp
from p8_vector_hadamard import series
from p8_vector_state import wkb

from . import spectral


def envelope(value):
    value_bound = wkb.box_bound(value)
    derivative_bound = wkb.box_bound(sp.factor(spectral.D0(value)))
    return (value_bound["absolute_upper"], derivative_bound["absolute_upper"]), (
        value_bound["reconstruction"], derivative_bound["reconstruction"])


def pair_add(*polynomials):
    out = {}
    for polynomial in polynomials:
        for degree, (value, derivative) in polynomial.items():
            old = out.get(degree, (sp.Integer(0), sp.Integer(0)))
            out[degree] = (old[0]+value, old[1]+derivative)
    return out


def pair_multiply(left, right):
    out = {}
    for j, (a, da) in left.items():
        for k, (b, db) in right.items():
            degree = j+k
            old = out.get(degree, (sp.Integer(0), sp.Integer(0)))
            out[degree] = (old[0]+a*b, old[1]+da*b+a*db)
    return out


def pair_scale(polynomial, factor):
    return {degree: (factor*value, factor*derivative)
            for degree, (value, derivative) in polynomial.items()}


@lru_cache(maxsize=None, typed=True)
def reference(kind, coefficient_order=4):
    if kind not in ("transverse", "longitudinal") or type(coefficient_order) is not int or coefficient_order not in (2, 4):
        raise ValueError("Require an actual vector mode and native WKB coefficient order 2 or 4")
    lam, mass = wkb.background()["lambda"], wkb.MASS_TIME_MIN
    S, R, T = {0: sp.Integer(1)}, {}, {}
    Sb, Rb, Tb = {0: (sp.Integer(1), sp.Integer(0))}, {}, {}
    reconstructions = []
    for order in range(1, coefficient_order+1):
        P = series.coefficient(kind, order)
        B = sp.factor(spectral.D0(P)-2*order*lam*P)
        E = sp.factor(-spectral.D0(B)/2+(sp.Rational(1, 2)+order)*lam*B)
        S[order], R[order], T[order] = P, B, E
        for exact, upper in ((S, Sb), (R, Rb), (T, Tb)):
            upper[order], checks = envelope(exact[order])
            reconstructions.extend(checks)
    S2 = spectral.multiply(S, S)
    Q = {degree-1: -value for degree, value in S2.items() if degree >= 2}
    numerator = spectral.add(spectral.multiply(Q, S2), spectral.multiply(T, S),
                             spectral.scale(spectral.multiply(R, R), sp.Rational(3, 4)))
    S2b = pair_multiply(Sb, Sb)
    Qb = {degree-1: pair for degree, pair in S2b.items() if degree >= 2}
    upper = pair_add(pair_multiply(Qb, S2b), pair_multiply(Tb, Sb),
                     pair_scale(pair_multiply(Rb, Rb), sp.Rational(3, 4)))
    low = {degree: sp.factor(numerator.get(degree, 0)) for degree in range(coefficient_order)}
    value_sum = sum(value/mass**(2*(degree-coefficient_order))
                    for degree, (value, _) in upper.items() if degree >= coefficient_order)
    derivative_sum = sum((derivative+4*degree*value)/mass**(2*(degree-coefficient_order))
                         for degree, (value, derivative) in upper.items() if degree >= coefficient_order)
    S_deviation = sum(pair[0]/mass**(2*degree) for degree, pair in Sb.items() if degree > 0)
    R_upper = sum(pair[0]/mass**(2*degree) for degree, pair in Rb.items())
    C, D = sp.ceiling(4*value_sum), sp.ceiling(4*derivative_sum+16*value_sum*R_upper)
    return {"coefficient_order": coefficient_order,
            "low_residual_numerator_coefficients": low,
            "box_reconstructions": reconstructions,
            "residual_numerator_coefficient_majorants": upper,
            "frequency_ratio_deviation_upper": S_deviation,
            "frequency_ratio_derivative_upper": R_upper,
            "residual_over_inverse_frequency_power_upper": C,
            "residual_derivative_over_inverse_frequency_power_upper": D,
            "reference_logarithmic_rate_upper": 2+2*R_upper,
            "oscillatory_evolution_beta_upper": 8*D+300*C,
            "proof_checks": {
                "all_low_residual_coefficients_vanish": all(value == 0 for value in low.values()),
                "frequency_ratio_between_one_half_and_three_halves": bool(S_deviation < sp.Rational(1, 2)),
                "logarithmic_rate_below_four": bool(2+2*R_upper < 4),
                "phase_frequency_above_one_quarter_omega": bool(C/mass**(2*coefficient_order+2) < sp.Rational(1, 4)),
                "phase_derivative_below_seven_omega": bool((D+4*C)/mass**(2*coefficient_order+2) < 1),
                "evolution_norm_integral_below_one_quarter": bool(2*C/mass**(2*coefficient_order+1) < sp.Rational(1, 4)),
                "feedback_term_below_claimed_allowance": bool(C/mass**(2*coefficient_order+1) < sp.Rational(1, 4)),
                "majorant_pairs_are_nonnegative": all(value >= 0 and derivative >= 0 for value, derivative in upper.values())}}


def algebra_checks():
    S, R, DR, lam, Dlam, U, t, P1 = sp.symbols("S R DR lam Dlam U t P1", nonzero=True)
    rate = lam+R/S
    D_rate = Dlam+DR/S-R**2/S**2
    residual = (1-S**2)/t-U-D_rate/2+rate**2/4
    reduced = (1-S**2)/t+2*P1+(-DR/2+lam*R/2)/S+3*R**2/(4*S**2)
    return {"exact_residual_denominator_identity": sp.factor(
        (residual-reduced).subs(P1, -U/2-Dlam/4+lam**2/8))}
