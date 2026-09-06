"""C1-potential bounds for the exact two-direction scattering split.

No U'' or U''' norm enters these estimates. The leading backward-history
term is controlled by Parseval in momentum, not a pointwise decay claim.
"""

from functools import cache

import sympy as sp
from p8a_existence.mode_lipschitz import rational
from p8a_preparation import bounds as prior


def nonnegative(value, *, positive=False):
    if not isinstance(positive, bool):
        raise TypeError("positive must be a boolean")
    result = rational(value)
    if result < 0 or (positive and result == 0):
        raise ValueError("a positive value is required" if positive else "a nonnegative value is required")
    return result


def majorants(history, potential_cap, derivative_cap, h0, h1, h2):
    """Return bounds conditional on smooth free-past modes and stated caps.

The frequency split is k=1. h0,h1,h2 bound the conformal Hubble function
and its first two derivatives on the SAMPLING interval only. Potential
and derivative caps hold on the full original history of length <=T.
    """
    t = nonnegative(history, positive=True)
    b, m, h0, h1, h2 = map(nonnegative, (potential_cap, derivative_cap, h0, h1, h2))
    if b*t**2 > 1:
        raise ValueError("the stated modulus cap requires B T^2 <= 1")
    q0, q1 = 2*b, 2*m+2*t*b**2
    a0, a1, a2 = 2*b*t, q0, q1
    forward = [(1+h0)*a0/2,
               ((1+h0)*a1+h1*a0)/2,
               ((1+h0)*a2+2*h1*a1+h2*a0)/2]
    backward_coefficient = (1+h0)/4
    local = [backward_coefficient*q0,
             backward_coefficient*q1+h1*q0/4]
    bminus1 = 2*t*b
    bprime = b+t*m+t**2*b**2
    derivative_remainder = m*bminus1+b*bprime
    infrared = (1+h0)*b*t**2+2*b*t
    return {"history": t, "potential_cap": b, "derivative_cap": m,
            "hubble_caps": [h0, h1, h2], "modulus_cap": sp.Integer(2),
            "modulus_exponent": b*t**2/2,
            "q_and_qprime": [q0, q1], "A_and_two_derivatives": [a0, a1, a2],
            "forward_two_derivatives": forward,
            "local_one_derivative": local,
            "backward_coefficient": backward_coefficient,
            "bminus1_inverse_k": bminus1, "bprime_inverse_k": bprime,
            "qprime_minus_uprime_inverse_k": derivative_remainder,
            "infrared_error": infrared}


@cache
def calibration():
    data = prior.calibration()
    geometry = data["geometry"]
    t, m = geometry["history"], geometry["history_uprime_cap"]
    h0, h1, h2 = geometry["h_max"], sp.Rational(1, 3), sp.Rational(1, 3)
    result = majorants(t, m*t, m, h0, h1, h2)
    margins = {"local_hprime": h1-h0**2-geometry["u_cap"],
               "local_hsecond": h2-m-2*h0*(h0**2+geometry["u_cap"]),
               "strict_modulus_exponent": sp.Rational(1, 2)-result["modulus_exponent"]}
    if any(value <= 0 for value in margins.values()):
        raise ValueError("the C1 scattering geometric caps failed")
    return {**result, "strict_margins": margins}


def identities():
    x = sp.Symbol("x", real=True)
    k = sp.Symbol("k", positive=True)
    h, aa, bb, q, jj = sp.symbols("h A B q J", complex=True)
    phase = sp.exp(2*sp.I*k*x)
    bminus1 = -(phase*bb-aa)/(2*sp.I*k)
    bprime = -phase*bb
    error = -sp.I*k*bminus1+bprime-h*bminus1
    forward, backward = -sp.Rational(1, 2)-h/(2*sp.I*k), -sp.Rational(1, 2)+h/(2*sp.I*k)
    after_ibp = forward*aa-backward*q/(2*sp.I*k)+backward*phase*jj/(2*sp.I*k)
    u, up, b, bp = sp.symbols("u up b bp", complex=True)
    qfun = sp.Function("q")(x)
    s = sp.Symbol("s", real=True)
    integral = sp.Integral(sp.exp(-2*sp.I*k*s)*qfun.subs(x, s), (s, 0, x))
    ibp = (-sp.exp(-2*sp.I*k*x)*qfun+qfun.subs(x, 0)
           +sp.Integral(sp.exp(-2*sp.I*k*s)*sp.diff(qfun, x).subs(x, s), (s, 0, x)))/(2*sp.I*k)
    return {
        "exact_forward_backward_split": sp.simplify(error-forward*aa-backward*phase*bb),
        "exact_backward_local_history_split": sp.simplify(
            error.subs(bb, (-sp.exp(-2*sp.I*k*x)*q+jj)/(2*sp.I*k))-after_ibp),
        "qprime_remainder_decomposition": sp.expand(up*b+u*bp-up-(up*(b-1)+u*bp)),
        "oscillatory_ibp_derivative": sp.simplify(sp.diff(integral-ibp, x)),
        "oscillatory_ibp_initial_endpoint": sp.simplify((integral-ibp).subs(x, 0).doit()),
        "real_history_half_parseval_factor": sp.pi*2/2/2-sp.pi/2,
    }


def controls():
    q0, k = sp.symbols("q0 k", positive=True)
    return {"omitted_original_endpoint": q0/(2*sp.I*k),
            "naive_backward_momentum_weight": "integral_1^infinity dk/k diverges",
            "first_Parseval_is_full_for_complex_products": True,
            "no_quantitative_u_second_or_third_derivative": True}
