"""Exact tensor polynomial division between two convergent subtractions."""

from functools import cache

import sympy as s

from . import tensor as ten

SIGMA, W, Q = s.symbols("sigma w q", positive=True)


@cache
def scalar_division():
    c = s.symbols("c0:4")
    polynomial = sum(c[j] * W**j for j in range(4))
    b = SIGMA + Q
    coefficients = tuple(
        sum(c[k] / b ** (j - k + 1) for k in range(j + 1)) for j in range(3)
    )
    truncated = sum(coefficients[j] * W**j for j in range(3))
    remainder = W**3 * polynomial.subs(W, b) / (b**3 * (b - W))
    return (
        polynomial,
        coefficients,
        s.factor(polynomial / (b - W) - truncated - remainder),
    )


@cache
def data():
    poly, coeff, residual = scalar_division()
    checks = {"generic_cubic_divided_difference": residual}
    for spin in (0, 2):
        cs = ten.centered_coefficients(spin)
        residue = sum((cs[j] * (ten.S + ten.Q) ** j for j in range(4)), s.zeros(9))
        checks[f"spin{spin}_full_tensor_spectral_residue"] = (
            residue - ten.S * ten.numerator(spin, ten.S)
        ).applyfunc(s.expand)
        checks[f"spin{spin}_zero_transfer_all_three_conversion_terms"] = sum(
            (cs[j].subs({p: 0 for p in ten.P}) for j in range(3)), s.zeros(9)
        )
    return {
        "centered_coefficients": "N_i(w-q,p)=sum_(k=0)^3 n_ik(p) w^k. Define J_in(q)=integral rho_i(s)/[s^3(s+q)^(n+1)]ds. Then A_r(p)=sum_i sum_(k=0)^r n_ik(p) J_i,r-k(q), for r=0,1,2.",
        "exact_conversion": "F(w-q,p)=A0(p)+w A1(p)+w^2 A2(p)+B6(w,p), where B6=w^3 integral sum_i rho_i(s)Q_i(s,p)/s^2 /[(s+q)^3(s+q-w)]ds. This is an exact identity of convergent analytic integrals off the cut.",
        "complete_tensor": "Polynomial division is applied to the entire spatial tensor numerator, not only TT, trace, or a fixed COM channel. The residue N_i(s)=s Q_i(s) yields precisely the actual spatial spectral tensor.",
        "time_operator": "In real time the finite conversion is A0 Gamma-A1 partial_t^2 Gamma+A2 partial_t^4 Gamma. Each coefficient is time local but generally nonpolynomial in spatial momentum.",
        "matching_boundary": "This is a conversion between explicitly defined nonlocal representatives. It neither fixes the remaining physical local tensor/contact polynomial nor equates the curved CD contact/endpoints to this flat benchmark.",
        "checks": checks,
        "gates": {
            "three_equal_time_coefficients": len(coeff) == 3,
            "cubic_numerator_not_dropped": s.degree(poly, W) == 3,
            "generic_division_exact": residual == 0,
        },
    }
