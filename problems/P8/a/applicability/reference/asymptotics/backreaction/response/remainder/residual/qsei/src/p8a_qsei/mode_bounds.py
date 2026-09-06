"""Phase-resolved mode majorants; coincident stress bounds alone are not used.

All variables are dimensionless: x=eta/eta_star, k=k_physical*eta_star,
U=eta_star^2*U_physical. The immutable preparation has a flat past at x=1.
"""

from functools import cache

import sympy as sp
from p8a_remainder import geometry


@cache
def calibration():
    prior = geometry.calibration()
    b = list(map(sp.Rational, prior["potential_bounds"][:4]))
    cap = sp.Rational(prior["delta_bar"])
    span = sp.Integer(3)
    if span**2*cap*b[0] > 1:
        raise ValueError("The derived Volterra exponential bound is outside its domain")
    # |b_k|<=exp(delta*b0*T^2/2)<2. All higher naive derivative bounds
    # are delta*n_j, by transferring derivatives onto the flat-past product U*b.
    n = [None, 2*span*b[0]]
    n.append(span*(2*b[1]+cap*b[0]*n[1]))
    n.append(span*(2*b[2]+cap*(2*b[1]*n[1]+b[0]*n[2])))
    product = [2*b[j]+cap*sum(sp.binomial(j, ell)*b[ell]*n[j-ell]
                             for ell in range(j)) for j in range(4)]
    # |b_k-1|<=delta*l0/k and |b_k^(j)|<=delta*lj/k, j=1,2,3.
    inverse_frequency = [2*span*b[0]]
    inverse_frequency.extend((product[j-1]+span*product[j])/2 for j in range(1, 4))
    h = [sp.Integer(2), 4+cap*b[0]]
    h.append(cap*b[1]+2*h[0]*h[1])
    error = [inverse_frequency[j]+inverse_frequency[j+1]
             +sum(sp.binomial(j, ell)*h[ell]*inverse_frequency[j-ell]
                  for ell in range(j+1)) for j in range(3)]
    infrared_b_minus_one = span**2*b[0]
    infrared_error = (1+h[0])*infrared_b_minus_one+n[1]
    return {"delta_bar": cap, "potential_bounds_through_three": b,
            "active_conformal_span": span, "volterra_modulus_cap": sp.Integer(2),
            "volterra_exponent_cap": cap*b[0]*span**2/2,
            "naive_b_derivatives_per_delta": n[1:],
            "product_Ub_derivatives_per_delta": product,
            "b_inverse_frequency_constants": inverse_frequency,
            "conformal_Hc_derivative_caps": h,
            "ultraviolet_error_derivatives_per_delta": error,
            "infrared_b_minus_one_per_delta": infrared_b_minus_one,
            "infrared_error_per_delta": infrared_error}


def identities():
    x = sp.Symbol("x", real=True)
    k = sp.Symbol("k", positive=True)
    b, u = sp.Function("b")(x), sp.Function("U")(x)
    w = sp.exp(-sp.I*k*x)*b
    results = {"exact_phase_resolved_mode_equation": sp.simplify(
        sp.exp(sp.I*k*x)*(sp.diff(w, x, 2)+(k**2+u)*w)
        -sp.diff(b, x, 2)+2*sp.I*k*sp.diff(b, x)-u*b)}
    for order in range(4):
        results[f"product_rule_order_{order}"] = sp.expand(sp.diff(u*b, x, order)
            -sum(sp.binomial(order, ell)*sp.diff(u, x, ell)*sp.diff(b, x, order-ell)
                 for ell in range(order+1)))
    # The endpoint term in one IBP is essential. Here q stands for (U*b)^(j-1).
    q = sp.Function("q")(x)
    s = sp.Symbol("s", real=True)
    x0 = sp.Symbol("x0", real=True)
    integral = sp.Integral(sp.exp(2*sp.I*k*(x-s))*q.subs(x, s), (s, x0, x))
    by_parts = (-q+sp.exp(2*sp.I*k*(x-x0))*q.subs(x, x0)
                +sp.Integral(sp.exp(2*sp.I*k*(x-s))*sp.diff(q, x).subs(x, s), (s, x0, x)))/(2*sp.I*k)
    results["integration_by_parts_differential_identity"] = sp.simplify(
        sp.diff(integral-by_parts, x)-2*sp.I*k*(integral-by_parts))
    results["integration_by_parts_initial_value"] = sp.simplify((integral-by_parts).subs(x, x0).doit())
    return results
