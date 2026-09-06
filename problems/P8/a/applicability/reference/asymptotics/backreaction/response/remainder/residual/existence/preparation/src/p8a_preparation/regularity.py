"""Tame highest-jet bounds for the smooth, flat-start Picard construction.

These exact checks accompany the induction in notes/regularity.md. They do
not say that the isolated logarithmic inverse is an endomorphism of C1.
"""

import sympy as sp
from p8a_existence.mode_lipschitz import exponential_upper, rational


def highest_jet_bound(derivative_cap, history, local_cap, inverse_cap):
    """Bound the same highest-jet block at every finite derivative order.

The quadratic polarized log kernel has mass T^2/2. Higher Dyson terms are
bounded by the larger A.10 sum with n^2 in place of n. All jets at the
original free-past join vanish, and future perturbations have a common
initial zero neighborhood. Without those hypotheses endpoint terms occur.
    """
    m, t, d, k = map(rational, (derivative_cap, history, local_cap, inverse_cap))
    if m < 0 or t <= 0 or d < 0 or k < 0:
        raise ValueError("nonnegative caps and a positive history are required")
    z = m*t**3
    if z > sp.Rational(1, 4):
        raise ValueError("the rational Dyson bound requires M T^3 <= 1/4")
    nonlinear = z+18*z**2*exponential_upper(2*z)
    result = k*(d+nonlinear)
    if result >= 1:
        raise ValueError("the highest-jet block must be a strict contraction")
    return {"strength": z, "quadratic_Frechet_cap": z,
            "higher_Frechet_cap": nonlinear-z, "Frechet_cap": nonlinear,
            "highest_jet_contraction": result, "highest_jet_margin": 1-result}


def insertion_counts(order, derivative_order):
    """Leibniz allocation count and its exact top-jet coefficient.

For an n-linear Dyson term, distributing m derivatives over n potential
factors has multinomial total n^m. There are n top-jet terms, with every
other term involving strictly lower derivatives when m>=1.
    """
    n, m = map(rational, (order, derivative_order))
    if n.q != 1 or m.q != 1 or n < 2 or m < 1:
        raise ValueError("integer order>=2 and derivative_order>=1 required")
    return {"total_multinomial_weight": n**m, "top_jet_terms": n,
            "lower_jet_weight": n**m-n}


def frechet_continuity_bound(derivative_cap, history):
    """Bound ||DR[u]-DR[v]|| per ||u-v|| for ||u||,||v||<=M T."""
    m, t = map(rational, (derivative_cap, history))
    if m < 0 or t <= 0 or m*t**3 > sp.Rational(1, 4):
        raise ValueError("nonnegative M, positive T and M T^3<=1/4 required")
    return t**2*(1+6*(exponential_upper(2*m*t**3)-1))


def identities():
    n = sp.Symbol("n", positive=True, integer=True)
    b, t = sp.symbols("B T", positive=True)
    k, d, z, e = sp.symbols("K d z E", nonnegative=True)
    u, v, up, vp = sp.symbols("u v up vp", real=True)
    return {
        "polarized_quadratic_kernel_mass": 2*b*t**2/2-b*t**2,
        "larger_Dyson_coefficient_gap": sp.expand(n**2-n-n*(n-1)),
        "highest_block_norm_allocation": sp.expand(k*d+k*z+18*k*z**2*e-k*(d+z+18*z**2*e)),
        "second_order_derivative_transfer": sp.expand(up*v+u*vp-(up*v+u*vp)),
    }
