"""Exact weighted integration gains and the actual shared-history response."""

import sympy as sp
from p8a_existence.mode_lipschitz import rational


def nonnegative(value, *, positive=False):
    if not isinstance(positive, bool):
        raise TypeError("positive must be a boolean")
    value = rational(value)
    if value < 0 or (positive and value == 0):
        raise ValueError("an exact positive value is required" if positive else "an exact nonnegative value is required")
    return value


def metric_gains(sigma):
    """Per weighted ||Delta X||, conditional on the stated positive a,h box."""
    alpha = 1/nonnegative(sigma, positive=True)
    return {"u": alpha, "h": alpha**2, "log_a": alpha**3,
            "a_squared": 18*alpha**3, "inverse_a_squared": alpha**3/2,
            "d": alpha**3/2, "dprime": alpha**2/2}


def prefix_response(derivative_cap, history, length, log_ceiling=20):
    """A10 actual nonlinear R' Lipschitz bound in the future weighted norm.

The exponential weight does not enlarge this causal prefix estimate.
The logarithm ceiling is verified using e>8/3, not accepted as an input
assumption. Both potentials have the same original past and U(0)=0.
    """
    m = nonnegative(derivative_cap)
    t, ell = (nonnegative(value, positive=True) for value in (history, length))
    n = nonnegative(log_ceiling, positive=True)
    if n.q != 1 or ell > t:
        raise ValueError("an integer log ceiling and length<=history are required")
    if sp.Rational(8, 3)**n < t/ell:
        raise ValueError("the proposed logarithm ceiling is not certified")
    z = m*t**3
    if z > sp.Rational(1, 4):
        raise ValueError("the A10 strength condition M T^3<=1/4 is required")
    quadratic = m*t*ell**2*(sp.Rational(5, 4)+n/2)
    higher = 36*m**2*t**5*ell
    return {"strength": z, "log_ratio_upper": n, "exponential_upper": sp.Integer(2),
            "quadratic": quadratic, "higher": higher, "total": quadratic+higher}


def identities():
    t, ell, sigma, history = sp.symbols("t ell sigma T", positive=True)
    s = sp.Symbol("s", nonnegative=True)
    weight = sp.integrate(sp.exp(-sigma*(t-s)), (s, 0, t))
    quadratic = ell**2*(sp.Rational(5, 4)+sp.log(history/ell)/2)
    return {
        "weighted_primitive_kernel": sp.simplify(weight-(1-sp.exp(-sigma*t))/sigma),
        "shared_quadratic_prefix_derivative": sp.simplify(
            sp.diff(quadratic, ell)-ell*(2+sp.log(history/ell))),
        "shared_quadratic_zero_prefix_limit": sp.limit(quadratic, ell, 0, dir="+"),
        "damped_metric_difference_has_two_gains": (1/sigma)*(1/sigma)-1/sigma**2,
    }
