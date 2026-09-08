"""An explicit smooth clock-norm window with finite derivative bounds."""
from functools import lru_cache

import sympy as sp
from p8_vector_hadamard import cutoffs


def cutoff(clock_norm):
    if type(clock_norm) is int:
        clock_norm = sp.Integer(clock_norm)
    return 1-cutoffs.turn_on(64*(clock_norm+1)**2)


@lru_cache(maxsize=None, typed=True)
def bump_bound(order):
    if type(order) is not int or not 0 <= order <= 5:
        raise ValueError("Require native derivative order 0..5")
    y, polynomial = cutoffs.bump_derivative_polynomial(order)
    return sum(abs(coefficient)*sp.factorial(power) for (power,), coefficient in sp.Poly(polynomial, y).terms())


@lru_cache(maxsize=None, typed=True)
def inverse_denominator_bound(order):
    bump_bound(order)
    if order == 0:
        return sp.Integer(8)
    return 16*sum(sp.binomial(order, k)*bump_bound(k)*inverse_denominator_bound(order-k)
                  for k in range(1, order+1))


@lru_cache(maxsize=None, typed=True)
def step_bound(order):
    bump_bound(order)
    if order == 0:
        return sp.Integer(1)
    return sum(sp.binomial(order, k)*bump_bound(k)*inverse_denominator_bound(order-k)
               for k in range(order+1))


@lru_cache(maxsize=None, typed=True)
def derivative_bound(order):
    bump_bound(order)
    if order == 0:
        return sp.Integer(1)
    return sum(sp.factorial(order)*step_bound(order-j)*32**(order-2*j)*64**j
               /(sp.factorial(order-2*j)*sp.factorial(j)) for j in range(order//2+1))


def proof_checks():
    exp_two_upper = sum(sp.Rational(2**n, sp.factorial(n)) for n in range(7))+sp.Rational(2**7, sp.factorial(7))/(1-sp.Rational(1, 4))
    return {"bump_denominator_above_one_eighth": bool(exp_two_upper < 8),
            "original_clock_tube_inside_flat_window": bool(64*sp.Rational(1, 10)**2 < 1),
            "support_inside_one_quarter_clock_norm_displacement": bool(sp.Rational(2, 64) < sp.Rational(1, 4)**2),
            "vacuum_clock_norm_in_open_zero_region": bool(sp.Integer(64) > 2),
            "step_derivative_bounds_are_finite_positive_integers": all(derivative_bound(n).is_Integer and derivative_bound(n) > 0 for n in range(6))}


def chain_rule_checks():
    s = sp.Symbol("clock_norm_displacement", real=True)
    f = sp.Function("step")
    out = {}
    for n in range(1, 6):
        target = sum(sp.factorial(n)/(sp.factorial(n-2*j)*sp.factorial(j))
                     *sp.Subs(sp.Derivative(f(sp.Symbol("z")), (sp.Symbol("z"), n-j)), sp.Symbol("z"), 64*s**2)
                     *(128*s)**(n-2*j)*64**j for j in range(n//2+1))
        out["quadratic_window_chain_rule_"+str(n)] = sp.simplify(sp.diff(f(64*s**2), s, n)-target)
    return out
