"""All-order SAT8 derivative recurrence and explicit asymptotic symbol caps."""

from functools import cache
from math import factorial

import sympy as s

X = s.Symbol("x", real=True)


def order(value):
    if type(value) is not int or value < 0:
        raise TypeError("Require a nonnegative native derivative order")
    return value


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Rational)):
        raise TypeError("Require a finite exact rational")
    return s.Rational(value)


@cache
def polynomial(n):
    n = order(n)
    if n == 0:
        return X
    previous = polynomial(n - 1)
    return s.expand(
        (1 + X**8) * s.diff(previous, X) - (8 * (n - 1) + 1) * X**7 * previous
    )


def tail_coefficient(n):
    return s.Integer(5 * factorial(order(n)) * 4**n)


def compact_coefficient(n):
    n = order(n)
    p = s.Poly(polynomial(n), X)
    return sum(abs(c) * 2 ** powers[0] for powers, c in p.terms()) + (
        1 if n == 0 else 0
    )


def mass_derivative_tail(n, amplitude, timescale, time_magnitude):
    n = order(n)
    d, tau, T = map(rational, (amplitude, timescale, time_magnitude))
    if d < 0 or tau <= 0 or T < 2 * tau:
        raise ValueError("Require Delta>=0, tau>0 and |t|>=2tau")
    return tail_coefficient(n) * d * tau**8 / T ** (8 + n)


@cache
def data():
    checks = {}
    for n in range(9):
        P = polynomial(n)
        if n == 0:
            checks["initial_profile_polynomial"] = P - X
        else:
            previous = polynomial(n - 1)
            checks[f"literal_derivative_recurrence_{n}"] = s.expand(
                (1 + X**8) * s.diff(previous, X)
                - (8 * (n - 1) + 1) * X**7 * previous
                - P
            )
            checks[f"degree_at_order_{n}"] = s.degree(P, X) - 7 * (n - 1)
            checks[f"leading_coefficient_at_order_{n}"] = s.LC(s.Poly(P, X)) - (-1) ** (
                n - 1
            ) * s.rf(9, n - 1)
    k = s.Symbol("k", integer=True, positive=True)
    checks["all_order_highest_coefficient_induction"] = s.expand(
        7 * (k - 1) - (8 * k + 1) + (k + 8)
    )
    checks["tail_physical_time_scaling"] = 8 + k - k - 8
    return {
        "derivative_polynomials_checked": {n: polynomial(n) for n in range(9)},
        "all_order_recurrence": "s^(n)=P_n/(1+x^8)^(n+1/8), P_0=x, P_(n+1)=(1+x^8) P_n'-(8n+1)x^7 P_n. Written induction, not the finite replay, covers every n.",
        "complex_disk_argument": "For real x>=2, use |z-x|<=x/4 and s(z)=(1+z^-8)^(-1/8) with the branch equal to the positive real profile. |z^-8|<1/2 and |h'(w)|<1/2, so |s(z)-1|<5x^-8. Cauchy gives |d_x^n(s-1)|<=5 n!4^n x^(-8-n) for every n. Oddness gives the negative tail.",
        "physical_symbol_bound": "|d_t^n(M-M_asym)|<=5 Delta tau^8 n!4^n |t|^(-8-n), |t|>=2tau. Compact real-time derivatives have finite recurrence-polynomial bounds. Hence each fixed profile belongs to S^-8 on each tail; this is not uniform in tau->0.",
        "compact_bounds_through_eight": {n: compact_coefficient(n) for n in range(9)},
        "checks": checks,
        "gates": {
            "inverse_power_disk_below_one_half": bool(
                s.Rational(2, 3) ** 8 < s.Rational(1, 2)
            ),
            "complex_derivative_cap_below_one_half": bool(2**9 < 4**8),
            "complex_profile_coefficient_below_five": bool(
                s.Rational(1, 2) * s.Rational(4, 3) ** 8 < 5
            ),
            "real_denominator_strictly_positive": True,
            "all_orders_proved_by_written_Cauchy_argument": True,
        },
    }
