"""Exact dimensionally regulated massive/massive/massless Gaussian masters."""

from functools import cache

import sympy as s

from .tensors import A, B, C, d, e, x, y, z


def shifted_gamma_ratio(v, n):
    if n >= 0:
        return s.prod((v + j for j in range(n)), start=s.Integer(1))
    return s.Integer(1) / s.prod((v + j for j in range(n, 0)), start=s.Integer(1))


@cache
def master(a, b, c):
    if a <= 0 or b <= 0:
        return s.Integer(0)
    return s.factor(
        shifted_gamma_ratio(2 - e, -c)
        * shifted_gamma_ratio(2 * e, a + b + c - 4)
        / shifted_gamma_ratio(2 * e, a + b + 2 * c - 4)
        * shifted_gamma_ratio(e, a + c - 2)
        * shifted_gamma_ratio(e, b + c - 2)
        / (s.factorial(a - 1) * s.factorial(b - 1))
    )


def integrate(expr):
    expr = s.expand(expr.subs({x: A - 1, y: B - 1, z: (A + B - C - 2) / 2}) / C)
    rows = {}
    for term in s.Add.make_args(expr):
        pd = term.as_powers_dict()
        powers = tuple(-int(pd.get(v, 0)) for v in (A, B, C))
        coeff = term.subs({A: 1, B: 1, C: 1})
        rows[powers] = rows.get(powers, 0) + coeff
    value = s.factor(
        sum(
            coeff.subs(d, 4 - 2 * e) * master(*powers) for powers, coeff in rows.items()
        )
    )
    return value, rows
