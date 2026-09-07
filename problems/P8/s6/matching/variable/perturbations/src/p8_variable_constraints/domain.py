"""Exact secondary-constraint denominator bound on the frozen clock box.

Bernstein coefficients are certificates for a polynomial on a continuous
rectangle, not samples. The c=2 edge is used only in a cleared polynomial;
no finite-coefficient parent solution is asserted at c=2.
"""

from functools import cache
from itertools import product

import sympy as sp

from . import action, reduction


def bernstein_coefficients(polynomial, variables):
    p = sp.Poly(polynomial, *variables)
    degrees = tuple(p.degree(v) for v in variables)
    coefficients = []
    for index in product(*(range(degree+1) for degree in degrees)):
        total = sp.S.Zero
        for power, coefficient in p.terms():
            if all(power[i] <= index[i] for i in range(len(variables))):
                total += coefficient*sp.prod(sp.binomial(index[i], power[i])/sp.binomial(degrees[i], power[i])
                                              for i in range(len(variables)))
        coefficients.append(total)
    return degrees, tuple(coefficients)


@cache
def derive():
    d = action.derive()
    u, c, dd = d["u"], d["c"], d["d"]
    F = 6400*(1-u**2)-800*c*(1-u**2)*dd**12-c*dd**2
    denominator = 24*c*(1-u**2)**2*dd**6*F
    margin_numerator = sp.cancel((-reduction.derive()["expected_secondary_D"]-sp.Rational(1, 8))*denominator)
    p = sp.Poly(margin_numerator, u, c)
    x, z = sp.symbols("x z", real=True)
    if any(power[0] % 2 for power in p.monoms()):
        raise ValueError("The secondary margin was not even in time")
    compact = sum(coefficient*(x/100)**(power[0]//2)*c**power[1] for power, coefficient in p.terms())
    degree_one, coefficients_one = bernstein_coefficients(compact.subs(c, 1), (x,))
    degree_positive, coefficients_positive = bernstein_coefficients(compact.subs(c, 2+2*z), (x, z))
    return {"F": F, "denominator": denominator, "numerator": margin_numerator,
            "clock_F_identity": sp.factor(F-100*c*dd**14*d["kbar"]),
            "margin_identity": sp.factor(margin_numerator/denominator
                                           +reduction.derive()["expected_secondary_D"]+sp.Rational(1, 8)),
            "c1_degree": degree_one, "c1_coefficients": coefficients_one,
            "positive_branch_degree": degree_positive, "positive_branch_coefficients": coefficients_positive}


def checks():
    d = derive()
    if d["clock_F_identity"] != 0 or d["margin_identity"] != 0:
        raise ValueError("The literal constraint/clock denominator bridge failed")
    groups = {"c1": d["c1_coefficients"], "2_to_4_cleared_polynomial": d["positive_branch_coefficients"]}
    if any(value.is_positive is not True for values in groups.values() for value in values):
        raise ValueError("The exact continuous secondary constraint margin was not strict")
    return {name: {"count": len(values), "minimum": min(values)} for name, values in groups.items()}
