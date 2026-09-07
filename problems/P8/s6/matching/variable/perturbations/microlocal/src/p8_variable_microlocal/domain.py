"""Exact continuous positivity of U-normalized candidate K/G matrices.

The rational principal coefficients extend to c=2 after the field weight;
this is only a coefficient extension, not a background/action there. Time
connections and exact Cauchy inverse thresholds are not uniform in that
singular limit.
"""

from functools import cache
from itertools import product

import sympy as sp

from . import leading


def bernstein(polynomial, variables):
    p = sp.Poly(polynomial, *variables)
    degrees = tuple(p.degree(x) for x in variables)
    values = []
    for index in product(*(range(n+1) for n in degrees)):
        value = sp.S.Zero
        for power, coefficient in p.terms():
            if all(i <= j for i, j in zip(power, index, strict=True)):
                value += coefficient*sp.prod(
                    sp.binomial(j, i)/sp.binomial(n, i)
                    for i, j, n in zip(power, index, degrees, strict=True))
        values.append(value)
    return degrees, tuple(values)


@cache
def derive():
    d = leading.derive()
    u, c, dd = d["u"], d["c"], d["d"]
    x, z = sp.symbols("x z", real=True)
    records = {}
    specifications = {
        "kinetic_trace_upper": (64-sp.trace(d["normalized_kinetic"][:2, :2]),
                                 100*c*(1-u**2)**2*dd**8),
        "gradient_trace_upper": (64-sp.trace(d["normalized_gradient"][:2, :2]),
                                  400*c*(1-u**2)**3*dd**12),
        "kinetic_minor_lower": (d["normalized_kinetic"][:2, :2].det()-16,
                                 200*c**2*(1-u**2)**2*dd**2),
        "gradient_minor_lower": (d["normalized_gradient"][:2, :2].det()-16,
                                  1600*c**2*(1-u**2)**3*dd**10),
    }
    for name, (margin, denominator) in specifications.items():
        numerator = sp.cancel(margin*denominator)
        polynomial = sp.Poly(numerator, u, c)
        if any(power[0] % 2 for power in polynomial.monoms()):
            raise ValueError("The principal margin was not an even time polynomial")
        compact = sum(coefficient*(x/100)**(power[0]//2)*(2+2*z)**power[1]
                      for power, coefficient in polynomial.terms())
        degrees, coefficients = bernstein(compact, (x, z))
        records[name] = {"denominator": denominator, "numerator": numerator,
                         "margin_identity": sp.cancel(margin-numerator/denominator),
                         "degrees": degrees, "coefficients": coefficients}
    return {"records": records,
            "clock_kinetic_identity": sp.factor(d["normalized_kinetic"][1, 1]-dd**6*d["kbar"]),
            "clock_gradient_identity": sp.factor(d["normalized_gradient"][1, 1]-dd**2*d["kbar"]),
            "claimed_matrix_bounds": "(1/4) I < normalized K_s, normalized G_s < 64 I",
            "domain": "u in [-1/10,1/10], c in [2,4] only for the extended rational coefficients"}


def checks():
    d = derive()
    if d["clock_kinetic_identity"] or d["clock_gradient_identity"]:
        raise ValueError("The canonical-clock diagonal bridge failed")
    for record in d["records"].values():
        if record["margin_identity"] != 0 or any(value <= 0 for value in record["coefficients"]):
            raise ValueError("The continuous principal coefficient positivity proof failed")
    return {name: {"count": len(record["coefficients"]), "minimum": min(record["coefficients"])}
            for name, record in d["records"].items()}
