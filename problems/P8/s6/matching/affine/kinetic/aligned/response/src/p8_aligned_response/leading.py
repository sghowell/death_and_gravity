"""Actual quadratic source and sufficient bounds from light-field time jets."""
from functools import cache

import sympy as sp
from p8_affine import dictionary, lower
from p8_affine_aligned import alignment
from p8_affine_retuned.bounds import exact

u = dictionary.u
N, K = sp.symbols("n delta_K_hat", real=True)


@cache
def source():
    x = dictionary.x
    bg = alignment.parent.old.background()
    # Q=Q_x=0 at x=-1, so differentiating Q_x+A Q=F gives Q_xx=F_x.
    qxx = sp.factor(sp.diff(lower.ode()["forcing"], x).subs(x, -1))
    Q2 = sp.factor(2*qxx*N**2)  # x+1=2*eps*n+O(eps²).
    actual = alignment.rolling()["second_order"].subs(sp.Symbol("Q_second_order", real=True), Q2)
    A = -2/bg["h"]
    C = -6*bg["H"]/bg["h"]-sp.Rational(27, 8)*bg["H"]/bg["h"]**2
    return {"Q_xx": qxx, "Q_second": Q2, "mixed_coefficient": sp.factor(A),
            "squared_lapse_coefficient": sp.factor(C), "actual": actual,
            "expected": A*N*K+C*N**2,
            "ODE_second_jet": sp.factor(qxx+3*sp.diff(bg["h"], u)/(4*bg["h"]**3)),
            "full_quadratic_source": sp.factor(actual-A*N*K-C*N**2)}


def compact_rational_bound(value):
    """Coefficient-sum numerator bound; require the actual positive denominator."""
    numerator, denominator = sp.fraction(sp.factor(value))
    poly = sp.Poly(numerator, u)
    den = sp.Poly(sp.expand(denominator), u)
    if den.TC() <= 0 or any(powers[0] % 2 or coefficient < 0 for powers, coefficient in den.terms()):
        raise ValueError("Denominator lacks the stated continuous positive-even witness")
    return sum(abs(coefficient)*sp.Rational(1, 2)**powers[0] for powers, coefficient in poly.terms())/den.TC()


@cache
def jet_constants():
    data = source()
    A = [compact_rational_bound(sp.diff(data["mixed_coefficient"], u, j)) for j in range(4)]
    C = [compact_rational_bound(sp.diff(data["squared_lapse_coefficient"], u, j)) for j in range(4)]
    mixed = [sum(sp.binomial(j, i)*A[i]*2**(j-i) for i in range(j+1)) for j in range(4)]
    squared = [sum(sp.binomial(j, i)*C[i]*2**(j-i) for i in range(j+1)) for j in range(4)]
    return {"A_coefficient_jets": A, "C_coefficient_jets": C,
            "mixed_source_jets": mixed, "squared_source_jets": squared,
            "mixed_envelope": max(mixed), "squared_envelope": max(squared)}


def source_envelope(lapse_jets, trace_jets):
    en, ek = exact(lapse_jets, "epsilon_n"), exact(trace_jets, "epsilon_K")
    if en.is_nonnegative is not True or ek.is_nonnegative is not True:
        raise ValueError("Require nonnegative bounds for all time jets through order three")
    data = jet_constants()
    return {"epsilon_n": en, "epsilon_K": ek,
            "source_jet_E": data["mixed_envelope"]*en*ek+data["squared_envelope"]*en**2,
            "requires_lapse_and_trace_jets_through_order_three": True,
            "sufficient_preparation": "n=n'=n''=0 initially; arbitrary finite trace jets",
            "bound_on_full_nonlinear_source_or_on_shell_evolution": False}


@cache
def checks():
    data = source()
    nf, kf = sp.Function("n")(u), sp.Function("K")(u)
    expr = data["expected"].subs({N: nf, K: kf})
    prepared = {nf: 0, sp.diff(nf, u): 0, sp.diff(nf, u, 2): 0}
    return {"original_ODE_second_jet": data["ODE_second_jet"],
            "actual_quadratic_source_not_arbitrary_Q2": data["full_quadratic_source"],
            "prepared_source_value": expr.xreplace(prepared),
            "prepared_source_first_jet": sp.diff(expr, u).xreplace(prepared),
            "prepared_source_second_jet": sp.diff(expr, u, 2).xreplace(prepared)}
