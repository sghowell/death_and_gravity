"""Exact volume-normalized frequencies and conservative all-time WKB bounds."""

from functools import cache

import sympy as sp
from p8_physical.vertices import chart_functions
from p8_s5 import compact

x, P = compact.x, compact.P
q = sp.Symbol("q", positive=True)
r = sp.Symbol("r", nonnegative=True)
lam = 1-2*(1-x*x)**3
R = lam**2-2*P*r


def derivative(value, weight, reciprocal=False):
    """ell^(2w+1) d(value/ell^(2w))/dt in a fixed local evaluation patch."""
    momentum_part = 2*x*r*sp.diff(value, r) if reciprocal else -2*x*q*sp.diff(value, q)
    return sp.cancel((1-x*x)*sp.diff(value, x)+momentum_part-2*weight*x*value)


@cache
def derive(chart):
    if chart == "tensor":
        K, B, C, f = sp.Rational(1, 2), sp.Integer(0), q, 3*x
    elif chart in ("unitary", "gamma"):
        functions = chart_functions(chart, q)
        K = functions["K"]
        B = -functions["beta"]/functions["alpha"]
        J, theta = 2*P, 2*x
        C = 2*(lam**2*q*q/J-q) if chart == "unitary" else 2*theta**2*q*q/J
        f = 3*x+functions["logZdot"]
    else:
        raise ValueError("Use unitary, gamma or tensor")
    b = sp.cancel(f+B)
    omega2 = sp.factor(C/(2*K)-b*b-derivative(b, sp.Rational(1, 2)))
    mass = sp.factor(omega2-q)
    return {"K": K, "B": sp.factor(B), "C": sp.factor(C), "f": sp.factor(f),
            "omega_squared": omega2, "mass_correction": mass}


def coefficient_bound(value, chart):
    """Exact polynomial L1 majorant on a proved covering domain, r<=1/68."""
    value = sp.cancel(value)
    numerator, denominator = sp.fraction(value)
    constant, factors = sp.factor_list(denominator, x, r)
    denominator_lower = abs(constant)
    encoded = []
    candidates = [(P, sp.Rational(1, 4), "P")]
    if chart == "unitary":
        # |x|>=1/sqrt(65)>1/9; this rational weakening avoids radicals.
        candidates.append((x, sp.Rational(1, 9), "abs(x)"))
    elif chart == "gamma":
        candidates.append((R, sp.Rational(15, 68), "Lambda^2-2*P*r"))
    elif chart != "tensor":
        raise ValueError("Unknown chart")
    for factor, power in factors:
        for base, lower, name in candidates:
            ratio = sp.cancel(factor/base)
            if not ratio.has(x, r) and ratio != 0:
                denominator_lower *= (abs(ratio)*lower)**power
                encoded.append({"factor": name, "power": int(power), "absolute_lower": str(lower)})
                break
        else:
            raise ValueError(f"Unproved oscillator denominator: {factor}")
    x_upper = sp.Rational(1, 4) if chart == "gamma" else sp.Integer(1)
    l1 = sum(abs(coefficient)*x_upper**powers[0]*sp.Rational(1, 68)**powers[1]
             for powers, coefficient in sp.Poly(numerator, x, r).terms())
    return {"expression": str(sp.factor(value)), "absolute_bound": str(l1/denominator_lower),
            "weighted_numerator_L1": str(l1), "denominator_absolute_lower": str(denominator_lower),
            "denominator_factors": encoded}


@cache
def bounds():
    data = {}
    all_bounds = [sp.Integer(10)**8]
    for chart in ("unitary", "gamma", "tensor"):
        mass = sp.cancel(derive(chart)["mass_correction"].subs(q, 1/r))
        first = derivative(mass, 1, reciprocal=True)
        second = derivative(first, sp.Rational(3, 2), reciprocal=True)
        entries = [coefficient_bound(value, chart) for value in (mass, first, second)]
        C0, C1, C2 = [sp.Rational(entry["absolute_bound"]) for entry in entries]
        data[chart] = {"mass_and_derivative_bounds": entries}
        all_bounds.extend([2*C0, C1, C2])
    threshold = sp.Integer(10)**8
    while threshold < max(all_bounds):
        threshold *= 10
    return {"charts": data, "q_threshold": str(threshold),
            "omega_squared_bounds": ["q/2", "3*q/2"],
            "first_adiabatic_ratio_bound": f"5*sqrt(2)/sqrt(q) <= 5*sqrt(2)/sqrt({threshold})",
            "second_adiabatic_ratio_bound": f"92/q <= 92/{threshold}",
            "warning": "Sufficient conservative linear-mode threshold, NOT a strong-coupling cutoff"}


@cache
def generic_checks():
    # Independent quadratic completion and integration-by-parts identity.
    K, B, C, f, bd, Qc, V = sp.symbols("K B C f bd Qc V")
    A = 1/(2*K)
    Q, Qdot = Qc/sp.sqrt(2*K), (V-f*Qc)/sp.sqrt(2*K)
    p0 = (Qdot-B*Q)/A
    L = sp.expand(p0*Qdot-(A*p0*p0/2+B*p0*Q+C*Q*Q/2))
    target = (V-(f+B)*Qc)**2/2-A*C*Qc*Qc/2
    completed = sp.simplify(L-target)
    after_ibp = sp.expand(target+(f+B)*Qc*V+bd*Qc*Qc/2)
    oscillator = (V*V-(A*C-(f+B)**2-bd)*Qc*Qc)/2
    tensor_mass = derive("tensor")["mass_correction"]
    unitary = derive("unitary")
    return {"quadratic_completion": completed, "oscillator_IBP": sp.expand(after_ibp-oscillator),
            "tensor_mass": sp.cancel(tensor_mass+3+3*x*x),
            "unitary_principal_bridge": sp.cancel(unitary["mass_correction"]+unitary["f"]**2
                                                   +derivative(unitary["f"], sp.Rational(1, 2)))}
