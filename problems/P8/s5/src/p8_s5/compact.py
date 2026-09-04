"""Background-scale compactification of the full D lapse Hamiltonian.

At each background time separately choose ell=tau*sqrt(1+u^2). This is a
constant local choice of units, NOT a global time-dependent canonical map.
No derivatives of a time-dependent mode normalization have been discarded.
"""

from functools import cache

import sympy as sp

from . import lapse_series
from . import nonlinear_d as m

x = sp.Symbol("x", real=True)
z = 1-x**2
P = 1-6*x**4+10*x**6-4*x**8


def normalized(expr, d_power, parity=0):
    """Map expr(u)*d^d_power to x=u/sqrt(d), keeping parity exactly."""
    d_power = sp.Rational(d_power)
    power = d_power + sp.Rational(parity, 2)
    if not power.is_Integer:
        raise ValueError("Incompatible parity and normalization weight")
    if expr == 0:
        return sp.Integer(0)
    num, den = sp.fraction(sp.cancel(expr/m.u**parity))
    y = sp.Symbol("y")
    def even(expr):
        poly = sp.Poly(expr, m.u)
        if any(i[0] % 2 for i, value in poly.terms() if value):
            raise ValueError("Coefficient has the wrong time parity")
        return sp.Add(*(value*y**(i[0]//2) for i, value in poly.terms()))
    return sp.factor(x**parity*(even(num)/even(den)).subs(y, x**2/z)/z**power)


@cache
def explicit_functions():
    """Local-unit covariant F jets, written independently as compact polynomials."""
    f0 = -4*(1+x**2)
    fx = 2*(2*x**2-1)*(1-6*z**3)
    fxx = -2*z**3*(45*x**8+122*x**6*z+96*x**4*z**2-12*x**2*z**3+5*z**4)
    omega = z**3*(m.N**-4-1)/2
    omega_time = -3*x*z**3*(m.N**-4-1)
    F = f0+fx*(m.X-1)+fxx*(m.X-1)**2/2
    return {"F0": f0, "FX": fx, "FXX": fxx, "F": F,
            "H": 2*x, "omega": omega, "omega_time": omega_time}


@cache
def hamiltonian():
    f = explicit_functions()
    Q = -3*f["H"]**2/2+f["H"]*m.sigma-m.sigma**2/6+m.shear2
    return (2*m.N*sp.exp(-3*f["omega"])*Q
            -2*f["omega_time"]*(-3*f["H"]+m.sigma)
            -m.N*sp.exp(3*f["omega"])*f["F"].subs(m.X, m.N**-2)
            -m.N*sp.exp(f["omega"])*m.rho/2)


@cache
def checks():
    f = m.functions()
    explicit = explicit_functions()
    out = {"H": sp.cancel(normalized(m.H, sp.Rational(1, 2), 1)-explicit["H"]),
           "omega": sp.cancel(normalized(f["omega"], 0)-explicit["omega"]),
           "omega_time": sp.cancel(normalized(f["omega_u"], sp.Rational(1, 2), 1)-explicit["omega_time"]),
           "weighted_J": sp.cancel(normalized(f["J"], 1)-2*P)}
    for k, key in enumerate(("F0", "FX", "FXX")):
        old = sp.diff(f["F"], m.X, k).subs(m.X, 1)
        out[key] = sp.cancel(normalized(old, 1)-explicit[key])
    # Independent compact-Hamiltonian Hessian at fixed normalized momenta.
    h = hamiltonian().subs({m.sigma: 0, m.rho: 0, m.shear2: 0})
    out["lapse_equation"] = sp.factor(sp.diff(h, m.N).subs(m.N, 1))
    out["lapse_hessian"] = sp.factor(sp.diff(h, m.N, 2).subs(m.N, 1)+4*P)
    lam = 1-2*z**3
    out["unitary_principal_kinetic"] = sp.cancel(normalized(f["J"]/f["Theta"]**2, 0)-P/(2*x**2))
    out["gamma_principal_kinetic_local_b"] = sp.cancel(normalized(f["J"]/f["Lambda"]**2, 1)-2*P/lam**2)
    return out


def bound_coefficient(value):
    """Explicit conservative bound on [-1,1]: coefficient L1 norm / (P>=1/4)^k."""
    numerator, denominator = sp.fraction(sp.cancel(value))
    degree = sp.degree(denominator, x)
    if degree % 8:
        raise ValueError("Unexpected compact denominator degree")
    power = int(degree/8)
    constant = sp.LC(sp.Poly(denominator, x))/(-4)**power
    if constant == 0 or sp.cancel(denominator-constant*P**power) != 0:
        raise ValueError("Unexpected compact denominator; a new domain proof is needed")
    l1 = sum(abs(c) for c in sp.Poly(numerator, x).all_coeffs())
    return {"expression": str(value), "absolute_bound": str(l1*4**power/abs(constant)),
            "denominator_P_power": power, "numerator_L1_norm": str(l1),
            "denominator_constant_abs": str(abs(constant)),
            "past_endpoint": str(value.subs(x, -1)), "future_endpoint": str(value.subs(x, 1))}


@cache
def coefficient_report():
    data = {}
    for group, orders in lapse_series.derive().items():
        data[group] = []
        for expr in orders:
            entries = {}
            for powers, coefficient in lapse_series.coefficients(expr).items():
                p, r, w = map(int, powers.split(","))
                weight = (1 if group == "hamiltonian" else 0)-sp.Rational(p, 2)-r-w
                value = normalized(coefficient, weight, p % 2)
                entries[powers] = bound_coefficient(value)
            data[group].append(entries)
    return data
