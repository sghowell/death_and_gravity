"""Pole-cancelled, delta-dependent loading by an actual conserved g TT source.

All source bounds are linear-probe bounds. The fixed support is not a low
temporal-frequency band, and the finite cutoff cost is intentionally coarse.
"""

from functools import cache
from math import comb, factorial

import sympy as sp

from .exact import DELTA_RADIUS, nonnegative, parameters


def loading_coefficients(coordinate, kg, kf, gg, inverse_spring, g_target, f_target, momentum):
    """Return coefficients of zeta^(1..4) in a^3*sigma/4.

    The target solves the full homogeneous equations. P=inverse_spring has
    a zero, not a pole, at (u,delta)=(0,0). No derivative of an undefined
    delta=0 center spring is used by this representation.
    """
    u = coordinate
    m = 2*kf*sp.diff(f_target, u)+sp.diff(kf, u)*f_target
    n = kf*f_target
    aa, bb = inverse_spring*m, inverse_spring*n
    coefficients = {
        1: 2*kg*sp.diff(g_target, u)+sp.diff(kg, u)*g_target
           +sp.diff(kg*sp.diff(aa, u), u)+gg*momentum*aa+m,
        2: kg*g_target+kg*sp.diff(aa, u)+sp.diff(kg*(aa+sp.diff(bb, u)), u)+gg*momentum*bb+n,
        3: kg*(aa+sp.diff(bb, u))+sp.diff(kg*bb, u),
        4: kg*bb,
    }
    return {"M": m, "N": n, "A": aa, "B": bb, "coefficients": coefficients}


def _jets(norm, order=3):
    # Shrink the coefficient-l1 u radius from 1/20 to 1/40. Delta radius
    # remains exactly v, so subsequent delta differences retain delta/v.
    return [sp.Rational(factorial(j)*40**j)*norm for j in range(order+1)]


def _product(left, right):
    return [sum(comb(n, j)*left[j]*right[n-j] for j in range(n+1))
            for n in range(min(len(left), len(right)))]


@cache
def coefficient_bounds():
    """Coefficient-l1 sigma bounds per |even amplitude|+|odd amplitude|."""
    kg = kf = gg = _jets(sp.Integer(1))
    p = _jets(sp.Rational(1, 100))
    g = f = _jets(sp.Integer(5))
    # M=2 Kf F'+Kf'F; N=Kf F, through two/three derivatives.
    m = [2*sum(comb(n, j)*kf[j]*f[n-j+1] for j in range(n+1))
         +sum(comb(n, j)*kf[j+1]*f[n-j] for j in range(n+1)) for n in range(3)]
    n = _product(kf, f)
    aa, bb = _product(p, m), _product(p, n)
    s1 = (2*kg[0]*g[1]+kg[1]*g[0]+kg[1]*aa[1]+kg[0]*aa[2]
          +4*gg[0]*aa[0]+m[0])
    s2 = (kg[0]*g[0]+kg[0]*aa[1]+kg[1]*(aa[0]+bb[1])
          +kg[0]*(aa[1]+bb[2])+4*gg[0]*bb[0]+n[0])
    s3 = kg[0]*(aa[0]+bb[1])+kg[1]*bb[0]+kg[0]*bb[1]
    s4 = kg[0]*bb[0]
    # ||4/a^3|| <=4/(1-v)^6 <5.
    return {1: 5*s1, 2: 5*s2, 3: 5*s3, 4: 5*s4}


@cache
def switch_bounds():
    """Explicit flat switch from e^(-1/x), with denominator >=e^-2>1/9.

    zeta(u)=s(200*(u+1/50)-1/2). The transition is contained strictly in
    (-1/50,-1/100). Its derivatives are NOT claimed small.
    """
    y = sp.Symbol("y")
    polynomial = sp.Integer(1)
    bump = []
    for _ in range(5):
        poly = sp.Poly(polynomial, y)
        bump.append(sum(abs(coefficient)*factorial(power[0]) for power, coefficient in poly.terms()))
        polynomial = sp.expand(y**2*(polynomial-sp.diff(polynomial, y)))
    inverse = [sp.Integer(9)]
    for order in range(1, 5):
        inverse.append(9*sum(comb(order, j)*2*bump[j]*inverse[order-j] for j in range(1, order+1)))
    switch = {0: sp.Integer(1)}
    for order in range(1, 5):
        switch[order] = 200**order*sum(comb(order, j)*bump[j]*inverse[order-j] for j in range(order+1))
    return {"bump_jets": bump, "reciprocal_jets": inverse, "zeta_jets": switch,
            "denominator_lower": sp.Rational(1, 9), "support_length": sp.Rational(1, 100)}


def source_sup_bound(amplitude_l1=1, cutoff_jets=None):
    amplitude_l1 = nonnegative(amplitude_l1, "amplitude l1")
    cutoff_jets = switch_bounds()["zeta_jets"] if cutoff_jets is None else cutoff_jets
    if (not isinstance(cutoff_jets, dict) or any(index not in cutoff_jets for index in range(1, 5))
            or any(isinstance(key, bool) or not isinstance(key, int) or key not in range(5)
                   for key in cutoff_jets)):
        raise ValueError("All four fixed cutoff derivative suprema are required")
    if 0 in cutoff_jets:
        nonnegative(cutoff_jets[0], "cutoff jet 0")
    return amplitude_l1*sum(coefficient_bounds()[index]*nonnegative(cutoff_jets[index], f"cutoff jet {index}")
                            for index in range(1, 5))


def source_difference_bound(delta, amplitude_l1=1, cutoff_jets=None, momentum_squared=1):
    delta, _ = parameters(delta, momentum_squared)
    return delta/DELTA_RADIUS*source_sup_bound(amplitude_l1, cutoff_jets)


@cache
def checks():
    u = sp.Symbol("u", real=True)
    kg, kf, gg, p, g, f, zeta = (sp.Function(name)(u) for name in ("kg", "kf", "gg", "P", "g", "f", "zeta"))
    momentum = sp.Symbol("K", real=True)
    d = loading_coefficients(u, kg, kf, gg, p, g, f, momentum)
    correction = p*(sp.diff(zeta, u)*d["M"]+sp.diff(zeta, u, 2)*d["N"])
    direct = (2*kg*sp.diff(zeta, u)*sp.diff(g, u)+sp.diff(kg*sp.diff(zeta, u), u)*g
              +sp.diff(kg*sp.diff(correction, u), u)+gg*momentum*correction
              +sp.diff(zeta, u)*d["M"]+sp.diff(zeta, u, 2)*d["N"])
    assembled = sum(value*sp.diff(zeta, u, order) for order, value in d["coefficients"].items())
    return {"full_pole_cancelled_commutator_loading": sp.expand(direct-assembled),
            "fixed_cutoff_transition_left": (sp.Rational(1, 2)/200-sp.Rational(1, 50))+sp.Rational(7, 400),
            "fixed_cutoff_transition_right": (sp.Rational(3, 2)/200-sp.Rational(1, 50))+sp.Rational(1, 80),
            "source_prefactor_margin": 5-4/(1-DELTA_RADIUS)**6}
