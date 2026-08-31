"""Shared symbols and exact-series helpers for the P9(b) S1 chain.

Conventions (FORMULATION 1.1): c = 1, signature (-,+,+,+), flat FLRW
ds^2 = -dt^2 + a^2 dx^i dx^i, H = adot/a, overdot = d/dt.  Scalar sector in
unitary gauge, GLV14 eq. (73): N = 1 + dN, N^i = delta^{ij} d_j psi,
h_ij = a^2 e^{2 zeta} delta_ij.  Single Fourier mode: every perturbation is
eps * f(t) * cos(k x); products are averaged over x (exact, quadratic order).
All expansions are exact Taylor polynomials in eps up to eps^2 ("jets").
"""
from __future__ import annotations

import sympy as sp

# --- coordinates, bookkeeping ------------------------------------------------------
t, x, y1, y2 = sp.symbols("t x y1 y2", real=True)
COORDS3 = (x, y1, y2)
COORDS4 = (t, x, y1, y2)
k = sp.Symbol("k", positive=True)
eps = sp.Symbol("epsilon", real=True)
cx, sx = sp.cos(k * x), sp.sin(k * x)

# --- background functions ----------------------------------------------------------
a = sp.Function("a", positive=True)(t)
H = sp.Function("H", positive=True)(t)
M2 = sp.Function("M2", positive=True)(t)        # effective M^2 (= BS14 M_*^2)
aK, aBg, aT, aHh = (sp.Function(n)(t) for n in ("alphaK", "alphaBglv", "alphaT", "alphaH"))
Ms2 = sp.Symbol("Mstar2", positive=True)        # GPV bare M_*^2 (constant)
MPl2 = sp.Symbol("MPl2", positive=True)         # frozen M_Pl^2 (family B)
f_, cc, Lam, M24, m33, m42, mt42 = (sp.Function(n)(t) for n in
                                    ("f", "c", "Lambda", "M24", "m33", "m42", "mt42"))
# matter k-essence P(Y): background Taylor coefficients P0=Pbar, P1=P_Y, ...
sgb = sp.Function("sigma_b")(t)                 # background matter field
P0, P1, P2, P3, P4 = (sp.Function(n)(t) for n in ("P0", "P1", "P2", "P3", "P4"))
# KGB functions K(X), G(X): background Taylor coefficients
phb = sp.Function("phi_b")(t)
K0, K1, K2, K3 = (sp.Function(n)(t) for n in ("K0", "K1", "K2", "K3"))
G0, G1, G2, G3 = (sp.Function(n)(t) for n in ("G0", "G1", "G2", "G3"))
# fluid translation symbols
rhom, pm, cm2 = sp.symbols("rho_m p_m c_m2", real=True)

# --- perturbation amplitude functions ----------------------------------------------
zf, dnf, psf, dsf = (sp.Function(n)(t) for n in ("zeta", "dN", "psi", "dsigma"))
PERT = (zf, dnf, psf, dsf)


def cut(e):
    """Drop O(eps^3) from a polynomial-in-eps expression (exact)."""
    e = sp.expand(e)
    if not e.has(eps):
        return e
    return sum(e.coeff(eps, n) * eps**n for n in range(3))


def mul(*fs):
    """Product of eps-jets, truncated at eps^2."""
    out = sp.Integer(1)
    for f in fs:
        out = cut(out * f)
    return out


def jetpow(base0, dbase, p):
    """(base0 + dbase)^p as an eps-jet; base0 eps-free, dbase = O(eps)."""
    r = cut(dbase / base0)
    return cut(base0**p * (1 + p * r + sp.Rational(p * (p - 1), 2) * r**2)) \
        if isinstance(p, int) else \
        cut(base0**p * (1 + p * r + p * (p - 1) / 2 * r**2))


# --- background derivative rules ---------------------------------------------------
_da1 = a * H
_da2 = sp.diff(_da1, t).subs(sp.Derivative(a, t), _da1)
_da3 = sp.diff(_da2, t).subs(sp.Derivative(a, t), _da1)
_da4 = sp.diff(_da3, t).subs(sp.Derivative(a, t), _da1)
_ARULES = {sp.Derivative(a, (t, 4)): _da4, sp.Derivative(a, (t, 3)): _da3,
           sp.Derivative(a, (t, 2)): _da2, sp.Derivative(a, t): _da1}

# chain rules: P_i(t) are P(Y) Taylor coefficients evaluated on Yb(t) = sgb'^2/2,
# K_i, G_i on Xb(t) = phb'^2/2.  d/dt P_i = P_{i+1} * Yb', etc.
Ybdot = sgb.diff(t) * sgb.diff(t, 2)
Xbdot = phb.diff(t) * phb.diff(t, 2)
_CHAIN = {sp.Derivative(P0, t): P1 * Ybdot, sp.Derivative(P1, t): P2 * Ybdot,
          sp.Derivative(P2, t): P3 * Ybdot, sp.Derivative(P3, t): P4 * Ybdot,
          sp.Derivative(K0, t): K1 * Xbdot, sp.Derivative(K1, t): K2 * Xbdot,
          sp.Derivative(K2, t): K3 * Xbdot,
          sp.Derivative(G0, t): G1 * Xbdot, sp.Derivative(G1, t): G2 * Xbdot,
          sp.Derivative(G2, t): G3 * Xbdot}


def bgsubs(e):
    """Apply a-dot and Taylor-coefficient chain rules until stable."""
    for _ in range(8):
        e2 = e.xreplace(_ARULES).xreplace(_CHAIN)
        if e2 == e:
            return e2
        e = e2
    raise RuntimeError("bgsubs did not stabilise")


def ddt(e):
    """d/dt with background rules applied."""
    return bgsubs(sp.diff(e, t))


def iszero(e):
    """Exact zero test for the rational-in-functions expressions used here."""
    return sp.cancel(sp.together(sp.expand(bgsubs(e)))) == 0


# --- spatial average over the mode -------------------------------------------------
_CS, _SS = sp.symbols("_CS _SS")
_AVG = {(0, 0): sp.Integer(1), (2, 0): sp.Rational(1, 2), (0, 2): sp.Rational(1, 2),
        (1, 1): sp.Integer(0), (1, 0): sp.Integer(0), (0, 1): sp.Integer(0)}


def xavg(e):
    """<e>_x for e polynomial of degree <= 2 in cos(kx), sin(kx) (exact)."""
    e = sp.expand(e).subs({cx: _CS, sx: _SS})
    if not (e.has(_CS) or e.has(_SS)):
        return e
    p = sp.Poly(e, _CS, _SS)
    out = sp.Integer(0)
    for (m, n), coef in zip(p.monoms(), p.coeffs()):
        if (m, n) not in _AVG:
            raise ValueError(f"trig degree too high: {(m, n)}")
        out += coef * _AVG[(m, n)]
    return sp.expand(out)


def subs_fun(e, m):
    """Substitute applied functions {F(t): expr} including their t-derivatives."""
    for fn, val in m.items():
        for order in (3, 2, 1):
            e = e.subs(sp.Derivative(fn, (t, order)), bgsubs(sp.diff(val, (t, order))))
        e = e.subs(fn, val)
    return e
