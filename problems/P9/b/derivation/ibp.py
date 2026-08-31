"""Integration-by-parts canonicaliser for quadratic Lagrangians in time.

A (x-averaged) quadratic Lagrangian L(t) in perturbation amplitudes v_i(t)
is brought to the canonical form

    L_c = sum T_ij vd_i vd_j  +  sum U_ij v_i v_j
          + sum_{i<j} W_ij (vd_i v_j - v_i vd_j),

with L = L_c + dF/dt and the certificate F returned explicitly, so the caller
can verify the identity L - L_c - dF/dt == 0 exactly (house style: every step
an identity test).  Rules used (each an exact IBP identity):

    q v''  w      = d/dt(q v' w) - q' v' w - q v' w'      (w any degree-1 factor)
    q v' v        = d/dt(q v^2 / 2) - (q'/2) v^2
    q (v' w + v w')= d/dt(q v w) - q' v w                 (symmetric part, i != j)
"""
from __future__ import annotations

import sympy as sp

from . import tools as T


def _deriv(v, n=1):
    return sp.Derivative(v, (T.t, n)) if n > 1 else sp.Derivative(v, T.t)


def canon(L, vars_):
    """Return (L_c, F) with L = L_c + dF/dt (exact, caller-verifiable)."""
    L = sp.expand(L)
    F = sp.Integer(0)
    # --- step 1: remove second derivatives -----------------------------------------
    for _ in range(6):
        found = False
        for v in vars_:
            vdd = _deriv(v, 2)
            c = sp.expand(L).coeff(vdd)
            if c == 0:
                continue
            found = True
            # q * vdd, q of degree 1 in fields: q*vdd = d/dt(q v') - q' v'
            L = sp.expand(L - c * vdd + (-T.ddt(c)) * _deriv(v))
            F += c * _deriv(v)
        if not found:
            break
    else:
        raise RuntimeError("second-derivative reduction did not terminate")
    # --- step 2: diagonal v' v terms -----------------------------------------------
    for v in vars_:
        vd = _deriv(v)
        c = sp.expand(L).coeff(vd).coeff(v)  # coefficient of vd*v
        if c != 0:
            L = sp.expand(L - c * vd * v - T.ddt(c) / 2 * v**2)
            F += c * v**2 / 2
    # --- step 3: symmetric off-diagonal v_i' v_j + v_i v_j' --------------------------
    n = len(vars_)
    for i in range(n):
        for j in range(i + 1, n):
            vi, vj = vars_[i], vars_[j]
            A = sp.expand(L).coeff(_deriv(vi)).coeff(vj)
            B = sp.expand(L).coeff(vi).coeff(_deriv(vj))
            s = sp.together((A + B) / 2)
            if s != 0:
                L = sp.expand(L - s * (_deriv(vi) * vj + vi * _deriv(vj))
                              - T.ddt(s) * vi * vj)
                F += s * vi * vj
    return sp.expand(L), F


def check_canon(L, Lc, F):
    """Exact certificate: L - Lc - dF/dt == 0."""
    return T.iszero(L - Lc - T.ddt(F))


def quad_matrices(Lc, vars_):
    """Extract (Tkin, U, W) from a canonical form; entries exact."""
    n = len(vars_)
    Tk = sp.zeros(n, n)
    U = sp.zeros(n, n)
    W = sp.zeros(n, n)
    Le = sp.expand(Lc)
    for i in range(n):
        for j in range(n):
            vi, vj = vars_[i], vars_[j]
            if i == j:
                Tk[i, i] = Le.coeff(_deriv(vi), 2)
                U[i, i] = Le.coeff(vi, 2).subs(
                    [(_deriv(w, 2), 0) for w in vars_] + [(_deriv(w), 0) for w in vars_])
            elif i < j:
                cij = Le.coeff(_deriv(vi)).coeff(_deriv(vj))
                Tk[i, j] = Tk[j, i] = cij / 2
                uij = Le.coeff(vi, 1).coeff(vj, 1).subs([(_deriv(w), 0) for w in vars_])
                U[i, j] = U[j, i] = uij / 2
                W[i, j] = Le.coeff(_deriv(vi)).coeff(vj)   # antisym part survives canon
                W[j, i] = -W[i, j]
    return Tk, U, W


def kleading(e, kvar):
    """Leading large-k value of a rational expression: lim_{k->oo} e (finite)."""
    e = sp.cancel(sp.together(e))
    num, den = sp.fraction(e)
    pn, pd = sp.Poly(num, kvar), sp.Poly(den, kvar)
    dn, dd = pn.degree(), pd.degree()
    if dn < dd:
        return sp.Integer(0)
    if dn == dd:
        return sp.cancel(pn.LC() / pd.LC())
    raise ValueError("expression diverges as k -> oo")
