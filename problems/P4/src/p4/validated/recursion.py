"""Certified order-by-order solution of the coefficient recursion.

For a ``PolySystem`` with residual F(u) = P(t,u) theta u - Q(t,u) and *level
equations*  E_n(u) = sum_(c, r, s) c F^{(r)}_{n+s}  (one list of
(coefficient, residual row, shift) triples per unknown), E_n is affine in u_n
with the exact matrix
    M_n = n D + E,
    D[row, l] = sum c P^{(r)}_{s+sigma}[., l],   E[row, l] = sum c (Psi^{(r)}_s - DQ^{(r)}_s)[., l],
P_k, Psi_k, DQ_k Taylor coefficients at u (they involve u_0..u_{s+sigma} only),
sigma = 1 for theta = d/dt and 0 for theta = t d/dt.  (Coefficient n of theta u is
s_n u_{n+sigma}, s_n = n + sigma; u_n enters F_{n+s} through s_{n-sigma} P_{s+sigma} u_n,
through DP(u_0) u_n paired with theta u, and through DQ u_n.)

Coefficients may be ``arb`` (dm = 1) or delta-polynomials (``Series`` with
cap = dm) carrying a parameter dependence; then M_n = sum_j M_{n,j} delta^j and
the level solve is triangular in the delta-degree.  ``arb_mat.solve`` succeeds
only if M_{n,0} is *proven* invertible: that is the certificate that the formal
series is uniquely determined at every order <= K, and the balls enclose the
exact coefficients.
"""
from __future__ import annotations

from dataclasses import dataclass

from flint import arb, arb_mat

from .arbseries import Series
from .biseries import BiSeries


@dataclass
class LevelEquations:
    rows: list                              # per unknown: list of (coef, r, shift)

    @property
    def max_shift(self):
        return max(s for row in self.rows for (_, _, s) in row)

    def values(self, F, n):
        """E_n as a list (arb or delta-Series) from residual series F."""
        out = []
        for row in self.rows:
            acc = None
            for c, r, s in row:
                if n + s >= 0:
                    term = F[r][n + s] * c
                    acc = term if acc is None else acc + term
            out.append(acc if acc is not None else arb(0))
        return out


def make_series(blocks, cap, dm):
    if dm == 1:
        return Series([b if not isinstance(b, Series) else b[0] for b in blocks], cap)
    return BiSeries.from_blocks(blocks, cap, dm)


def series_from_coefs(coefs, d, cap=None, dm=1, extra_zero=True):
    """d component series from ``coefs`` = list over n of d-lists (u_n = 0 appended)."""
    out = []
    for i in range(d):
        c = [coefs[n][i] for n in range(len(coefs))]
        if extra_zero:
            c = c + [arb(0)]
        out.append(make_series(c, cap, dm))
    return out


def var_series(cap, dm):
    return Series.var(cap) if dm == 1 else BiSeries.var(cap, dm)


def dmat_list(M, d, dm):
    """d x d entries (arb or delta-Series) -> list over delta-degree of arb_mat."""
    def entry(e, j):
        if isinstance(e, Series):
            return e[j]
        return e if j == 0 else arb(0)
    return [arb_mat([[entry(M[i][l], j) for l in range(d)] for i in range(d)]) for j in range(dm)]


def structure_matrices(sys, eqs, coefs, dm=1):
    """(D, E) as lists of arb_mat over the delta-degree, from the low orders only."""
    d, sigma = sys.d, sys.sigma
    kmax = eqs.max_shift + sigma
    cap = kmax + 2
    u = series_from_coefs(coefs[:kmax + 1], d, cap=cap, dm=dm, extra_zero=False)
    t = var_series(cap, dm)
    Pser = sys.P_series(u, cap=cap, t=t)
    Psi = sys.psi_series(u, cap=cap, t=t)
    dQ = sys.dQ_series(u, cap=cap, t=t)
    D = [[arb(0)] * d for _ in range(d)]
    E = [[arb(0)] * d for _ in range(d)]
    for i, row in enumerate(eqs.rows):
        for c, r, s in row:
            for l in range(d):
                D[i][l] = D[i][l] + Pser[r][l][s + sigma] * c
                if s >= 0:
                    E[i][l] = E[i][l] + (Psi[r][l][s] - dQ[r][l][s]) * c
    return dmat_list(D, d, dm), dmat_list(E, d, dm)


def level_matrix(D, E, n):
    return [Dj * n + Ej for Dj, Ej in zip(D, E)]


def solve_level(M, rhs, dm):
    """Solve  (sum_j M_j delta^j) u = -rhs  for u (delta-triangular); rhs entries arb or Series."""
    d = M[0].nrows()

    def rhs_vec(k):
        return arb_mat([[-(x[k] if isinstance(x, Series) else (x if k == 0 else arb(0)))] for x in rhs])

    us = []
    for k in range(dm):
        b = rhs_vec(k)
        for j in range(1, k + 1):
            b = b - M[j] * us[k - j]
        us.append(M[0].solve(b))                          # raises if not provably invertible
    if dm == 1:
        return [us[0][i, 0] for i in range(d)]
    return [Series([us[k][i, 0] for k in range(dm)], dm) for i in range(d)]


def solve_recursion(sys, eqs, coefs_init, K, dm=1, DE=None):
    """Extend ``coefs_init`` (orders 0..n0-1, each a d-list) to order K.

    Raises ZeroDivisionError if some M_{n,0} cannot be proven invertible (resonance)."""
    d = sys.d
    coefs = [list(c) for c in coefs_init]
    D, E = structure_matrices(sys, eqs, coefs, dm) if DE is None else DE
    for n in range(len(coefs), K + 1):
        cap = n + eqs.max_shift + 1
        u = series_from_coefs(coefs, d, cap=cap, dm=dm)          # u_n = 0 placeholder
        Fres = sys.residual(u, cap=cap, t=var_series(cap, dm))
        rhs = eqs.values(Fres, n)
        coefs.append(solve_level(level_matrix(D, E, n), rhs, dm))
    return coefs, (D, E)


def extract_level_matrix(sys, eqs, coefs, n, dm=1):
    """M_n by finite differences of E_n in u_n (exact: E_n is affine in u_n)."""
    d = sys.d
    cap = n + eqs.max_shift + 1
    base = coefs[:n]

    def En(un):
        u = series_from_coefs(base + [un], d, cap=cap, dm=dm)
        return eqs.values(sys.residual(u, cap=cap, t=var_series(cap, dm)), n)

    f0 = En([arb(0)] * d)
    cols = []
    for i in range(d):
        e = [arb(0)] * d
        e[i] = arb(1)
        fi = En(e)
        cols.append([fi[k] - f0[k] for k in range(d)])
    M = [[cols[l][k] for l in range(d)] for k in range(d)]
    return dmat_list(M, d, dm), f0


def level_residuals(sys, eqs, coefs, n_from=0, dm=1):
    """E_n(u-bar), n_from <= n <= K, of the truncated series (all must contain 0)."""
    d = sys.d
    K = len(coefs) - 1
    cap = K + eqs.max_shift + 2
    u = series_from_coefs(coefs, d, cap=cap, dm=dm, extra_zero=False)
    Fres = sys.residual(u, cap=cap, t=var_series(cap, dm))
    return [eqs.values(Fres, n) for n in range(n_from, K + 1)]


def contains_zero(x):
    """True if the arb, or every delta-coefficient of the Series, contains 0."""
    if isinstance(x, Series):
        return all(c.contains(arb(0)) for c in x.coeffs())
    return x.contains(arb(0))
