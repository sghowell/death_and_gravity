"""High-order Taylor expansions at the sonic point x = 0.

The sonic point is a regular singular point of the (background and linearized)
fluid rows: row3/row4 read  M(y) z' + ... = 0  with  det M(y(0)) = 0,
z = (W, V).  With ell the left null vector of M_0, the order-n coefficients
z_n are fixed by the two scalar equations
    [Res]_{n-1} = n M_0 z_n + (known)      (rank one in z_n),
    ell . [Res]_n = 0                       (z_{n+1} drops out),
which we assemble numerically by linear extraction (evaluate the residual
series at z_n = 0, e_1, e_2) and solve by least squares.  A vanishing singular
value of that 2x2 system signals a resonance (second Frobenius exponent equal
to an integer).
"""
from __future__ import annotations

import numpy as np

from . import css
from .tps import TPS, Dual


def _null_left(q0):
    return np.array([q0["c"], -q0["a"]])


def background_series(V0, branch, K=40, return_info=False):
    """Taylor coefficients (4, K) of (A, N, W, V) at the sonic point for the given
    first-order branch (0 or 1).  Also returns the smallest singular value of the
    order-n system for each n (resonance diagnostic) if ``return_info``."""
    y0, br = css.sonic_branches(V0)
    y1 = br[branch]
    coef = np.zeros((4, K))
    coef[:, 0], coef[:, 1] = y0, y1
    ell = _null_left(css.coeffs(*y0))
    svmin = np.full(K, np.nan)

    def resid(cc, n):
        ser = [TPS(cc[i]) for i in range(4)]
        r3, r4 = css.fluid_residuals(*ser, ser[2].deriv(), ser[3].deriv())
        return np.array([r3.c[n - 1], r4.c[n - 1], ell[0] * r3.c[n] + ell[1] * r4.c[n]])

    for n in range(2, K):
        ser = [TPS(coef[i]) for i in range(4)]
        q = css.coeffs(*ser)
        coef[0, n] = (ser[0] * q["FA"]).c[n - 1] / n
        coef[1, n] = (ser[1] * q["FN"]).c[n - 1] / n
        f0 = resid(coef, n)
        cc = coef.copy(); cc[2, n] = 1.0
        f1 = resid(cc, n) - f0
        cc = coef.copy(); cc[3, n] = 1.0
        f2 = resid(cc, n) - f0
        Mn = np.column_stack([f1, f2])
        zn, _, _, sv = np.linalg.lstsq(Mn, -f0, rcond=None)
        svmin[n] = sv[-1] / sv[0] if sv[0] > 0 else np.nan
        coef[2, n], coef[3, n] = zn
    if return_info:
        return coef, svmin
    return coef


def eval_series(coef, x, deriv=False):
    """Evaluate Taylor series (4, K) (and optionally its x-derivative) at x."""
    n = np.arange(coef.shape[1])
    val = coef @ x ** n
    if not deriv:
        return val
    dval = coef[:, 1:] @ (n[1:] * x ** (n[1:] - 1))
    return val, dval


def radius_estimate(coef):
    """Crude radius-of-convergence estimate from the tail coefficients (root test)."""
    K = coef.shape[1]
    tail = np.abs(coef[:, K // 2:]).max(axis=0)
    n = np.arange(K // 2, K)
    good = tail > 0
    return float(np.exp(-np.max(np.log(tail[good]) / n[good]))) if good.any() else np.inf


# ----------------------------------------------------------------------------
# linear perturbations  h = H_ss(x) + eps h_p(x) e^{kappa s}
# ----------------------------------------------------------------------------
def _lin_residuals(bg, pc, kappa):
    """Series of the linearized rows (3, 4), of d/dx-rows (1, 2) and of the
    linearized momentum constraint, given background series ``bg`` (4, K) and
    perturbation series ``pc`` (4, K, complex)."""
    ser = [Dual(TPS(bg[i]), TPS(pc[i])) for i in range(4)]
    Wx = Dual(TPS(bg[2]).deriv(), TPS(pc[2]).deriv())
    Vx = Dual(TPS(bg[3]).deriv(), TPS(pc[3]).deriv())
    Ws = Dual(TPS(np.zeros_like(bg[2])), kappa * TPS(pc[2]))
    Vs = Dual(TPS(np.zeros_like(bg[3])), kappa * TPS(pc[3]))
    r3, r4 = css.fluid_residuals(*ser, Wx, Vx, Ws, Vs)
    q = css.coeffs(*ser)
    dA = ser[0] * q["FA"]                       # A_x  (dual)
    dN = ser[1] * q["FN"]                       # N_x  (dual)
    cons = ser[0] * (q["G"] - q["FA"])          # A_s = A (G - F_A); linearized: kappa A_p
    return r3.b, r4.b, dA.b, dN.b, cons.b - kappa * TPS(pc[0])


def perturbation_series(bg, kappa, K=None, return_info=False):
    """Taylor coefficients (4, K) of the perturbation (A_p, N_p, W_p, V_p), analytic
    at the sonic point, gauge N_p(0) = 0, normalised A_p(0) = 1, and satisfying the
    linearized momentum constraint.  ``bg`` = background series (4, K_bg)."""
    K = bg.shape[1] if K is None else K
    bg = bg[:, :K]
    pc = np.zeros((4, K), dtype=complex)
    y0 = bg[:, 0]
    ell = _null_left(css.coeffs(*y0))

    def cond0(w, v):
        cc = pc.copy(); cc[0, 0], cc[1, 0], cc[2, 0], cc[3, 0] = 1.0, 0.0, w, v
        r3, r4, _, _, cons = _lin_residuals(bg, cc, kappa)
        return np.array([ell[0] * r3.c[0] + ell[1] * r4.c[0], cons.c[0]])

    f0 = cond0(0, 0)
    M0 = np.column_stack([cond0(1, 0) - f0, cond0(0, 1) - f0])
    pc[0, 0], pc[1, 0] = 1.0, 0.0
    pc[2, 0], pc[3, 0] = np.linalg.solve(M0, -f0)
    svmin = np.full(K, np.nan)

    def resid(cc, n):
        r3, r4, _, _, _ = _lin_residuals(bg, cc, kappa)
        return np.array([r3.c[n - 1], r4.c[n - 1], ell[0] * r3.c[n] + ell[1] * r4.c[n]])

    for n in range(1, K):
        _, _, dA, dN, _ = _lin_residuals(bg, pc, kappa)
        pc[0, n] = dA.c[n - 1] / n
        pc[1, n] = dN.c[n - 1] / n
        f0 = resid(pc, n)
        cc = pc.copy(); cc[2, n] = 1.0
        f1 = resid(cc, n) - f0
        cc = pc.copy(); cc[3, n] = 1.0
        f2 = resid(cc, n) - f0
        Mn = np.column_stack([f1, f2])
        zn, _, _, sv = np.linalg.lstsq(Mn, -f0, rcond=None)
        svmin[n] = sv[-1] / sv[0] if sv[0] > 0 else np.nan
        pc[2, n], pc[3, n] = zn
    if return_info:
        _, _, _, _, cons = _lin_residuals(bg, pc, kappa)
        return pc, svmin, cons.c
    return pc


def resonant_kappas(bg, nmax=12):
    """kappa values at which the order-n perturbation system is singular
    (second Frobenius exponent sigma(kappa) = n): det_n(kappa) is affine in kappa."""
    return np.array([_det_root(bg, n) for n in range(1, nmax + 1)])


def _det_root(bg, n):
    """Root in kappa of the (affine) determinant of the order-n 2x2 system."""
    def det(kap):
        K = n + 2
        pc = perturbation_series(bg, kap, K=K)      # fills orders < n correctly
        ell = _null_left(css.coeffs(*bg[:, 0]))
        cc = pc.copy(); cc[2, n:] = 0; cc[3, n:] = 0
        def resid(c):
            r3, r4, _, _, _ = _lin_residuals(bg[:, :K], c, kap)
            return np.array([r3.c[n - 1], r4.c[n - 1], ell[0] * r3.c[n] + ell[1] * r4.c[n]])
        f0 = resid(cc)
        c1 = cc.copy(); c1[2, n] = 1.0
        c2 = cc.copy(); c2[3, n] = 1.0
        Mn = np.column_stack([resid(c1) - f0, resid(c2) - f0])
        # use the two rows (row of M_0 with larger norm, ell-row)
        i = 0 if abs(Mn[0]).sum() > abs(Mn[1]).sum() else 1
        return Mn[i, 0] * Mn[2, 1] - Mn[i, 1] * Mn[2, 0]
    d0, d1 = det(0.0), det(1.0)
    return -d0 / (d1 - d0)
