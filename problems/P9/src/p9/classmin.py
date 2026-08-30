"""True (non-relaxed) chi2 minimization over the class: a smooth nonconvex problem in the
node values u, solved locally (SLSQP with analytic gradients) and projected back onto the class
so that the reference point is exactly feasible. Any class point is a valid reference for T;
a better one only tightens the bound."""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from .model import Frozen

LN10 = np.log(10.0)


def chi2_and_grad(fr: Frozen, u: np.ndarray):
    """chi2 with Mp profiled analytically, and its gradient wrt u."""
    D = fr.A_sn @ u
    y = fr.sn.m - 5.0 * np.log10(D) - fr.sn_offset
    ones = np.ones_like(y)
    Cinv_y = fr.Wsn.T @ (fr.Wsn @ y)
    Cinv_1 = fr.Wsn.T @ (fr.Wsn @ ones)
    Mp = (Cinv_1 @ y) / (Cinv_1 @ ones)
    r = y - Mp
    Cinv_r = fr.Wsn.T @ (fr.Wsn @ r)          # profiling makes d chi2/d Mp = 0, so the envelope theorem applies
    pred, J = fr.bao_pred(u, jac=True)          # J = P except on D_V rows (nonlinear)
    rb = fr.bao.value - pred
    Cb_rb = fr.Wb.T @ (fr.Wb @ rb)
    chi2 = float(r @ Cinv_r + rb @ Cb_rb)
    # d r_j / d u = -(5/ln10) A_j / D_j
    grad = -2.0 * (5.0 / LN10) * (fr.A_sn.T @ (Cinv_r / D)) - 2.0 * (J.T @ Cb_rb)
    return chi2, grad, Mp


def project_to_class(fr: Frozen, u: np.ndarray, shrink: float = 1e-9) -> np.ndarray:
    """Forward sweep enforcing |u_{k+1}-u_k| <= L h_k min(u_k,u_{k+1}) and the box, exactly."""
    sp = fr.spec
    ulo, uhi = sp.u_box
    v = np.clip(u.copy(), ulo, uhi)
    if np.isfinite(sp.L):
        for k in range(sp.n_seg):
            Lh = sp.L * sp.hs[k] * (1.0 - shrink)
            hi = v[k] * (1.0 + Lh)           # then min(u_k,u_{k+1}) = u_k
            lo = v[k] / (1.0 + Lh)           # then min = u_{k+1}: u_k - u_{k+1} <= Lh u_{k+1}
            v[k + 1] = np.clip(v[k + 1], lo, hi)
    return v


def in_class(fr: Frozen, u: np.ndarray, tol: float = 1e-12) -> bool:
    sp = fr.spec
    ulo, uhi = sp.u_box
    if np.any(u < ulo - tol) or np.any(u > uhi + tol):
        return False
    if np.isfinite(sp.L):
        d = np.abs(np.diff(u))
        return bool(np.all(d <= sp.L * sp.hs * np.minimum(u[:-1], u[1:]) + tol))
    return True


def minimize_chi2_over_class(fr: Frozen, u0: np.ndarray, maxiter: int = 500, verbose: bool = False):
    """Local minimization of the true chi2 over the class. Returns (u*, Mp*, chi2*) with u* in class."""
    sp = fr.spec
    N = sp.n_seg
    ulo, uhi = sp.u_box
    cons = []
    if np.isfinite(sp.L):
        Lh = sp.L * sp.hs
        # d_k - Lh u_k >= 0 etc. written as g(u) >= 0
        def g(u):
            d = u[1:] - u[:-1]
            return np.concatenate([Lh * u[:-1] - d, Lh * u[1:] - d, Lh * u[:-1] + d, Lh * u[1:] + d])

        M = np.zeros((4 * N, N + 1))
        for k in range(N):
            M[k, k] = Lh[k] + 1.0; M[k, k + 1] = -1.0
            M[N + k, k] = 1.0; M[N + k, k + 1] = Lh[k] - 1.0
            M[2 * N + k, k] = Lh[k] - 1.0; M[2 * N + k, k + 1] = 1.0
            M[3 * N + k, k] = -1.0; M[3 * N + k, k + 1] = Lh[k] + 1.0
        cons.append(dict(type="ineq", fun=g, jac=lambda u: M))
    f = lambda u: chi2_and_grad(fr, u)[0]
    jac = lambda u: chi2_and_grad(fr, u)[1]
    res = minimize(f, project_to_class(fr, u0), jac=jac, method="SLSQP", constraints=cons,
                   bounds=[(ulo, uhi)] * (N + 1), options=dict(maxiter=maxiter, ftol=1e-12))
    u = project_to_class(fr, res.x)
    assert in_class(fr, u)
    chi2, _, Mp = chi2_and_grad(fr, u)
    if verbose:
        print(f"  class-min: SLSQP {res.message} nit={res.nit} chi2={chi2:.4f} (start {f(u0):.4f})")
    return u, float(Mp), chi2
