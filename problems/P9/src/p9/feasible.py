"""Feasible (primal) points with large H0: a lower bound on the true maximum over F.

Local nonconvex optimization of min u_0 subject to chi2(u) <= T and the class constraints;
the returned point is projected onto the class and its chi2 re-evaluated exactly, so
"feasible" means feasible for the true statement, not for a relaxation."""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from . import C_KM_S
from .classmin import chi2_and_grad, in_class, project_to_class
from .model import Frozen


def max_H0_point(fr: Frozen, T: float, u0: np.ndarray, maxiter: int = 1000, verbose: bool = True):
    """Returns (u, Mp, chi2, H0_at_r_lo) for a class point with chi2 <= T and small u_0."""
    sp = fr.spec
    N = sp.n_seg
    ulo, uhi = sp.u_box
    cons = []
    if np.isfinite(sp.L):
        Lh = sp.L * sp.hs
        M = np.zeros((4 * N, N + 1))
        for k in range(N):
            M[k, k] = Lh[k] + 1.0; M[k, k + 1] = -1.0
            M[N + k, k] = 1.0; M[N + k, k + 1] = Lh[k] - 1.0
            M[2 * N + k, k] = Lh[k] - 1.0; M[2 * N + k, k + 1] = 1.0
            M[3 * N + k, k] = -1.0; M[3 * N + k, k + 1] = Lh[k] + 1.0
        cons.append(dict(type="ineq", fun=lambda u: M @ u, jac=lambda u: M))
    # chi2 constraint with a small safety margin so that the projected point stays feasible
    margin = 1e-3
    cons.append(dict(type="ineq", fun=lambda u: T - margin - chi2_and_grad(fr, u)[0],
                     jac=lambda u: -chi2_and_grad(fr, u)[1]))
    e0 = np.zeros(N + 1); e0[0] = 1.0
    res = minimize(lambda u: u[0], project_to_class(fr, u0), jac=lambda u: e0, method="SLSQP",
                   constraints=cons, bounds=[(ulo, uhi)] * (N + 1), options=dict(maxiter=maxiter, ftol=1e-12))
    u = project_to_class(fr, res.x)
    chi2, _, Mp = chi2_and_grad(fr, u)
    ok = in_class(fr, u) and chi2 <= T
    if verbose:
        print(f"  feasible-point search: {res.message}; u0={u[0]:.5f} chi2={chi2:.4f} (T={T:.4f}) "
              f"in_class={in_class(fr, u)} feasible={ok} -> H0 = {C_KM_S/(sp.r_lo*u[0]):.3f} at r_lo", flush=True)
    return u, float(Mp), chi2, (C_KM_S / (sp.r_lo * u[0]) if ok else None)
