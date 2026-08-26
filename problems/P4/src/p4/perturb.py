"""Stage C: linear perturbations  h(s,x) = H_ss(x) + eps h_p(x) e^{kappa s}  (KHA95 sec. 5).

The time-dependent KHA system is linearized exactly with dual numbers
(``p4.tps.Dual``) around the CSS background; the linearized momentum constraint
    kappa A_p = dC . y_p,   C(y) := A (G - F_A)   (so that A_s = C),
is solved algebraically for A_p, and the reduced linear system for
(N_p, W_p, V_p) is co-integrated with the (constraint-reduced, centre-scaled)
background from the sonic point toward the centre for a whole batch of kappa
values at once.  Sonic-point data: analytic (Taylor) with gauge N_p(0) = 0,
normalisation A_p(0) = 1.  Matching function
    E(kappa) = A_p(x_end) e^{x_end}      (amplitude of the irregular
                                          central-mass mode A_p ~ e^{-x});
its zeros are the eigenvalues (in the sonic-point gauge).
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from . import css, taylor
from .tps import Dual


def _split(u, nk):
    """State vector (or dense (n, nt) array) -> complex (N_p, W_p, V_p), each (nk[, nt])."""
    p = u[3:].reshape((3, 2, nk) + u.shape[1:])
    return p[0, 0] + 1j * p[0, 1], p[1, 0] + 1j * p[1, 1], p[2, 0] + 1j * p[2, 1]


def linearized_rhs(x, u, kappas):
    """RHS of background (scaled, reduced) + perturbations (plain, complex) for all kappas."""
    nk = len(kappas)
    Nh, Wh, Vh = u[0], u[1], u[2]
    e_x = np.exp(x)
    e_2x = e_x * e_x
    N, W, V = Nh / e_x, Wh * e_2x, Vh * e_x
    Am1 = css.Am1_constraint(N, W, V)
    A = 1 + Am1
    q = css.coeffs(A, N, W, V, Am1=Am1)
    Delta, P, Q = css.det_and_numerators(q)
    Wx, Vx = P / Delta, Q / Delta
    dbg = [Nh * (Am1 - 2 * W / 3), Wx / e_2x - 2 * Wh, Vx / e_x - Vh]
    Np, Wp, Vp = _split(u, nk)
    # linearized constraint  kappa A_p = C_A A_p + C_rest
    yA = [Dual(A, 1.0), Dual(N, 0.0), Dual(W, 0.0), Dual(V, 0.0)]
    qA = css.coeffs(*yA, Am1=Dual(Am1, 1.0))
    C_A = (yA[0] * (qA["G"] - qA["FA"])).b
    z = np.zeros(nk, dtype=complex)
    yR = [Dual(A, z), Dual(N, Np), Dual(W, Wp), Dual(V, Vp)]
    qR = css.coeffs(*yR, Am1=Dual(Am1, z))
    C_R = (yR[0] * (qR["G"] - qR["FA"])).b
    Ap = C_R / (kappas - C_A)
    # linearized fluid rows with z_p' = 0, then solve M z_p' = -(R3, R4)
    yp = [Dual(A, Ap), Dual(N, Np), Dual(W, Wp), Dual(V, Vp)]
    qp = css.coeffs(*yp, Am1=Dual(Am1, Ap))
    Ws, Vs = Dual(0.0, kappas * Wp), Dual(0.0, kappas * Vp)
    R3 = (qp["sa"] * Ws + qp["sb"] * Vs + qp["a"] * Wx + qp["b"] * Vx + qp["e"]).b
    R4 = (qp["sc"] * Ws + qp["sd"] * Vs + qp["c"] * Wx + qp["d"] * Vx + qp["f"]).b
    dWp = (-R3 * q["d"] + R4 * q["b"]) / Delta
    dVp = (-q["a"] * R4 + q["c"] * R3) / Delta
    dNp = (yp[1] * qp["FN"]).b
    out = np.empty(3 + 6 * nk)
    out[:3] = dbg
    out[3:] = np.concatenate([dNp.real, dNp.imag, dWp.real, dWp.imag, dVp.real, dVp.imag])
    return out


def A_p_from_constraint(x, yh3, Np, Wp, Vp, kappas):
    """A_p at one point from the linearized momentum constraint."""
    e_x = np.exp(x)
    N, W, V = yh3[0] / e_x, yh3[1] * e_x**2, yh3[2] * e_x
    Am1 = css.Am1_constraint(N, W, V)
    A = 1 + Am1
    yA = [Dual(A, 1.0), Dual(N, 0.0), Dual(W, 0.0), Dual(V, 0.0)]
    qA = css.coeffs(*yA, Am1=Dual(Am1, 1.0))
    C_A = (yA[0] * (qA["G"] - qA["FA"])).b
    z = np.zeros_like(np.asarray(Np, dtype=complex))
    yR = [Dual(A, z), Dual(N, Np), Dual(W, Wp), Dual(V, Vp)]
    qR = css.coeffs(*yR, Am1=Dual(Am1, z))
    C_R = (yR[0] * (qR["G"] - qR["FA"])).b
    return C_R / (kappas - C_A)


class Problem:
    """Eigenvalue problem on the EC background (V0, branch)."""

    def __init__(self, V0, branch=1, K=36, delta=None, x_end=-9.0, rtol=1e-12, atol=1e-14):
        self.V0, self.branch, self.K = V0, branch, K
        self.bg = taylor.background_series(V0, branch, K=K)
        R = taylor.radius_estimate(self.bg)
        self.delta = float(np.clip(R / 6, 1e-3, 0.1)) if delta is None else delta
        self.x_end, self.rtol, self.atol = x_end, rtol, atol
        self.y_start = taylor.eval_series(self.bg, -self.delta)      # plain (A, N, W, V)

    def sonic_perturbation(self, kappa):
        """(N_p, W_p, V_p) at x = -delta and the Taylor coefficients."""
        pc = taylor.perturbation_series(self.bg, kappa, K=self.K)
        yp = taylor.eval_series(pc, -self.delta)
        return yp[1:], pc

    def integrate(self, kappas, x_end=None, dense=False):
        kappas = np.atleast_1d(np.asarray(kappas, dtype=complex))
        nk = len(kappas)
        x_end = self.x_end if x_end is None else x_end
        yh0 = css.to_scaled(-self.delta, self.y_start)[1:]
        pstart = np.array([self.sonic_perturbation(k)[0] for k in kappas])   # (nk, 3)
        u0 = np.concatenate([yh0, pstart[:, 0].real, pstart[:, 0].imag, pstart[:, 1].real,
                             pstart[:, 1].imag, pstart[:, 2].real, pstart[:, 2].imag])
        sol = solve_ivp(linearized_rhs, (-self.delta, x_end), u0, args=(kappas,), method="DOP853",
                        rtol=self.rtol, atol=self.atol, dense_output=dense)
        if sol.status != 0:
            raise RuntimeError(sol.message)
        return sol, kappas

    def E(self, kappas, x_end=None):
        """Matching function E(kappa) = A_p(x_end) e^{x_end} (vectorised over kappas)."""
        sol, kappas = self.integrate(kappas, x_end)
        x = sol.t[-1]
        Np, Wp, Vp = _split(sol.y[:, -1], len(kappas))
        Ap = A_p_from_constraint(x, sol.y[:3, -1], Np, Wp, Vp, kappas)
        return Ap * np.exp(x)

    def perturbation_profile(self, kappa, x_end=None):
        """Dense (x, background scaled, A_p, N_p, W_p, V_p) for one kappa."""
        sol, kappas = self.integrate([kappa], x_end, dense=True)
        Np, Wp, Vp = _split(sol.y, 1)
        Ap = np.array([A_p_from_constraint(sol.t[i], sol.y[:3, i], Np[0, i], Wp[0, i], Vp[0, i],
                                           kappas)[0] for i in range(len(sol.t))])
        return sol.t, sol.y[:3], Ap, Np[0], Wp[0], Vp[0]


def centre_exponents(prob, kappa, x=-8.0):
    """Eigenvalues of the constant-coefficient limit of the scaled perturbation system
    (n, w, v) = (N_p e^x, W_p e^{-2x}, V_p e^{-x}) as x -> -inf, on the asymptotic
    background (Nh, Wh, Vh)_inf taken from a regular shot."""
    from . import shoot
    sh = shoot.shoot(prob.V0, prob.branch, x_end=-8.0, keep_sol=True)
    Nh, Wh, Vh = sh.sol.sol(-8.0)
    kap = np.array([kappa], dtype=complex)
    e_x = np.exp(x)
    S = np.array([1 / e_x, e_x**2, e_x])                      # plain = S * scaled
    def f(nwv):
        plain = S * nwv
        u = np.concatenate([[Nh, Wh, Vh], [plain[0].real], [plain[0].imag], [plain[1].real],
                            [plain[1].imag], [plain[2].real], [plain[2].imag]])
        d = linearized_rhs(x, u, kap)
        dp = np.array([d[3] + 1j * d[4], d[5] + 1j * d[6], d[7] + 1j * d[8]])
        return dp / S - np.array([-1, 2, 1]) * nwv              # d/dx of scaled variables
    J = np.column_stack([f(np.eye(3)[:, i].astype(complex)) for i in range(3)])
    return np.linalg.eigvals(J), J
