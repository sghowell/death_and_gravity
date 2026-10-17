"""Independent float QNM method for KAT-7: hyperboloidal Chebyshev collocation (in the manner of
Jansen, arXiv:1709.09178).  Evidence only, never a certificate.

Shares no code with qnm_float (Frobenius series + Taylor continuation of the t-slicing ODE).
Equation: the hyperboloidal radial operator of p20_identity C1,
    E[u] = (Delta u')' + i w (2 r^2 eta u' + (r^2 eta)' u) + w^2 r^2 (1 - eta^2)/F u - L u = 0,
eta = (2r - r_+ - r_c)/(r_c - r_+), which is regular-singular at r_+ and r_c with exponents
{0, i w/kappa}: polynomial collocation without boundary conditions selects the smooth solutions.
Variable: x = log((r - r_-)/(r_+ - r_-)) in [0, X_c] (pushes the nearby singular point r_- to
x = -infinity; the remaining singularities sit at Im x = +-pi), Chebyshev-Gauss-Lobatto nodes.
Dividing E[u] by 1/s with s = r - r_- = (r_+ - r_-) e^x gives the quadratic pencil
    T(w) u = [ (Delta/s)(u_xx - u_x) + Delta' u_x - L s u ]
             + w i [ 2 r^2 eta u_x + s (2 r eta + r^2 eta') u ] + w^2 [ s * r^2 (1 - eta^2)/F ] u.
Full spectrum: numpy eig of the companion linearisation (double precision, all eigenvalues,
no seeds).  Refinement: Newton on the bordered system [T(w) T'(w)u; c^T 0] with acb_mat at
high precision (flint), starting from the double-precision eigenvalue.
"""
from __future__ import annotations

import numpy as np
from flint import acb, acb_mat, arb, ctx

from .arbutil import precision


class Collocation:
    def __init__(self, geom, l, N, prec=256):
        """geom: object with rn, rm, rp, rc, Lam, Delta(), dDelta() as mpmath values (qnm_float.Geom)."""
        self.l, self.N, self.prec = l, N, prec
        with precision(prec):
            self._build(geom)

    def _build(self, g):
        cv = lambda x: arb(str(x))
        rn, rm, rp, rc, Lam = (cv(v) for v in (g.rn, g.rm, g.rp, g.rc, g.Lam))
        self.km_arb = cv(g.km)
        N = self.N
        Xc = ((rc - rm) / (rp - rm)).log()
        xi = [(arb.pi() * j / N).cos() for j in range(N + 1)]          # 1 ... -1
        x = [Xc * (1 - t) / 2 for t in xi]                              # 0 ... Xc  (x_0 = 0 at r_+)
        cfac = [arb(2) if j in (0, N) else arb(1) for j in range(N + 1)]
        D = [[arb(0)] * (N + 1) for _ in range(N + 1)]
        for i in range(N + 1):
            for j in range(N + 1):
                if i != j:
                    D[i][j] = (cfac[i] / cfac[j]) * (-1) ** (i + j) / (xi[i] - xi[j])
            D[i][i] = -sum((D[i][j] for j in range(N + 1) if j != i), arb(0))
        # d/dx = (dxi/dx) d/dxi,  xi = 1 - 2x/Xc  ->  dxi/dx = -2/Xc
        Dx = [[D[i][j] * (-2 / Xc) for j in range(N + 1)] for i in range(N + 1)]
        Dxx = [[sum((Dx[i][k] * Dx[k][j] for k in range(N + 1)), arb(0)) for j in range(N + 1)]
               for i in range(N + 1)]
        L = self.l * (self.l + 1)
        s = [(rp - rm) * xx.exp() for xx in x]
        r = [rm + si for si in s]
        lam3 = Lam / 3
        A0, B0, C0, B1, C1, W2 = [], [], [], [], [], []
        for ri, si in zip(r, s):
            Delta_over_s = -lam3 * (ri - rn) * (ri - rp) * (ri - rc)
            dDelta = -4 * lam3 * ri ** 3 + 2 * ri - 2
            et = (2 * ri - rp - rc) / (rc - rp)
            det = 2 / (rc - rp)
            A0.append(Delta_over_s)
            B0.append(dDelta - Delta_over_s)                          # coefficient of u_x in L0
            C0.append(-L * si)
            B1.append(2 * ri * ri * et)
            C1.append(si * (2 * ri * et + ri * ri * det))
            W2.append(si * 4 * ri ** 4 / (lam3 * (rc - rp) ** 2 * (ri - rn) * (ri - rm)))
        n = N + 1
        self.L0 = [[A0[i] * Dxx[i][j] + B0[i] * Dx[i][j] + (C0[i] if i == j else 0) for j in range(n)]
                   for i in range(n)]
        self.L1 = [[B1[i] * Dx[i][j] + (C1[i] if i == j else 0) for j in range(n)] for i in range(n)]   # times i
        self.W2 = W2
        self.r_nodes = r

    # -------------------------------------------------------------------------------- spectrum
    def eig_double(self):
        """All eigenvalues (double precision) of the companion linearisation, in kappa_- units."""
        n = self.N + 1
        L0 = np.array([[float(v.mid()) for v in row] for row in self.L0])
        L1 = 1j * np.array([[float(v.mid()) for v in row] for row in self.L1])
        w2 = np.array([float(v.mid()) for v in self.W2])
        Ci = np.zeros((2 * n, 2 * n), dtype=complex)
        Ci[:n, n:] = np.eye(n)
        Ci[n:, :n] = -L0 / w2[:, None]
        Ci[n:, n:] = -L1 / w2[:, None]
        ev = np.linalg.eigvals(Ci)
        return ev / float(self.km_arb.mid())

    def _T(self, w):
        n = self.N + 1
        i = acb(0, 1)
        return acb_mat([[self.L0[a][b] + w * i * self.L1[a][b] + (w * w * self.W2[a] if a == b else 0)
                         for b in range(n)] for a in range(n)])

    def _dT(self, w):
        n = self.N + 1
        i = acb(0, 1)
        return acb_mat([[i * self.L1[a][b] + (2 * w * self.W2[a] if a == b else 0) for b in range(n)]
                        for a in range(n)])

    def refine(self, w0_over_km, iters=40, tol_bits=None):
        """Newton on the nonlinear eigenproblem; returns (omega/kappa_-, |last step|/|omega|)."""
        n = self.N + 1
        with precision(self.prec):
            km = self.km_arb
            w = acb(arb(str(w0_over_km.real)), arb(str(w0_over_km.imag))) * km
            # initial eigenvector: solve T(w) u = e (inverse iteration step)
            T = self._T(w)
            rhs = acb_mat([[acb(1)] for _ in range(n)])
            u = T.solve(rhs, nonstop=True)
            k0 = max(range(n), key=lambda k: abs(u[k, 0].mid()))
            nrm = u[k0, 0]
            u = acb_mat([[u[k, 0] / nrm] for k in range(n)])
            last = None
            tol = arb(2) ** (-(tol_bits or (self.prec - 40)))
            for _ in range(iters):
                T, dT = self._T(w), self._dT(w)
                Tu = T * u
                dTu = dT * u
                M = acb_mat([[T[a, b] for b in range(n)] + [dTu[a, 0]] for a in range(n)]
                            + [[acb(1) if b == k0 else acb(0) for b in range(n)] + [acb(0)]])
                rhs = acb_mat([[-Tu[a, 0]] for a in range(n)] + [[-(u[k0, 0] - 1)]])
                d = M.solve(rhs, nonstop=True)
                dw = d[n, 0]
                u = acb_mat([[u[k, 0] + d[k, 0]] for k in range(n)])
                w = w + dw
                # strip accumulated ball radii (this is a float method; rigor is not claimed)
                w = acb(w.real.mid(), w.imag.mid())
                u = acb_mat([[acb(u[k, 0].real.mid(), u[k, 0].imag.mid())] for k in range(n)])
                last = abs(dw.mid()) / abs(w.mid())
                if last < tol:
                    break
            res = w / km
            return complex(float(res.real.mid()), float(res.imag.mid())), res, float(last)


def spectrum_double(geom, l, N):
    return Collocation(geom, l, N, prec=128).eig_double()


def refine_mode(geom, l, N, seed_over_km, prec=256):
    with precision(prec):
        c = Collocation(geom, l, N, prec)
        ev = c.eig_double()
        k = int(np.argmin(np.abs(ev - seed_over_km)))
        w_dbl = ev[k]
        _, wball, last = c.refine(w_dbl)
        return w_dbl, wball, last


def acb_str(z, digits=25):
    old = ctx.prec
    try:
        return [z.real.mid().str(digits, radius=False), z.imag.mid().str(digits, radius=False)]
    finally:
        ctx.prec = old
