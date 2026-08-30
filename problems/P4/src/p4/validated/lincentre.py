"""Theorem B, Stage 2 (b): the regular (exponent-0) solutions of the scaled 4D linearised system
at the regular centre, as certified Frobenius series in t = e^x, for complex kappa.

Along the certified A2 centre series z(t) = (n, w, v, t) the system  P(t) theta p~ = G(t) p~
(theta = t d/dt = d/dx; P, G from ``linscaled.full_system``, kappa-graded) is Fuchsian at t = 0:
P_0 is invertible for every kappa and B_0 = P_0^{-1} G_0 has spectrum {-3, -3, 0, 0} with
rank G_0 = 2 (rows A and W of G_0 vanish identically -- checked as polynomials), so the
exponent-0 family is 2-dimensional and analytic in kappa, no positive integer is an exponent, and
    (n P_0 - G_0) p_n = sum_{k=1}^n [G_k - (n-k) P_k] p_{n-k},   G_0 p_0 = 0.
Basis (analytic in kappa): p_0 with (w~, v~) = (1, 0) resp. (0, 1) and (A^, n~) solved from rows
{A, W} of G_0 (its 2x2 minor in the columns (A^, n~) is kappa-independent and invertible).
Tail (n > K), affine contraction on l^1_nu exactly as Stage 1 / notes s2-theorem-b.md section 1.3
in Euler form:  ||M_n^{-1}|| <= c/n,  c = ||P_0^{-1}||/(1 - g/(K+1)),  g = ||P_0^{-1} G_0|| < K+1;
    Z = c sum_{k>=1} (||G_k||/(K+1) + ||P_k||) nu^k  (polynomial part along the truncated background)
      + c (||inc G||_nu/(K+1) + ||inc P||_nu)      (background tail: Banach-algebra majorants with
        ||z - z_bar||_nu <= eps_v, ||theta z - theta z_bar||_nu <= eps_theta from the A2 certificate),
    Y = sum_{K<n<=K+deg} nu^n |M_n^{-1} sum_{m<=K} B_{n,m} p_m| + Z_2 ||p_bar||_nu,
Z < 1  =>  sum_{n>K} |p_n| t^n <= eps (t/nu)^{K+1},  eps = Y/(1-Z), for every kappa in the box.
The same code runs on the kappa-derivative augmentation (d = 8, ``linscaled.augmented``): the
exponent-0 family with p^_0 = (p_0, dp_0/dkappa) gives (r_i, dr_i/dkappa) with certified tails.
"""
from __future__ import annotations

from flint import acb, acb_mat, arb

from . import linscaled
from .arbseries import abs_upper, precision
from .linstep import norm_inf
from .linsys import to_acb


def _rows_vanish_at_centre(S):
    """Rows of G(t=0, z'=0) that are identically zero as polynomials (exact)."""
    n, w, v, T, dn, dw, dv = linscaled.CTX.gens()[1:]
    zero = linscaled.CTX.constant(0)
    out = []
    for r in range(S.d):
        if all(M[r][c].compose(zero, n, w, v, zero, zero, zero, zero).is_zero() for M in S.G for c in range(S.d)):
            out.append(r)
    return out


def coefficient_series(S, ce, cap=None):
    """(P_j(t), G_j(t)) as lists over j of d x d lists of real Series along the centre expansion."""
    zs = ce.series()
    t = zs[0].zero().__class__.var(cap)
    z = list(zs) + [t]
    dz = [s.euler() for s in zs]
    return S.series_matrices(z, dz, cap=cap)


class RegularFamily:
    def __init__(self, S, ce, kappa, K=50, prec=256):
        self.S, self.ce, self.kappa, self.K = S, ce, to_acb(kappa), K
        self.d = S.d
        self.zero_rows = _rows_vanish_at_centre(S)
        with precision(prec):
            self._run()

    def _combine(self, mats):
        kp = [acb(1), self.kappa, self.kappa * self.kappa]
        return [acb_mat([[sum((kp[j] * mats[j][r][c][k] for j in range(len(mats))), acb(0))
                          for c in range(self.d)] for r in range(self.d)]) for k in range(self.deg + 1)]

    def _run(self):
        d, K = self.d, self.K
        Pser, Gser = coefficient_series(self.S, self.ce)
        self.deg = max(len(m) for M in Pser + Gser for row in M for m in row) - 1
        self.Pk, self.Gk = self._combine(Pser), self._combine(Gser)
        P0, G0 = self.Pk[0], self.Gk[0]
        rows = [r for r in range(d) if r not in self.zero_rows]
        free = [c for c in range(d) if c % 4 >= 2]                 # (w~, v~) (and their derivatives)
        solved = [c for c in range(d) if c % 4 < 2]                # (A^, n~)
        if len(rows) != len(solved):
            raise ValueError("unexpected kernel structure of G_0")
        M = acb_mat([[G0[r, c] for c in solved] for r in rows])
        self.basis = []
        for e in range(d // 4 * 0 + 2):                           # two regular directions
            rhs = acb_mat([[-G0[r, free[e]]] for r in rows])
            sol = M.solve(rhs)
            p0 = [acb(0)] * d
            for i, c in enumerate(solved):
                p0[c] = sol[i, 0]
            p0[free[e]] = acb(1)
            self.basis.append(self._recur(p0))

    def _recur(self, p0):
        d, K = self.d, self.K
        P, G = self.Pk, self.Gk
        p = [acb_mat([[x] for x in p0])]
        for n in range(1, K + 1):
            rhs = acb_mat(d, 1)
            for k in range(1, min(n, self.deg) + 1):
                rhs += (G[k] - P[k] * (n - k)) * p[n - k]
            p.append((P[0] * n - G[0]).solve(rhs))
        return p

    def certify(self, nu, prec=256):
        """Affine-contraction tail bound for every basis solution; returns (ok, eps_list, details)."""
        S, K, d = self.S, self.K, self.d
        with precision(prec):
            nu = arb(nu)
            P, G = self.Pk, self.Gk
            Pinv0 = P[0].inv()
            g = norm_inf((Pinv0 * G[0]).tolist())
            if not g < K + 1:
                return False, None, dict(g=float(g))
            c = norm_inf(Pinv0.tolist()) / (1 - g / (K + 1))
            Z1 = c * sum((abs_upper(norm_inf(G[k].tolist())) / (K + 1) + norm_inf(P[k].tolist())) * nu**k
                         for k in range(1, self.deg + 1))
            # background tail: l^1_nu majorants (T = t exactly: eps 0)
            cert = self.ce.cert
            zs = self.ce.series()
            l1 = lambda s: sum((abs_upper(s[k]) * nu**k for k in range(len(s))), arb(0))
            a = [l1(s) for s in zs] + [nu] + [l1(s.euler()) for s in zs]
            ev, eth = cert.tail_bound(nu), abs_upper(nu * cert.deriv_tail_bound(nu))
            incP, incG = S.increments(a, [ev] * 3 + [arb(0)] + [eth] * 3)
            ka = [arb(1), abs_upper(self.kappa), abs_upper(self.kappa) ** 2]
            nP = sum((ka[j] * norm_inf(incP[j]) for j in range(len(incP))), arb(0))
            nG = sum((ka[j] * norm_inf(incG[j]) for j in range(len(incG))), arb(0))
            Z2 = c * (nG / (K + 1) + nP)
            Z = Z1 + Z2
            if not Z < 1:
                return False, None, dict(Z1=float(Z1), Z2=float(Z2), c=float(c), g=float(g))
            eps = []
            for p in self.basis:
                Y1 = arb(0)
                for n in range(K + 1, K + self.deg + 1):
                    rhs = acb_mat(d, 1)
                    for k in range(n - K, min(n, self.deg) + 1):
                        rhs += (G[k] - P[k] * (n - k)) * p[n - k]
                    Y1 += norm_inf((P[0] * n - G[0]).solve(rhs).tolist()) * nu**n
                qn = sum((norm_inf(p[m].tolist()) * nu**m for m in range(K + 1)), arb(0))
                eps.append(abs_upper((Y1 + Z2 * qn) / (1 - Z)))
            self.nu, self.eps = nu, eps
            return True, eps, dict(Z1=float(Z1), Z2=float(Z2), c=float(c), g=float(g), eps=[float(e) for e in eps])

    def eval(self, x, prec=256):
        """Basis values r_i(t = e^x) (lists of d acb balls, tails included; needs certify())."""
        with precision(prec):
            t = arb(x).exp()
            q = (t / self.nu).abs_upper()
            out = []
            for p, e in zip(self.basis, self.eps):
                acc = acb_mat(self.d, 1)
                for n in reversed(range(self.K + 1)):
                    acc = acc * t + p[n]
                tb = abs_upper(e * q ** (self.K + 1))
                out.append([acc[i, 0] + acb(arb(0, tb), arb(0, tb)) for i in range(self.d)])
            return out
