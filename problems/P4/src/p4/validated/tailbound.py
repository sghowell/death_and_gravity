"""Rigorous geometric tail bound for the coefficient recursion (Banach fixed point).

Setting (see recursion.py): the exact coefficients u_n, n > K, are the unique
solution of  M_n u_n = -R_n(u_{<n}),  M_n = n D + E.  Write u = u-bar + v with
u-bar the computed truncation (balls, orders <= K, zero beyond) and v the tail.
On X = { v : ||v||_nu := sum_{n>K} |v_n|_inf nu^n < inf }  define
    (T v)_n = v_n - M_n^{-1} E_n(u-bar + v),     n > K,
whose fixed points are exactly the tails solving the recursion.  With
    ||M_n^{-1}||_inf <= c/n              (c = ||D^{-1}|| / (1 - ||D^{-1}E||/(K+1)),  n >= K+1),
    Y  >= ||T(0)||_nu = sum_{n>K} ||M_n^{-1} E_n(u-bar)|| nu^n,
    Z(eps) >= sup_{||v|| <= eps} ||DT(v)||   (operator norm on X),
the conditions  Z(eps) < 1  and  Y + Z(eps) eps <= eps  make T a contraction
of the closed ball B_eps into itself (mean value inequality), so the true tail
satisfies  sum_{n>K} |u_n|_inf nu^n <= eps,  i.e.
    |u_n|_inf <= eps nu^{-n}   (n > K)   and   sum_{n>K} |u_n| r^n <= eps (r/nu)^{K+1}   (r <= nu).

DT(v) = -M^{-1} B(u-bar + v) with the off-diagonal blocks (m < n)
    B_{nm} = sum_(c,r,s) c [ s_{m-sigma} P^{(r)}_{n-m+s+sigma} + Psi^{(r)}_{n-m+s} - DQ^{(r)}_{n-m+s} ],
bounded columnwise in l^1_nu using s_{m-sigma}/n <= 1 and 1/n <= 1/(K+2):
    Z_1 = c sum_{k>=1} || B~_k ||_inf nu^k     (u-bar part; B~_k entrywise bounds, finite sum),
    Z_2(eps) = c || [ nu^{-(s+sigma)} inc(P) + nu^{-s}/(K+2) (sum_i ||theta u-bar_i|| inc(dP_i) + inc(DQ))
                      + nu^{-(s+sigma)} eps sum_i |dP_i|^abs(|u-bar| + eps) ] ||_inf,
where inc(p) = p^abs(|u-bar|_nu + eps) - p^abs(|u-bar|_nu) (Banach-algebra Lipschitz bounds)
and the last term is the theta v contribution (its factor s_j <= n cancels the 1/n).
"""
from __future__ import annotations

from dataclasses import dataclass, field

from flint import arb, arb_mat

from . import recursion
from .arbseries import abs_upper, l1nu, to_arb
from .polysys import abs_eval, abs_eval_increment

NU_CANDIDATES = ("0.5", "0.4", "0.3", "0.25", "0.2", "0.16", "0.13", "0.1", "0.08",
                 "0.065", "0.05", "0.04", "0.03", "0.02", "0.015", "0.01", "0.005", "0.002", "0.001")


def norm_inf(M):
    """Max row sum of |entries| (upper bound) of an arb_mat."""
    best = arb(0)
    for i in range(M.nrows()):
        s = arb(0)
        for j in range(M.ncols()):
            s += abs_upper(M[i, j])
        best = best.max(s)
    return best


@dataclass
class Certificate:
    K: int
    nu: arb
    eps: arb
    Y: arb
    Z1: arb
    Z2: arb
    c: arb
    ok: bool
    details: dict = field(default_factory=dict)

    def coef_bound(self, n):
        """|u_n|_inf <= eps nu^{-n} for n > K."""
        return self.eps * self.nu ** (-n)

    def tail_bound(self, r):
        """sum_{n>K} |u_n|_inf r^n <= eps (r/nu)^{K+1}  for 0 <= r <= nu."""
        q = to_arb(r) / self.nu
        if not q <= 1:
            raise ValueError("tail bound only valid for |x| <= nu")
        return (self.eps * q ** (self.K + 1)).abs_upper()

    def deriv_tail_bound(self, r):
        """sum_{n>K} n |u_n|_inf r^{n-1} <= (eps/nu) q^K (K+1 - K q)/(1-q)^2,  q = r/nu < 1."""
        q = to_arb(r) / self.nu
        if not q < 1:
            raise ValueError("derivative tail bound needs |x| < nu")
        K = self.K
        return (self.eps / self.nu * q**K * (K + 1 - K * q) / (1 - q) ** 2).abs_upper()

    def __repr__(self):
        return (f"Certificate(K={self.K}, nu={float(self.nu):.4g}, eps={float(self.eps):.3e}, "
                f"Y={float(self.Y):.3e}, Z1={float(self.Z1):.3f}, Z2={float(self.Z2):.3e}, "
                f"c={float(self.c):.3f}, ok={self.ok})")


class _TailData:
    """nu-independent precomputation: untruncated series of the truncation u-bar."""

    def __init__(self, sys, eqs, coefs, D, E):
        d, K = sys.d, len(coefs) - 1
        self.sys, self.eqs, self.d, self.K = sys, eqs, d, K
        self.sigma = sys.sigma
        D = D[0] if isinstance(D, list) else D
        E = E[0] if isinstance(E, list) else E
        Dinv = D.inv()                                    # raises if not provably invertible
        self.g = norm_inf(Dinv * E)
        self.Dinv_norm = norm_inf(Dinv)
        self.u = recursion.series_from_coefs(coefs, d, cap=None, extra_zero=False)
        self.theta_u = [sys.theta(ui) for ui in self.u]
        self.Pser = sys.P_series(self.u)
        self.Psi = sys.psi_series(self.u)
        self.dQser = sys.dQ_series(self.u)
        self.Fres = sys.residual(self.u)
        nmax = max(len(f) for f in self.Fres) + 2
        self.Evals = []                                   # |E_n(u-bar)|_inf for n > K
        for n in range(K + 1, nmax + 1):
            vals = eqs.values(self.Fres, n)
            self.Evals.append(max((abs_upper(v) for v in vals), key=lambda a: a.mid()))
        # entrywise bound matrices B~_k, k >= 1
        kmax = max(len(s) for row in self.Pser + self.Psi + self.dQser for s in row) + 2
        self.Bt = []
        inv_K2 = arb(1) / (K + 2)
        for k in range(1, kmax + 1):
            Bk = [[arb(0)] * d for _ in range(d)]
            for i, row in enumerate(eqs.rows):
                for c, r, s in row:
                    ac = abs_upper(to_arb(c))
                    for l in range(d):
                        Bk[i][l] += ac * (abs_upper(self.Pser[r][l][k + s + self.sigma])
                                          + (abs_upper(self.Psi[r][l][k + s])
                                             + abs_upper(self.dQser[r][l][k + s])) * inv_K2)
            self.Bt.append(Bk)

    def c_bound(self):
        if not self.g < self.K + 1:
            return None
        return self.Dinv_norm / (1 - self.g / (self.K + 1))

    def Y(self, nu, c):
        tot = arb(0)
        for j, ev in enumerate(self.Evals):
            n = self.K + 1 + j
            tot += ev * nu**n / n
        return c * tot

    def Z1(self, nu, c):
        tot = arb(0)
        for k, Bk in enumerate(self.Bt, start=1):
            tot += norm_inf(arb_mat(Bk)) * nu**k
        return c * tot

    def Z2(self, nu, c, eps):
        d, K, sigma = self.d, self.K, self.sigma
        sys = self.sys
        rn = [l1nu(ui, nu) for ui in self.u]
        args = [nu] + rn
        args_eps = [nu] + [r + eps for r in rn]
        thn = [l1nu(t, nu) for t in self.theta_u]
        dP, dQ = sys.dP(), sys.dQ()
        inv_K2 = arb(1) / (K + 2)
        Z = [[arb(0)] * d for _ in range(d)]
        for i, row in enumerate(self.eqs.rows):
            for c_, r, s in row:
                ac = abs_upper(to_arb(c_))
                wP = nu ** (-(s + sigma))
                wQ = nu ** (-s)
                for l in range(d):
                    term = wP * abs_eval_increment(sys.P[r][l], args, eps)
                    term += wQ * inv_K2 * abs_eval_increment(dQ[r][l], args, eps)
                    for ip in range(d):
                        if dP[r][ip][l].is_zero():
                            continue
                        term += wQ * inv_K2 * thn[ip] * abs_eval_increment(dP[r][ip][l], args, eps)
                        term += wP * eps * abs_eval(dP[r][ip][l], args_eps)
                    Z[i][l] += ac * term
        return c * norm_inf(arb_mat(Z))


def certify_tail(sys, eqs, coefs, D, E, nu=None, margins=(1.5, 2, 4, 10, 100), verbose=False):
    """Certificate for the tail of ``coefs`` (balls up to K); largest successful nu if nu is None."""
    data = _TailData(sys, eqs, coefs, D, E)
    c = data.c_bound()
    K = data.K
    if c is None:
        return Certificate(K, arb(0), arb(0), arb(0), arb(0), arb(0), arb(0), False,
                           dict(reason=f"||D^-1 E|| = {data.g} not < K+1"))
    cands = [to_arb(nu)] if nu is not None else [to_arb(s) for s in NU_CANDIDATES]
    last = None
    for nu_ in cands:
        Y, Z1 = data.Y(nu_, c), data.Z1(nu_, c)
        if verbose:
            print(f"nu={float(nu_):.4g}: Y={float(Y):.3e} Z1={float(Z1):.4f}")
        if not Z1 < 1:
            last = Certificate(K, nu_, arb(0), Y, Z1, arb(0), c, False, dict(reason="Z1 >= 1"))
            continue
        for mg in margins:
            eps = Y / (1 - Z1) * mg
            Z2 = data.Z2(nu_, c, eps)
            Z = Z1 + Z2
            if Z < 1 and (Y + Z * eps) < eps:
                return Certificate(K, nu_, eps, Y, Z1, Z2, c, True,
                                   dict(margin=mg, Dinv_norm=data.Dinv_norm, g=data.g))
        last = Certificate(K, nu_, eps, Y, Z1, Z2, c, False, dict(reason="Z2 too large"))
    return last
