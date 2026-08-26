"""Tail bound for a parameter derivative y = du/dlambda from the certificate of u.

Setting (see variational.py): the augmented level equations split into the u-rows
E^u_n(u) and the y-rows E^y_n(u, y), and the y-rows are *affine in y* with
    dE^y_n / dy_n = M_n = n D + E          (the u-block's level matrix, exactly),
    d(E^y)/d(y-tail) = the u-block's off-diagonal structure (P, Psi, DQ at u),
because P~_{yy} = P, Psi~_{yy} = Psi, DQ~_{yy} = DQ.  Let v be the u-tail
(||v||_nu <= eps_u, from ``tailbound.certify_tail`` on the u-block, whose
contraction constant is Z = Z1 + Z2(eps_u) with ||M_n^{-1}|| <= c/n) and z the
unknown y-tail.  For fixed v the map
    (T z)_n = z_n - M_n^{-1} E^y_n(u-bar + v, y-bar + z),      n > K,
is affine with Lipschitz constant <= Z < 1 (same estimate as for the u-block,
since the yy-blocks coincide and the Lipschitz increments over the u-ball are
the same Banach-algebra bounds), so its unique fixed point -- the true y-tail --
satisfies  ||z||_nu <= ||T(0)||_nu / (1 - Z)  with
    ||T(0)||_nu <= Y_y + Gamma eps_u,
    Y_y   = c sum_{n>K} |E^y_n(u-bar, y-bar)|_inf nu^n / n            (finite: polynomial residual),
    Gamma = c || [ nu^{-(s+sigma)} ( P~^abs_{y,u} + eps_u sum_{i<d} dP~^abs_{y,i}/du )
                   + nu^{-s}/(K+2) ( sum_i ||theta z-bar_i||_nu dP~^abs_{y,i}/du + DQ~^abs_{y,u} ) ] ||_inf
evaluated at (nu; ||u-bar||_nu + eps_u; ||y-bar||_nu): the mean-value bound of
E^y_n(u-bar + v, y-bar) - E^y_n(u-bar, y-bar) over the u-ball (the eps_u dP term is
dP paired with theta v, where s_j <= n cancels the 1/n; the others use 1/n <= 1/(K+2)).
The y-rows may contain extra u-only combinations (the l' F_n terms of the sonic
level equations); they are covered by the same formula through their triples.

Result: |y_n|_inf <= eps_y nu^{-n} (n > K) and sum_{n>K} |y_n| r^n <= eps_y (r/nu)^{K+1}.
"""
from __future__ import annotations

from flint import arb, arb_mat

from . import recursion
from .arbseries import abs_upper, l1nu, to_arb
from .polysys import abs_eval
from .tailbound import Certificate, norm_inf


def linear_tail_certificate(sys_aug, eqs_aug, coefs_aug, d, cert_u):
    """Certificate for the y-tail (components d..2d-1 of ``coefs_aug``) given the
    u-certificate ``cert_u`` (K, nu, eps, Z1, Z2, c of the d-dimensional u-block)."""
    if not cert_u.ok:
        raise ValueError("u-certificate not valid")
    K, nu, eps_u, c = cert_u.K, cert_u.nu, cert_u.eps, cert_u.c
    Z = cert_u.Z1 + cert_u.Z2
    if not Z < 1:
        raise ValueError("u-certificate contraction constant not < 1")
    sigma = sys_aug.sigma
    D2 = sys_aug.d                                    # d u-components + d_y <= d y-components
    assert d < D2 <= 2 * d and len(coefs_aug) == K + 1
    rows_y = eqs_aug.rows[d:]
    z = recursion.series_from_coefs(coefs_aug, D2, cap=None, extra_zero=False)
    theta_z = [sys_aug.theta(zi) for zi in z]
    Fres = sys_aug.residual(z)
    nmax = max(len(f) for f in Fres) + 2
    # Y_y: residual of the truncation in the y-level rows beyond K
    Yy = arb(0)
    for n in range(K + 1, nmax + 1):
        vals = []
        for row in rows_y:
            acc = arb(0)
            for c_, r, s in row:
                if n + s >= 0:
                    acc += Fres[r][n + s] * c_
            vals.append(abs_upper(acc))
        Yy += max(vals, key=lambda a: a.mid()) * nu**n / n
    Yy = c * Yy
    # Gamma: majorant of d E^y / d(u-tail) over the u-ball, in l^1_nu operator norm
    rn = [l1nu(zi, nu) for zi in z]
    args = [nu] + [rn[i] + eps_u for i in range(d)] + [rn[i] for i in range(d, D2)]
    thn = [l1nu(t, nu) for t in theta_z]
    dP, dQ = sys_aug.dP(), sys_aug.dQ()
    inv_K2 = arb(1) / (K + 2)
    G = [[arb(0)] * d for _ in range(len(rows_y))]
    for i, row in enumerate(rows_y):
        for c_, r, s in row:
            ac = abs_upper(to_arb(c_))
            wP = nu ** (-(s + sigma))
            wQ = nu ** (-s)
            for l in range(d):
                term = wP * abs_eval(sys_aug.P[r][l], args)
                term += wQ * inv_K2 * abs_eval(dQ[r][l], args)
                for ip in range(D2):
                    if dP[r][ip][l].is_zero():
                        continue
                    dPabs = abs_eval(dP[r][ip][l], args)
                    term += wQ * inv_K2 * thn[ip] * dPabs          # dP paired with theta z-bar
                    if ip < d:
                        term += wP * eps_u * dPabs                  # dP paired with theta v (u-tail)
                G[i][l] += ac * term
    Gamma = c * norm_inf(arb_mat(G))
    eps_y = ((Yy + Gamma * eps_u) / (1 - Z)).abs_upper()
    return Certificate(K, nu, eps_y, Yy, cert_u.Z1, cert_u.Z2, c, True,
                       dict(Gamma=Gamma, eps_u=eps_u))
