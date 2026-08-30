"""Theorem B, Stage 2: one validated Taylor step of a linear system  P(s) y' = G(s) y  (complex
kappa) along a certified background step, and complex Lohner sets for its solutions.

Step x_k -> x_k - h (s in [-h, 0]).  Background data (``lintube.StepData``, kappa-independent):
the truncated scaled background z_K(s) (exact polynomial through the reference point), the tube
|z(s) - z_K(s)| <= eps_z, |z'(s) - z_K'(s)| <= eps_dz for every solution in the tube, the
kappa-graded coefficient matrices as polynomials in s (P_j(s), G_j(s), truncated at degree K with
coefficient-sum tails dP_j, dG_j), their box values on nsub sub-intervals and the increment
majorants inc P_j, inc G_j over the tube.  For a fixed kappa:
 1. Y_K(s) = sum_{i<=K} Y_i s^i, Y_0 = I, from the recursion (i+1) P_0 Y_{i+1} = sum G_i Y - sum P_i (.) Y;
 2. defect of Y_K along the *true* background: E = Y_K' - A(z) Y_K = P(z)^{-1} [ D_K(s) + (P(z) - P(z_K)) Y_K'
    - (G(z,z') - G(z_K,z_K')) Y_K ],  D_K = P(z_K) Y_K' - G(z_K, z_K') Y_K an exact polynomial (ball
    coefficients: orders < K are rounding only, orders >= K the truncation), so
    sup|E| <= Pinv [ Dsup + dP Y'sup + dG Ysup + |inc P| Y'sup + |inc G| Ysup ],
    Pinv, L = sup |P^{-1}|, sup |P^{-1} G| over the sub-boxes (infinity norms);
 3. Groenwall: ||Y(-h) - Y_K(-h)|| <= h e^{L h} sup|E| =: bound  (operator norm, valid for every
    solution: y(-h) = Y_K(-h) y(0) + e, |e|_inf <= bound |y(0)|_inf);
 4. Lohner update  c + A r  ->  Y_K(-h)(c + A r) + box(bound |y|_inf),  complex QR re-orthogonalisation.
"""
from __future__ import annotations

import numpy as np
from flint import acb, acb_mat, acb_poly, arb

from .arbseries import abs_upper


def _cbox(rad):
    return acb(arb(0, rad), arb(0, rad))


def _rad_box(z):
    """acb box centred at 0 containing the ball z - mid(z)."""
    return acb(arb(0, abs_upper(z.real - z.real.mid())), arb(0, abs_upper(z.imag - z.imag.mid())))


def norm_inf(rows):
    """Max row sum of upper bounds of |entries| (rows: lists of acb/arb)."""
    best = arb(0)
    for row in rows:
        best = best.max(sum((abs_upper(e) for e in row), arb(0)))
    return best


def sup_poly(p, h):
    """sum_k |p_k| h^k >= sup_{|s|<=h} |p(s)| for an acb_poly."""
    tot, hp = arb(0), arb(1)
    for c in p.coeffs():
        tot += abs_upper(c) * hp
        hp *= h
    return tot


class LohnerSet:
    """{c + A r : |Re r_j|, |Im r_j| <= rho_j} in C^d (c exact midpoints, A acb_mat, r acb boxes)."""

    def __init__(self, c, radii=None, A=None):
        self.d = len(c)
        self.c = [acb(z.real.mid(), z.imag.mid()) for z in c]
        self.A = A if A is not None else acb_mat(self.d, self.d, [acb(int(i == j)) for i in range(self.d) for j in range(self.d)])
        self.r = [(_rad_box(z) if radii is None else _cbox(radii[j])) for j, z in enumerate(c)]

    @classmethod
    def from_balls(cls, balls):
        return cls(balls)

    def hull(self):
        v = self.A * acb_mat([[r] for r in self.r])
        return [self.c[i] + v[i, 0] for i in range(self.d)]

    def max_abs(self):
        return max((abs_upper(z) for z in self.hull()), default=arb(0))

    def propagate(self, J, bound, extra=None, Sc=None):
        """New set containing J (c + A r) + e, |e_i| <= bound * max_j Sc_j |y_j| / Sc_i  (+ extra boxes);
        Sc: weights of the norm in which ``bound`` holds (default: plain infinity norm)."""
        d = self.d
        Sc = Sc or [1.0] * d
        ymax = max((abs_upper(z) * Sc[i] for i, z in enumerate(self.hull())), default=arb(0))
        C = J * self.A
        Cm = np.array([[complex(C[i, j].mid()) for j in range(d)] for i in range(d)])
        size = np.array([np.linalg.norm(Cm[:, j]) * float(abs_upper(self.r[j])) for j in range(d)])
        Q, _ = np.linalg.qr(Cm[:, np.argsort(-size)])
        Anew = acb_mat([[acb(float(Q[i, j].real), float(Q[i, j].imag)) for j in range(d)] for i in range(d)])
        Ainv = Anew.inv()
        Jc = J * acb_mat([[z] for z in self.c])
        cnew = [acb(Jc[i, 0].real.mid(), Jc[i, 0].imag.mid()) for i in range(d)]
        pert = [Jc[i, 0] - cnew[i] + _cbox(bound * ymax / Sc[i]) + (extra[i] if extra else 0) for i in range(d)]
        v = (Ainv * C) * acb_mat([[r] for r in self.r]) + Ainv * acb_mat([[p] for p in pert])
        self.A, self.c = Anew, cnew
        self.r = [_rad_box(v[i, 0]) for i in range(d)]


def taylor_matrix(P, G, K, d):
    """Y_0..Y_K (acb_mat) with sum_i P_i s^i (sum_n n Y_n s^{n-1}) = sum_i G_i s^i sum_n Y_n s^n;
    P, G: lists over s-order of d x d acb_mat (P[0] invertible)."""
    Y = [acb_mat(d, d, [acb(int(i == j)) for i in range(d) for j in range(d)])]
    for n in range(K):
        rhs = acb_mat(d, d)
        for i in range(min(n, len(G) - 1) + 1):
            rhs += G[i] * Y[n - i]
        for i in range(1, min(n + 1, len(P) - 1) + 1):
            rhs -= P[i] * Y[n + 1 - i] * (n + 1 - i)
        Y.append(P[0].solve(rhs) / (n + 1))
    return Y


def poly_products(Pp, Gp, Y, d):
    """(D(s), Y(s), Y'(s)) as d x d acb_poly lists: D = P(s) Y'(s) - G(s) Y(s) (exact products)."""
    Yp = [[acb_poly([Y[n][r, c] for n in range(len(Y))]) for c in range(d)] for r in range(d)]
    dYp = [[acb_poly([Y[n][r, c] * n for n in range(1, len(Y))]) for c in range(d)] for r in range(d)]
    D = []
    for r in range(d):
        row = []
        for c in range(d):
            acc = acb_poly([])
            for l in range(d):
                acc += Pp[r][l] * dYp[l][c] - Gp[r][l] * Yp[l][c]
            row.append(acc)
        D.append(row)
    return D, Yp, dYp


def block_inv(M, d0=None):
    """Ball inverse; for the kappa-derivative augmentation [[A, 0], [B, A]] (d0 = size of A)
    blockwise as [[A^-1, 0], [-A^-1 B A^-1, A^-1]] (an exact zero block; interval LU would not)."""
    if d0 is None:
        return M.inv()
    d, L = M.nrows(), M.tolist()
    A = acb_mat([[L[r][c] for c in range(d0)] for r in range(d0)]).inv()
    C = -(A * acb_mat([[L[r][c] for c in range(d0)] for r in range(d0, d)]) * A)
    rows = [[A[r, c] if c < d0 else acb(0) for c in range(d)] for r in range(d0)]
    rows += [[C[r - d0, c] if c < d0 else A[r - d0, c - d0] for c in range(d)] for r in range(d0, d)]
    return acb_mat(rows)


def wnorm(rows, Sc):
    """Weighted infinity norm max_r sum_c |M_rc| Sc_r / Sc_c (upper bound)."""
    return norm_inf([[e * (Sc[r] / Sc[c]) for c, e in enumerate(row)] for r, row in enumerate(rows)])


def step_bound(sd, kap, K, d, Pp, Gp, Pbox, Gbox, incP, incG, dPt, dGt, Sc=None, boxdata=None, d0=None):
    """Y_K(-h) (acb_mat) and the Groenwall bound of the step for the kappa-combined data
    (Pp, Gp: s-polynomial matrices as acb_poly; Pbox/Gbox: lists over sub-boxes of acb_mat;
    incP/incG: arb majorant matrices; dPt/dGt: arb tails of the truncated s-polynomials), all in
    the weighted norm ``Sc``.  ``boxdata`` = (Pbox_b, Gbox_b, dP, dG) for a kappa *box* around the
    point kap: P, G over the box and P(box) - P(kap), G(box) - G(kap) per sub-box; then the bound
    also covers  Phi(kappa) - Phi(kap)  for every kappa in the box (Groenwall on
    Delta' = A(kappa) Delta + (A(kappa) - A(kap)) Phi(kap), Delta(0) = 0)."""
    Sc = Sc or [1.0] * d
    wr = max(Sc) / min(Sc)
    nrm = lambda rows: wnorm(rows, Sc)                          # noqa: E731  (weighted throughout)
    h = arb(sd.h)
    Pc = [acb_mat([[Pp[r][c][i] for c in range(d)] for r in range(d)]) for i in range(K + 1)]
    Gc = [acb_mat([[Gp[r][c][i] for c in range(d)] for r in range(d)]) for i in range(K + 1)]
    Y = taylor_matrix(Pc, Gc, K, d)
    D, Yp, dYp = poly_products(Pp, Gp, Y, d)
    Dsup = nrm([[sup_poly(D[r][c], h) for c in range(d)] for r in range(d)])
    Ysup = nrm([[sup_poly(Yp[r][c], h) for c in range(d)] for r in range(d)])
    dYsup = nrm([[sup_poly(dYp[r][c], h) for c in range(d)] for r in range(d)])
    Pinv, Lm, Lb, dA = arb(0), arb(0), arb(0), arb(0)
    for b, (Pb, Gb) in enumerate(zip(Pbox, Gbox)):
        Pi = block_inv(Pb, d0)
        Ac = Pi * Gb
        Pinv = Pinv.max(nrm(Pi.tolist()))
        Lm = Lm.max(nrm(Ac.tolist()))
        if boxdata is not None:
            Pib = block_inv(boxdata[0][b], d0)
            Lb = Lb.max(nrm((Pib * boxdata[1][b]).tolist()))
            dA = dA.max(nrm(Pib.tolist()) * (nrm(boxdata[3][b].tolist())
                                                    + nrm(boxdata[2][b].tolist()) * nrm(Ac.tolist())))
    E = Pinv * (Dsup + (dPt * wr + nrm(incP)) * dYsup + (dGt * wr + nrm(incG)) * Ysup)
    bound = abs_upper(h * (Lm * h).exp() * E)
    if boxdata is not None:
        bound = abs_upper(bound + h * (Lb * h).exp() * dA * (Ysup + bound))
    Jm = acb_mat(d, d)
    for n in reversed(range(K + 1)):
        Jm = Jm * (-h) + Y[n]
    return Jm, bound, dict(L=float(Lm), Pinv=float(Pinv), Dsup=float(Dsup), Ysup=float(Ysup), bound=float(bound))
