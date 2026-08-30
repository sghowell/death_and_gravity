"""Theorem B, Stage 3: validated propagation of a *Taylor model in kappa* of the 4D linearised
solution along the certified background tube (``lintube.Tube``, system ``linscaled.full_system``).

For kappa = kappa_c + delta, |delta| <= r, the scaled 4D system  P(s) y' = (G_c(s) + delta G_1(s)) y
(P kappa-free, G affine in kappa: G_c = G_0 + kappa_c G_1) is propagated as
    y(x; delta) = sum_{k<=m} y_k(x) delta^k + R(x; delta),
with the polynomial part (y_0, ..., y_m) -- the kappa-Taylor coefficients y_k = d^k y/dkappa^k / k! --
carried as m+1 complex Lohner sets in C^4 (disc radii) and the remainder R as one more Lohner set.  Per step:
 1. block-Toeplitz Taylor recursion (the order-m kappa-derivative system is lower block-triangular
    Toeplitz with A_c on the diagonal and A_1 below; its fundamental matrix has blocks Y_k(s) =
    d^k Y/dkappa^k / k!):  (n+1) P_0 Y_k^{(n+1)} = sum_i G_{c,i} Y_k^{(n-i)} + sum_i G_{1,i} Y_{k-1}^{(n-i)}
    - sum_{i>=1} (n+1-i) P_i Y_k^{(n+1-i)},  Y_0(0) = I, Y_k(0) = 0 (k >= 1);
 2. Groenwall bound of the block system in the weighted norm ||y^||_lam = max_k lam^k |y_k|_inf
    (lintube/linstep argument: defect D_k = P Y_k' - G_c Y_k - G_1 Y_{k-1} along the true background,
    ||A^||_lam <= L + lam L_1, L = sup||P^-1 G_c||, L_1 = sup||P^-1 G_1|| over the tube sub-boxes):
    ||Y^(-h) - Y^_K(-h)||_lam <= bound, so block k of Y^ y^ has error <= bound ||y^||_lam / lam^k;
 3. Cauchy remainder of the one-step map Y(-h; kappa) (entire in kappa): sup_{|delta|<=rho} ||Y(-h; kappa)||
    <= e^{h (L + rho L_1)} =: M, hence  ||Y(-h; delta) - sum_{k<=m} Y_k(-h) delta^k|| <= M (r/rho)^{m+1}
    / (1 - r/rho) =: eps_Y for |delta| <= r, with rho = (m+1)/(h L_1) (the optimum);
 4. composition: (sum_k Y_k delta^k + R_Y)(sum_l y_l delta^l + R) -> block-Toeplitz map of the
    polynomial part (Lohner sets, sets k-j entering block k as balls), the degrees m+1..2m and
    R_Y y go to the thin set S_0 (delta-independent).
Everything is an arb upper bound; the kappa-independent step data are ``Tube.system_data``.
"""
from __future__ import annotations

import time

import numpy as np
from flint import acb, acb_mat, acb_poly, arb

from .arbseries import abs_upper, precision
from .linstep import norm_inf, sup_poly
from .linsys import to_acb


def _cbox(rad):
    return acb(arb(0, rad), arb(0, rad))


def _box_of(z):
    """acb box centred at 0 containing the ball z (midpoint included in the radius: safe)."""
    return acb(arb(0, abs_upper(z.real)), arb(0, abs_upper(z.imag)))


def _col(v):
    return acb_mat([[x] for x in v])


def _mat(rows):
    return acb_mat([[acb(e) for e in row] for row in rows])


class LSet:
    """Complex Lohner set {c + A r : |r_j| <= rho_j} in C^d with *disc* radii (rotation-invariant: an
    axis-aligned complex box would lose up to sqrt(2) per step under the rotating flows at large Im kappa),
    exact-midpoint update, extra balls accepted."""

    def __init__(self, balls):
        self.d = len(balls)
        self.c = [acb(z.real.mid(), z.imag.mid()) for z in balls]
        self.A = acb_mat(self.d, self.d, [acb(int(i == j)) for i in range(self.d) for j in range(self.d)])
        self.rho = [abs_upper(z - self.c[i]) for i, z in enumerate(balls)]

    def radii(self):
        """Disc radius of every component of the hull: R_i = sum_j |A_ij| rho_j."""
        return [sum((abs_upper(self.A[i, j]) * self.rho[j] for j in range(self.d)), arb(0)) for i in range(self.d)]

    def hull(self):
        return [self.c[i] + _cbox(R) for i, R in enumerate(self.radii())]

    def max_abs(self):
        return max((abs_upper(self.c[i]) + R for i, R in enumerate(self.radii())), default=arb(0))

    def propagate(self, J, extra=None, rad=None):
        """New set containing  J (c + A r) + extra_i + disc(rad_i)  (extra: acb balls, rad: arb moduli)."""
        d = self.d
        Jc = J * _col(self.c)
        vc = [Jc[i, 0] + (extra[i] if extra is not None else 0) for i in range(d)]
        cnew = [acb(v.real.mid(), v.imag.mid()) for v in vc]
        pert = [abs_upper(vc[i] - cnew[i]) + (rad[i] if rad is not None else 0) for i in range(d)]
        C = J * self.A
        Cm = np.array([[complex(C[i, j].mid()) for j in range(d)] for i in range(d)])
        size = np.array([np.linalg.norm(Cm[:, j]) * float(self.rho[j]) for j in range(d)])
        Q, _ = np.linalg.qr(Cm[:, np.argsort(-size)])
        Anew = acb_mat([[acb(float(Q[i, j].real), float(Q[i, j].imag)) for j in range(d)] for i in range(d)])
        Ainv = Anew.inv()
        M = Ainv * C
        self.rho = [abs_upper(sum((abs_upper(M[i, j]) * self.rho[j] + abs_upper(Ainv[i, j]) * pert[j] for j in range(d)), arb(0)))
                    for i in range(d)]
        self.A, self.c = Anew, cnew


def combine4(data, kc, K):
    """kappa_c-combined step data of the 4D system: s-order matrices and polynomials of P, G_c, G_1,
    their sub-box values, increment majorants and truncation tails."""
    kc = to_acb(kc)
    ka = abs_upper(kc)
    (Pco, Pt), (G0co, G0t), (G1co, G1t) = data["P"][0], data["G"][0], data["G"][1]
    d = len(Pco)
    P = [_mat([[Pco[r][c][i] for c in range(d)] for r in range(d)]) for i in range(K + 1)]
    G1 = [_mat([[G1co[r][c][i] for c in range(d)] for r in range(d)]) for i in range(K + 1)]
    Gc = [_mat([[G0co[r][c][i] + kc * G1co[r][c][i] for c in range(d)] for r in range(d)]) for i in range(K + 1)]
    poly = lambda M: [[acb_poly([M[i][r, c] for i in range(K + 1)]) for c in range(d)] for r in range(d)]   # noqa: E731
    nb = len(data["Pbox"][0])
    Pbox = [_mat(data["Pbox"][0][b]) for b in range(nb)]
    G1box = [_mat(data["Gbox"][1][b]) for b in range(nb)]
    Gcbox = [_mat([[data["Gbox"][0][b][r][c] + kc * data["Gbox"][1][b][r][c] for c in range(d)] for r in range(d)])
             for b in range(nb)]
    incP = norm_inf(data["incP"][0])
    incG1 = norm_inf(data["incG"][1])
    incGc = norm_inf([[data["incG"][0][r][c] + ka * data["incG"][1][r][c] for c in range(d)] for r in range(d)])
    return dict(P=P, Gc=Gc, G1=G1, Pp=poly(P), Gcp=poly(Gc), G1p=poly(G1), Pbox=Pbox, Gcbox=Gcbox, G1box=G1box,
                incP=incP, incGc=incGc, incG1=incG1, dPt=Pt, dGct=G0t + ka * G1t, dG1t=G1t, d=d)


def block_taylor(cd, K, m):
    """Y[k][n]: s-order-n coefficient of the kappa-Taylor block Y_k(s), k <= m, n <= K."""
    P, Gc, G1, d = cd["P"], cd["Gc"], cd["G1"], cd["d"]
    Id = acb_mat(d, d, [acb(int(i == j)) for i in range(d) for j in range(d)])
    Y = [[Id if k == 0 else acb_mat(d, d)] for k in range(m + 1)]
    for n in range(K):
        for k in range(m + 1):
            rhs = acb_mat(d, d)
            for i in range(n + 1):
                rhs += Gc[i] * Y[k][n - i]
                if k:
                    rhs += G1[i] * Y[k - 1][n - i]
            for i in range(1, n + 2):
                rhs -= P[i] * Y[k][n + 1 - i] * (n + 1 - i)
            Y[k].append(P[0].solve(rhs) / (n + 1))
    return Y


def tm_step(sd, cd, K, m, lam, r):
    """One step: (Y_k(-h) for k <= m, bound, eps_Y, info) for the block system of order m."""
    h, d = arb(sd.h), cd["d"]
    lam, r = arb(lam), arb(r)
    Y = block_taylor(cd, K, m)
    Pp, Gcp, G1p = cd["Pp"], cd["Gcp"], cd["G1p"]
    Ysup, dYsup, Dsup, Yprev = [], [], [], None
    for k in range(m + 1):
        Yp = [[acb_poly([Y[k][n][rr, c] for n in range(K + 1)]) for c in range(d)] for rr in range(d)]
        dYp = [[acb_poly([Y[k][n][rr, c] * n for n in range(1, K + 1)]) for c in range(d)] for rr in range(d)]
        D = []
        for rr in range(d):
            row = []
            for c in range(d):
                acc = acb_poly([])
                for l in range(d):
                    acc += Pp[rr][l] * dYp[l][c] - Gcp[rr][l] * Yp[l][c]
                    if k:
                        acc -= G1p[rr][l] * Yprev[l][c]
                row.append(sup_poly(acc, h))
            D.append(row)
        Dsup.append(norm_inf(D))
        Ysup.append(norm_inf([[sup_poly(Yp[rr][c], h) for c in range(d)] for rr in range(d)]))
        dYsup.append(norm_inf([[sup_poly(dYp[rr][c], h) for c in range(d)] for rr in range(d)]))
        Yprev = Yp
    Pinv, L, L1 = arb(0), arb(0), arb(0)
    for Pb, Gb, G1b in zip(cd["Pbox"], cd["Gcbox"], cd["G1box"]):
        Pi = Pb.inv()
        Pinv = Pinv.max(norm_inf(Pi.tolist()))
        L = L.max(norm_inf((Pi * Gb).tolist()))
        L1 = L1.max(norm_inf((Pi * G1b).tolist()))
    lp = [lam**k for k in range(m + 1)]
    wsum = lambda v: sum((lp[k] * v[k] for k in range(m + 1)), arb(0))          # noqa: E731
    E = Pinv * (wsum(Dsup) + (cd["dPt"] + cd["incP"]) * wsum(dYsup)
                + ((cd["dGct"] + cd["incGc"]) + lam * (cd["dG1t"] + cd["incG1"])) * wsum(Ysup))
    bound = abs_upper(h * ((L + lam * L1) * h).exp() * E)
    rho = arb(m + 1) / (h * L1)
    if not rho > 2 * r:
        rho = 2 * r
    q = r / rho
    epsY = abs_upper((h * (L + rho * L1)).exp() * q ** (m + 1) / (1 - q))
    Yend = []
    for k in range(m + 1):
        acc = acb_mat(d, d)
        for n in reversed(range(K + 1)):
            acc = acc * (-h) + Y[k][n]
        Yend.append(acc)
    return Yend, bound, epsY, dict(L=float(L), L1=float(L1), Pinv=float(Pinv), Dsup=float(Dsup[0]),
                                    bound=float(bound), epsY=float(epsY), h=sd.h, x=float(sd.x) - sd.h)


class TMState:
    """Polynomial part of the kappa-Taylor model: m+1 Lohner sets S_k in C^d with
    y(x; delta) in S_0 + sum_{k>=1} S_k delta^k for every |delta| <= r.  Every delta-independent error
    (Cauchy remainder of the step map, degrees > m of the composition, Groenwall) is added to the *thin*
    set S_0 -- an isotropic remainder set of its own would wrap much faster under the rotating flow."""

    def __init__(self, coefs):
        self.m = len(coefs) - 1
        self.sets = [LSet(c) for c in coefs]

    def hulls(self):
        return [s.hull() for s in self.sets]

    def step(self, Yend, bound, epsY, lam, r):
        m, d = self.m, self.sets[0].d
        hulls = self.hulls()
        absh = [s.max_abs() for s in self.sets]
        lam, r = arb(lam), arb(r)
        lp = [lam**k for k in range(m + 1)]
        rp = [r**k for k in range(2 * m + 1)]
        ymax = max((lp[j] * absh[j] for j in range(m + 1)), key=lambda a: a.mid())
        Yn = [norm_inf(Yend[k].tolist()) + bound / lp[k] for k in range(m + 1)]
        extras, rads = [], []
        for k in range(m + 1):
            ex = [acb(0)] * d
            for j in range(1, k + 1):
                v = Yend[j] * _col(hulls[k - j])
                ex = [e + v[i, 0] for i, e in enumerate(ex)]
            extras.append(ex)
            rads.append(abs_upper(bound * ymax / lp[k]))
        T = arb(0)
        for j in range(m + 1, 2 * m + 1):
            T += rp[j] * sum((Yn[k] * absh[j - k] for k in range(j - m, m + 1)), arb(0))
        S = sum((absh[l] * rp[l] for l in range(m + 1)), arb(0))
        rads[0] = abs_upper(rads[0] + epsY * S + T)
        for k in range(m + 1):
            self.sets[k].propagate(Yend[0], extras[k], [rads[k]] * d)

    def width(self):
        return max(float(R) for s in self.sets for R in s.radii())


def propagate_tm(tube, S, kappa_c, r, m, x_stop=None, K=None, prec=256, nsub=8, lam=None, verbose=False):
    """Taylor model of the fundamental matrix Phi(x_end, x_start; kappa_c + delta), |delta| <= r:
    the d unit vectors are propagated as TMStates with exact initial data (no wrapping of an initial
    remainder).  Returns (Phi = [Phi_k acb_mat, k <= m] with all remainders inside Phi_0's balls, None
    (kept for the signature), log); Phi(delta) y in sum_k Phi_k delta^k y for |delta| <= r."""
    kc = to_acb(kappa_c)
    lam = 1.0 if lam is None else lam            # block weights lam^k: lam = 1 (the kappa-scale of the solution) beats lam = r
    d = S.d
    states = [TMState([[acb(int(i == j)) for i in range(d)]] + [[acb(0)] * d] * m) for j in range(d)]
    log, t0 = [], time.time()
    with precision(prec):
        for sd in tube.steps:
            if x_stop is not None and float(sd.x) <= x_stop + 1e-15:
                break
            data = tube.system_data(sd, S, nsub)
            Kl = len(sd.co) - 1 if K is None else min(K, len(sd.co) - 1)
            cd = combine4(data, kc, Kl)
            Yend, bound, epsY, info = tm_step(sd, cd, Kl, m, lam, r)
            for st in states:
                st.step(Yend, bound, epsY, lam, r)
            info.update(width=max(st.width() for st in states), rem=0.0)
            log.append(info)
            if verbose:
                print(f"x={info['x']:+.4f} h={info['h']:.4f} L={info['L']:.3g} L1={info['L1']:.3g} bound={info['bound']:.1e} "
                      f"epsY={info['epsY']:.1e} width={info['width']:.1e} rem={info['rem']:.1e}", flush=True)
        hulls = [st.hulls() for st in states]
        Phi = [acb_mat([[hulls[j][k][i] for j in range(d)] for i in range(d)]) for k in range(m + 1)]
    log.append(dict(time=time.time() - t0))
    return Phi, None, log


def tm_apply(Phi, Rcols, coefs, rem, r):
    """Compose: (sum_k Phi_k delta^k + R_Phi)(sum_l c_l delta^l + rem) truncated at degree m = len(Phi)-1.
    coefs: m+1 lists of d balls; rem: d arb radii (|delta| <= r).  Returns (coefs', rem')."""
    m, d = len(Phi) - 1, Phi[0].nrows()
    r = arb(r)
    rp = [r**k for k in range(2 * m + 1)]
    prod = {}
    for k in range(m + 1):
        for l in range(m + 1):
            v = Phi[k] * _col(coefs[l])
            prod[(k, l)] = [v[i, 0] for i in range(d)]
    out = [[sum((prod[(k, j - k)][i] for k in range(j + 1)), acb(0)) for i in range(d)] for j in range(m + 1)]
    trunc = [arb(0)] * d
    for j in range(m + 1, 2 * m + 1):
        for i in range(d):
            trunc[i] += rp[j] * abs_upper(sum((prod[(k, j - k)][i] for k in range(j - m, m + 1)), acb(0)))
    absPhi = [[sum((abs_upper(Phi[k][i, c]) * rp[k] for k in range(m + 1)), arb(0)) for c in range(d)] for i in range(d)]
    B = [sum((abs_upper(coefs[l][c]) * rp[l] for l in range(m + 1)), arb(0)) + rem[c] for c in range(d)]
    rem_out = []
    for i in range(d):
        acc = trunc[i]
        for c in range(d):
            acc += absPhi[i][c] * rem[c] + (B[c] * abs_upper(Rcols[c][i]) if Rcols is not None else 0)
        rem_out.append(abs_upper(acc))
    return out, rem_out


def tm_eval(coefs, rem, delta):
    """Balls of the model at delta (acb, possibly a box): sum_k c_k delta^k + box(rem)."""
    out = []
    for i in range(len(rem)):
        acc = acb(0)
        for k in reversed(range(len(coefs))):
            acc = acc * delta + coefs[k][i]
        out.append(acc + _cbox(rem[i]))
    return out
