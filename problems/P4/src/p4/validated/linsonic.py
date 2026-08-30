"""Sonic-point series of the linear perturbation for complex kappa in a box (Theorem B, Stage 1).

System: the 3D linearised system of ``linsys``,  Pk(x; kappa) q' = Gk(x; kappa) q,
q = (N_p, W_p, V_p), along the certified A1 background series u(x) (``sonic.SonicExpansion``,
point run + tail certificate).  The sonic point is a regular singular point: Pk(0) has the
fluid null vector l of the background, and the level equations are, as in A1,
    E_n = (F^{N}_{n-1}, F^{fl}_{n-1}, l.F_n),   F = Pk theta q - Gk q,
    M_n(kappa) = n D(kappa) + E(kappa),   D = (kappa - A_0) S_0 D^,   E = -(0, 0, l.Gk_0).
Order 0: N_p(0) = 0 (gauge), A_p(0) = 1 (normalisation, via the linearised constraint
A_0 (C~_W W_p0 + C~_V V_p0) = (kappa - A_0) S_0) and l.Gk_0 q_0 = 0 (row proportionality):
one free amplitude, as in S1/KHA99.  Orders n >= 1: one ball solve each (M_n provably
invertible = non-resonance certificate for every kappa in the box).

kappa-dependence: Taylor model in delta = kappa - kappa_c (kappa_c real):
    q_n(kappa) = sum_{k<=m} q_{n,k} delta^k + R_n,   |R_n| <= sup_box |q_{n,m+1}(xi)| |delta|^{m+1},
polynomial part from the recursion at the exact centre (delta a formal variable; the level
solve is triangular in the delta-degree, so every retained coefficient is exact), the
(m+1)-th coefficient from the same recursion with the whole box as base point (integral form
of Taylor's remainder along the segment [kappa_c, kappa], which lies in the convex box).
Naive box propagation would blow up ~5x per order (A1 note sec. 3).

Tail (``certify_linear``): affine contraction on the tail in l^1_nu, nu < nu_u of the background
certificate, with the background tail entering through Banach-algebra increments of the
polynomial coefficient matrices (see notes/s2-theorem-b.md sec. 1).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from flint import acb, acb_mat, arb

from . import recursion
from .arbseries import Series, l1nu, precision, to_arb
from .linsys import LinSystem, abs_up, kappa_box, to_acb
from .polysys import _PolyEvaluator, abs_eval, abs_eval_increment
from .shootsys import eval_box
from .tailbound import Certificate

NU_CANDIDATES = ("0.09", "0.08", "0.065", "0.05")


def _col(v):
    return acb_mat([[to_acb(x)] for x in v])


def _zero_col():
    return acb_mat([[acb(0)], [acb(0)], [acb(0)]])


def norm_inf(M):
    """Max row sum of exact upper bounds |entries| of an acb_mat (arb)."""
    best = arb(0)
    for i in range(M.nrows()):
        s = arb(0)
        for j in range(M.ncols()):
            s += abs_up(M[i, j])
        best = best.max(s)
    return best


def _vec_norm(v):
    return max((abs_up(v[i, 0]) for i in range(v.nrows())), key=lambda a: a.mid())


class _Structure:
    """acb coefficient matrices of Pk, Gk along the truncated background for kappa = base + delta:
    P[dg][k], G[dg][k] = coefficient of x^k delta^dg (dg <= 1 for Pk, <= 2 for Gk)."""

    def __init__(self, Pser, Gser, base):
        k = to_acb(base)
        lp = max(len(s) for M in Pser for row in M for s in row)
        lg = max(len(s) for M in Gser for row in M for s in row)

        def mats(ser, pw, length):
            out = []
            for n in range(length):
                out.append(acb_mat([[sum((pw[j] * ser[j][r][q][n] for j in range(len(ser))), acb(0))
                                     for q in range(3)] for r in range(3)]))
            return out
        one, zero = acb(1), acb(0)
        self.P = [mats(Pser, [one, k], lp), mats(Pser, [zero, one], lp)]
        self.G = [mats(Gser, [one, k, k * k], lg), mats(Gser, [zero, one, 2 * k], lg), mats(Gser, [zero, zero, one], lg)]
        self.lp, self.lg = lp, lg

    def DE(self, rows):
        """(D, E) as lists over the delta-degree of acb_mat (D: 0..1, E: 0..2)."""
        D = [[[acb(0)] * 3 for _ in range(3)] for _ in range(2)]
        E = [[[acb(0)] * 3 for _ in range(3)] for _ in range(3)]
        for i, row in enumerate(rows):
            for c, r, s in row:
                for dg in range(2):
                    for l in range(3):
                        D[dg][i][l] += c * self.P[dg][s + 1][r, l]
                if s >= 0:
                    for dg in range(3):
                        for l in range(3):
                            E[dg][i][l] -= c * self.G[dg][s][r, l]
        return [acb_mat(d) for d in D], [acb_mat(e) for e in E]

    def residual(self, q, j, upto, dm):
        """delta-polynomial (list over dg < dm of 3x1) of F_j = [Pk theta q - Gk q]_j using q_m, m <= upto."""
        out = [_zero_col() for _ in range(dm)]
        for i in range(j + 1):
            if i + 1 <= upto and j - i < self.lp:
                for a in range(2):
                    for b in range(dm - a):
                        out[a + b] += self.P[a][j - i] * q[i + 1][b] * (i + 1)
            if i <= upto and j - i < self.lg:
                for a in range(3):
                    for b in range(dm - a):
                        out[a + b] -= self.G[a][j - i] * q[i][b]
        return out

    def level_values(self, q, n, rows, upto, dm):
        """E_n(q) as a delta-polynomial (list over dg of 3x1), q_m used for m <= upto only."""
        F = {}
        for row in rows:
            for c, r, s in row:
                if n + s >= 0 and n + s not in F:
                    F[n + s] = self.residual(q, n + s, upto, dm)
        out = [_zero_col() for _ in range(dm)]
        for i, row in enumerate(rows):
            for c, r, s in row:
                if n + s >= 0:
                    for dg in range(dm):
                        out[dg][i, 0] += c * F[n + s][dg][r, 0]
        return out


# ----------------------------------------------------------------------------
# the recursion (one base point, delta a formal variable truncated at dm = m + 2)
# ----------------------------------------------------------------------------
def _run(L, st, u0, rows, base, K, dm):
    """Coefficients q_n (list over n of lists over dg < dm of 3x1 acb_mat) for kappa = base + delta."""
    args = [arb(0)] + list(u0)
    A0, S0 = u0[0], eval_box(L.S, args)
    CW, CV = eval_box(L.dC[2], args), eval_box(L.dC[3], args)
    ell = [c for (c, _, _) in rows[2]]
    D, E = st.DE(rows)
    # order 0: l.Gk_0 q_0 = 0 (columns W, V; N_p0 = 0) and A_0 (C~_W W + C~_V V) = (kappa - A_0) S_0
    G0 = st.G
    M0 = [acb_mat([[ell[0] * G0[dg][0][1, l] + ell[1] * G0[dg][0][2, l] for l in (1, 2)],
                   [A0 * CW if dg == 0 else acb(0), A0 * CV if dg == 0 else acb(0)]]) for dg in range(3)]
    rhs = [acb_mat([[acb(0)], [(to_acb(base) - A0) * S0]]), acb_mat([[acb(0)], [to_acb(S0)]])]
    wv = []
    for k in range(dm):
        b = rhs[k] if k < 2 else acb_mat([[acb(0)], [acb(0)]])
        for j in range(1, min(k, 2) + 1):
            b = b - M0[j] * wv[k - j]
        wv.append(M0[0].solve(b))
    q = [[acb_mat([[acb(0)], [wv[k][0, 0]], [wv[k][1, 0]]]) for k in range(dm)]]
    for n in range(1, K + 1):
        R = st.level_values(q, n, rows, n - 1, dm)
        M = [D[0] * n + E[0], D[1] * n + E[1], E[2]]
        qn = []
        for k in range(dm):
            b = -R[k]
            for j in range(1, min(k, 2) + 1):
                b = b - M[j] * qn[k - j]
            qn.append(M[0].solve(b))
        q.append(qn)
    return q, D, E


@dataclass
class LinSonicExpansion:
    kappa_c: acb                    # centre (real)
    width: float                    # box half-width (both parts); kappa = kappa_c + [-w,w] + i[-w,w]
    m: int
    K: int
    L: LinSystem
    bg: object                      # certified SonicExpansion (point run)
    rows: list
    point: dict                     # recursion at the centre (delta-polynomials)
    box: dict = None                # recursion with the box as base point
    balls: list = None              # list over n of [N_p, W_p, V_p] acb enclosures over the box
    Ap: list = None                 # list over n of A_p_n enclosures (constraint, series division)
    cert: object = None
    info: dict = field(default_factory=dict)

    @property
    def kappa(self):
        return kappa_box(self.kappa_c, self.width) if self.width else self.kappa_c

    def floats(self):
        """(4, K+1) complex midpoints of (A_p, N_p, W_p, V_p)."""
        rows = [self.Ap] + [[b[i] for b in self.balls] for i in range(3)]
        return np.array([[complex(float(z.real.mid()), float(z.imag.mid())) for z in r] for r in rows])

    def radii(self):
        rows = [self.Ap] + [[b[i] for b in self.balls] for i in range(3)]
        return np.array([[float(abs_up(z - z.mid())) for z in r] for r in rows])

    def certify(self, nu=None, prec=256):
        with precision(prec):
            self.cert = certify_linear(self, nu=nu)
        return self.cert

    def eval(self, x, with_tail=True):
        """Balls (A_p, N_p, W_p, V_p)(x) over the box (tail included if certified, |x| <= nu)."""
        x = to_arb(x)
        q = []
        for i in range(3):
            acc = acb(0)
            for n in reversed(range(self.K + 1)):
                acc = acc * x + self.balls[n][i]
            q.append(acc)
        if with_tail:
            if self.cert is None:
                raise ValueError("not certified; call certify() or use with_tail=False")
            tb = self.cert.tail_bound(abs(x).abs_upper())
            q = [v + acb(arb(0, tb), arb(0, tb)) for v in q]
        u = self.bg.eval(x, with_tail=with_tail)
        return [self.L.A_p(u, self.kappa, q)] + q


# ----------------------------------------------------------------------------
# driver, A_p series (constraint), row-1 consistency check
# ----------------------------------------------------------------------------
def _dmul(a, b, dm):
    """Product of delta-polynomials (lists of acb) truncated at degree dm - 1."""
    out = [acb(0)] * dm
    for i, ai in enumerate(a):
        for j in range(min(len(b), dm - i)):
            out[i + j] += ai * b[j]
    return out


def _Ap_delta(L, ubar, base, q, K, dm):
    """A_p_n(delta), n <= K, of A (C~_N N_p + C~_W W_p + C~_V V_p) / ((kappa - A) S), kappa = base + delta,
    by series division in x with delta-polynomial coefficients (division by the delta-polynomial
    den_0(delta) = (base - A_0) S_0 + delta S_0 through its truncated inverse)."""
    ev = _PolyEvaluator([Series.var(K + 1)] + [s.with_cap(K + 1) for s in ubar])
    Cq = [ev(L.dC[l + 1]) for l in range(3)]
    A, S, AS = ubar[0], ev(L.S), ev(L.Apoly * L.S)
    kap = to_acb(base)
    zero = [acb(0)] * dm

    def cmul(c, v):                                  # real scalar times delta-polynomial
        return [to_acb(c) * x for x in v]

    def add(u, v):
        return [x + y for x, y in zip(u, v)]

    num = []
    for n in range(K + 1):
        acc = list(zero)
        for l in range(3):
            for j in range(n + 1):
                acc = add(acc, cmul(Cq[l][n - j], [q[j][dg][l, 0] for dg in range(dm)]))
        num.append(acc)
    num2 = []
    for n in range(K + 1):
        acc = list(zero)
        for j in range(n + 1):
            acc = add(acc, cmul(A[n - j], num[j]))
        num2.append(acc)
    den = [[kap * S[n] - AS[n], to_acb(S[n])] + [acb(0)] * (dm - 2) for n in range(K + 1)]
    d0 = den[0]
    inv0 = [1 / d0[0]]                               # 1/(d00 + delta d01) truncated
    for k in range(1, dm):
        inv0.append(-inv0[k - 1] * d0[1] / d0[0])
    a = []
    for n in range(K + 1):
        s = num2[n]
        for j in range(1, n + 1):
            s = [x - y for x, y in zip(s, _dmul(den[j], a[n - j], dm))]
        a.append(_dmul(s, inv0, dm))
    return a


def row1_residual(ex):
    """rho_n, n < K, of  S A_p' - 2 V A' V_p - DQ_0 . p  (row 1 of the 4D linearisation, not used by
    the 3D system): all contain 0 iff the constraint-eliminated series solves the full linearised
    system to order K (Bianchi identity)."""
    L, K = ex.L, ex.K
    ubar = ex.bg.series()
    ev = _PolyEvaluator([Series.var(K + 1)] + [s.with_cap(K + 1) for s in ubar])
    S, VA = ev(L.S), ubar[3] * ubar[0].deriv()
    dQ0 = [ev(L.sys.dQ()[0][l]) for l in range(4)]
    p = [ex.Ap] + [[b[i] for b in ex.balls] for i in range(3)]
    Apd = [ex.Ap[n + 1] * (n + 1) for n in range(K)]
    out = []
    for n in range(K):
        s = sum((to_acb(S[n - j]) * Apd[j] for j in range(n + 1)), acb(0))
        s -= 2 * sum((to_acb(VA[n - j]) * p[3][j] for j in range(n + 1)), acb(0))
        for l in range(4):
            s -= sum((to_acb(dQ0[l][n - j]) * p[l][j] for j in range(n + 1)), acb(0))
        out.append(s)
    return out


def linear_sonic_expansion(bg, kappa_c, width=0.0, m=5, K=40, prec=256, check=True):
    """Perturbation series at the sonic point for kappa in kappa_c + [-w, w] + i[-w, w] (w = ``width``;
    w = 0: point).  ``bg``: certified A1 point expansion with bg.K >= K + 1 (q_K needs u_{K+1})."""
    if bg.K < K + 1:
        raise ValueError("background expansion must have K_bg >= K + 1")
    with precision(prec):
        L = LinSystem()
        kc = to_acb(kappa_c)
        ubar = bg.series()
        Pser, Gser = L.series_parts(ubar)
        ell, row = bg.point["ell"], bg.point["row"]
        rows = [[(arb(1), 0, -1)], [(arb(1), row - 1, -1)], [(ell[0], 1, 0), (ell[1], 2, 0)]]
        u0, dm = bg.balls[0], m + 2
        st_pt = _Structure(Pser, Gser, kc)
        q_pt, D_pt, E_pt = _run(L, st_pt, u0, rows, kc, K, dm)
        a_pt = _Ap_delta(L, ubar, kc, q_pt, K, dm)
        point = dict(coefs=q_pt, Ap=a_pt, D=D_pt, E=E_pt, structure=st_pt)
        box, rem = None, None
        if width:
            kb = kappa_box(kc, width)
            st_bx = _Structure(Pser, Gser, kb)
            q_bx, D_bx, E_bx = _run(L, st_bx, u0, rows, kb, K, dm)
            a_bx = _Ap_delta(L, ubar, kb, q_bx, K, dm)
            box = dict(coefs=q_bx, Ap=a_bx, D=D_bx, E=E_bx, structure=st_bx)
            dw = acb(arb(0, width), arb(0, width))
            wpow = abs_up(dw) ** (m + 1)

            def model(poly, top):                    # Taylor model in delta -> ball over the box
                top = abs_up(top) * wpow
                acc = acb(0)
                for k in reversed(range(m + 1)):
                    acc = acc * dw + poly[k]
                return acc + acb(arb(0, top), arb(0, top)), float(top)
            balls, Ap, rem = [], [], []
            for n in range(K + 1):
                bn, rn = [], []
                for i in range(3):
                    b, r = model([q_pt[n][k][i, 0] for k in range(dm)], q_bx[n][m + 1][i, 0])
                    bn.append(b)
                    rn.append(r)
                a, r = model(a_pt[n], a_bx[n][m + 1])
                balls.append(bn)
                Ap.append(a)
                rem.append([r] + rn)
        else:
            balls = [[q_pt[n][0][i, 0] for i in range(3)] for n in range(K + 1)]
            Ap = [a_pt[n][0] for n in range(K + 1)]
        ex = LinSonicExpansion(kc, float(width), m, K, L, bg, rows, point, box, balls, Ap, info=dict(rem=rem))
        if check:
            rho = row1_residual(ex)
            ex.info["row1_contains_zero"] = all(r.contains(acb(0)) for r in rho)
            ex.info["row1_max_abs"] = max(float(abs_up(r)) for r in rho)
        return ex


# ----------------------------------------------------------------------------
# tail certificate: affine contraction on the tail of the linear recursion in l^1_nu
# ----------------------------------------------------------------------------
def _arb_norm_inf(rows):
    best = arb(0)
    for row in rows:
        best = best.max(sum(row, arb(0)))
    return best


def certify_linear(ex, nu=None):
    """Certificate for the tail of ``ex.balls`` (n > K), valid for every kappa in the box.

    True coefficients q = q-bar + z.  For n > K,  M_n z_n = -[E_n(q-bar) + sum_{K<m<n} M_nm z_m]  with
    ||M_n^{-1}|| <= c/n (c = ||D^{-1}||/(1 - ||D^{-1}E||/(K+1))), so T(z) = z - M^{-1} E(q-bar + z) is affine
    with linear part of l^1_nu-norm <= Z and ||T(0)|| <= Y; if Z < 1 the unique fixed point (the true
    tail, by uniqueness of the recursion) has ||z||_nu <= Y/(1-Z) =: eps.  The coefficient matrices
    along the *true* background u = u-bar + v enter as (polynomial part along u-bar: finite exact sums,
    Z1/Y1) + (increments over the background tail, |v|_nu <= eps_v, |theta v|_nu <= eps_th from the A1
    certificate at nu < nu_u, via Banach-algebra majorants: Z2/Y2)."""
    L, K, rows = ex.L, ex.K, ex.rows
    run = ex.box or ex.point
    st, D, E = run["structure"], run["D"][0], run["E"][0]
    cert_u = ex.bg.cert
    if cert_u is None or not cert_u.ok:
        raise ValueError("background expansion must be certified first")
    Dinv = D.inv()
    g, Dn = norm_inf(Dinv * E), norm_inf(Dinv)
    if not g < K + 1:
        return Certificate(K, arb(0), arb(0), arb(0), arb(0), arb(0), arb(0), False, dict(reason="||D^-1 E|| >= K+1"))
    c = Dn / (1 - g / (K + 1))
    ubar = ex.bg.series()
    th_ubar = [s.deriv() for s in ubar]
    P0, G0 = st.P[0], st.G[0]
    kabs = abs_up(ex.kappa)
    q1 = [[_col(b)] for b in ex.balls]
    jmax = st.lp + K + 1
    F = [st.residual(q1, j, K, 1)[0] for j in range(jmax + 1)]

    def En(n):
        out = [arb(0)] * 3
        for i, row in enumerate(rows):
            for cc, r, s in row:
                if 0 <= n + s <= jmax:
                    out[i] += abs_up(cc) * abs_up(F[n + s][r, 0])
        return max(out, key=lambda a: a.mid())

    cands = [to_arb(nu)] if nu is not None else [to_arb(s) for s in NU_CANDIDATES]
    last = None
    for nu_ in cands:
        if not nu_ < cert_u.nu:
            continue
        Z1 = arb(0)
        for k in range(1, max(st.lp, st.lg)):
            Bk = [[arb(0)] * 3 for _ in range(3)]
            for i, row in enumerate(rows):
                for cc, r, s in row:
                    for l in range(3):
                        v = abs_up(P0[k + s + 1][r, l]) if k + s + 1 < st.lp else arb(0)
                        if k + s < st.lg:
                            v += abs_up(G0[k + s][r, l]) / (K + 2)
                        Bk[i][l] += abs_up(cc) * v
            Z1 += _arb_norm_inf(Bk) * nu_**k
        Z1 *= c
        Y1 = c * sum((En(n) * nu_**n / n for n in range(K + 1, jmax + 2)), arb(0))
        eps_v, eps_th = cert_u.tail_bound(nu_), cert_u.deriv_tail_bound(nu_)
        ru = [l1nu(s, nu_) for s in ubar]
        thu = [l1nu(s, nu_) for s in th_ubar]
        args = [nu_] + ru
        args_e = [nu_] + [r + eps_v for r in ru]
        incP = [[sum((kabs**j * abs_eval_increment(L.Pk[j][r][q], args, eps_v) for j in range(2)), arb(0))
                 for q in range(3)] for r in range(3)]
        incG = [[arb(0)] * 3 for _ in range(3)]
        for r in range(3):
            for q in range(3):
                for j in range(3):
                    term = abs_eval_increment(L.Ga[j][r][q], args, eps_v)
                    for i in (range(4) if j < 2 else ()):
                        if not L.Gb[j][r][q][i].is_zero():
                            term += abs_eval_increment(L.Gb[j][r][q][i], args, eps_v) * thu[i]
                            term += abs_eval(L.Gb[j][r][q][i], args_e) * eps_th
                    incG[r][q] += kabs**j * term
        Zm = [[arb(0)] * 3 for _ in range(3)]
        Y2 = arb(0)
        nq = [sum((abs_up(ex.balls[n][l]) * nu_**n for n in range(K + 1)), arb(0)) for l in range(3)]
        nth = [sum((abs_up(ex.balls[n][l]) * n * nu_ ** (n - 1) for n in range(1, K + 1)), arb(0)) for l in range(3)]
        for i, row in enumerate(rows):
            for cc, r, s in row:
                ac = abs_up(cc)
                for l in range(3):
                    Zm[i][l] += ac * (nu_ ** (-(s + 1)) * incP[r][l] + nu_ ** (-s) * incG[r][l] / (K + 2))
                    Y2 += ac * nu_ ** (-s) * (incP[r][l] * nth[l] + incG[r][l] * nq[l])
        Z2 = c * _arb_norm_inf(Zm)
        Y2 = c * Y2 / (K + 1)
        Z, Y = Z1 + Z2, Y1 + Y2
        det = dict(g=g, Dinv_norm=Dn, eps_v=eps_v, eps_th=eps_th, Y1=Y1, Y2=Y2)
        if Z < 1:
            return Certificate(K, nu_, (Y / (1 - Z)).abs_upper(), Y, Z1, Z2, c, True, det)
        last = Certificate(K, nu_, arb(0), Y, Z1, Z2, c, False, dict(reason="Z >= 1", **det))
    return last


def gauge_eigenvalue(bg):
    """kappa-bar = 2 - A_0 + 2 W_0 / 3 (pure-gauge mode of the sonic-point gauge) as an arb ball."""
    A0, W0 = bg.balls[0][0], bg.balls[0][2]
    return 2 - A0 + W0 * 2 / 3
