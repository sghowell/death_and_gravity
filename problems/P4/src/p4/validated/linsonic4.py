"""Theorem B, Stage 3: sonic-point series of the *4D* linearised system (A_p kept) for complex
kappa -- point, box, or delta-Taylor model -- with the tail certificate.

System (``linsys.LinSystem`` data, plain variables u = (A, N, W, V), x around the sonic point):
    P(u) p' = [DQ(u) - Psi(u, u') - kappa P_s(u)] p =: G(x; kappa) p,   P kappa-free, G affine in kappa.
P(u_0) has rank 3 (fluid null vector ell); the level equations are exactly A1's
    E_n = (F^A_{n-1}, F^N_{n-1}, F^{row}_{n-1}, ell . F_n),   F = P p' - G p,
so M_n(kappa) = n D + E(kappa) with D kappa-*independent* (block-triangular, invertible) and E of
rank <= 1 (only the ell-row), affine in kappa: M_n is singular iff n = sigma(kappa) :=
-tr(D^{-1} E(kappa)), an affine function -- the single Frobenius exponent (S1's resonances
kappa = -0.099 - 1.099 n).  Unlike the constraint-reduced 3D form (``linsonic``, D ~ (kappa - A_0),
ghost roots accumulating at A_0) this recursion is regular on the whole rectangle R.
Order 0: A_p(0) = 1 (normalisation), N_p(0) = 0 (gauge), ell . G_0 p_0 = 0, linearised constraint
(kappa - A_0) S_0 A_p(0) = A_0 (C~_W W_p(0) + C~_V V_p(0)); the constraint surface is invariant, so
the series lies on it (``constraint_residual`` checks this as a ball identity).  For non-resonant
kappa the analytic solution with these normalisations is unique, hence it coincides with the
Stage-1 (3D) series wherever both exist.
Taylor model in delta = kappa - kappa_c exactly as in Stage 1 (point run in delta-polynomial
arithmetic, remainder from the (m+1)-th coefficient of the run with the box as base point);
tail certificate as in Stage 1 with the rank-1 sharpening ||M_n^{-1}|| <= c/n,
c = ||D^{-1}|| (1 + ||D^{-1}E|| / (K + 1 - max Re sigma)), valid for every kappa in the box.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from flint import acb, acb_mat, arb

from .arbseries import Series, l1nu, precision, to_arb
from .linsonic import NU_CANDIDATES, _arb_norm_inf, norm_inf
from .linsys import LinSystem, abs_up, kappa_box, to_acb
from .polysys import _PolyEvaluator, abs_eval, abs_eval_increment
from .shootsys import eval_box
from .tailbound import Certificate


def _col(v):
    return acb_mat([[to_acb(x)] for x in v])


def _zero_col():
    return acb_mat(4, 1)


def _upper(a):
    """arb >= sup of the real ball a."""
    return arb(a.mid()) + arb(a.rad())


def polys4(L):
    """Exact 4D data: P (kappa^0), Ga = [DQ, -P_s] (kappa^0, kappa^1), Gb[r][c][i] = -dP_{ri}/du_c."""
    sys = L.sys
    dP, dQ = sys.dP(), sys.dQ()
    P = [[sys.P[r][c] for c in range(4)] for r in range(4)]
    Ga = [[[dQ[r][c] for c in range(4)] for r in range(4)], [[-L.Ps[r][c] for c in range(4)] for r in range(4)]]
    Gb = [[[-dP[r][i][c] for i in range(4)] for c in range(4)] for r in range(4)]
    return P, Ga, Gb


def series_parts4(polys, u, cap=None):
    """(Pser, [G0ser, G1ser]) 4x4 real Series matrices along the background series u (4 Series)."""
    P, Ga, Gb = polys
    u = [ui.with_cap(cap) for ui in u]
    ev = _PolyEvaluator([Series.var(cap)] + u)
    th = [ui.deriv() for ui in u]
    Pser = [[ev(P[r][c]) for c in range(4)] for r in range(4)]
    G0 = []
    for r in range(4):
        row = []
        for c in range(4):
            acc = ev(Ga[0][r][c])
            for i in range(4):
                if not Gb[r][c][i].is_zero():
                    acc = acc + ev(Gb[r][c][i]) * th[i]
            row.append(acc)
        G0.append(row)
    G1 = [[ev(Ga[1][r][c]) for c in range(4)] for r in range(4)]
    return Pser, [G0, G1]


class _Structure4:
    """P[k] (acb_mat, x-order k) and G[dg][k] (delta-degree dg <= 1) for kappa = base + delta."""

    def __init__(self, Pser, Gser, base):
        k = to_acb(base)
        self.lp = max(len(s) for row in Pser for s in row)
        self.lg = max(len(s) for M in Gser for row in M for s in row)
        ent = lambda M, n: acb_mat([[M[r][c][n] for c in range(4)] for r in range(4)])   # noqa: E731
        self.P = [ent(Pser, n) for n in range(self.lp)]
        G0 = [ent(Gser[0], n) for n in range(self.lg)]
        G1 = [ent(Gser[1], n) for n in range(self.lg)]
        self.G = [[G0[n] + k * G1[n] for n in range(self.lg)], G1]

    def DE(self, rows):
        D = [[acb(0)] * 4 for _ in range(4)]
        E = [[[acb(0)] * 4 for _ in range(4)] for _ in range(2)]
        for i, row in enumerate(rows):
            for c, r, s in row:
                for l in range(4):
                    D[i][l] += c * self.P[s + 1][r, l]
                if s >= 0:
                    for dg in range(2):
                        for l in range(4):
                            E[dg][i][l] -= c * self.G[dg][s][r, l]
        return acb_mat(D), [acb_mat(e) for e in E]

    def residual(self, q, j, upto, dm):
        """delta-polynomial (list over dg < dm of 4x1) of [P p' - G p]_j using p_m, m <= upto."""
        out = [_zero_col() for _ in range(dm)]
        for i in range(j + 1):
            if i + 1 <= upto and j - i < self.lp:
                for b in range(dm):
                    out[b] += self.P[j - i] * q[i + 1][b] * (i + 1)
            if i <= upto and j - i < self.lg:
                for a in range(2):
                    for b in range(dm - a):
                        out[a + b] -= self.G[a][j - i] * q[i][b]
        return out

    def level_values(self, q, n, rows, upto, dm):
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


def _run4(L, st, u0, rows, base, K, dm):
    """Coefficients p_n (list over n of lists over dg < dm of 4x1 acb_mat) for kappa = base + delta."""
    args = [arb(0)] + list(u0)
    A0, S0 = u0[0], eval_box(L.S, args)
    CW, CV = eval_box(L.dC[2], args), eval_box(L.dC[3], args)
    ell = [c for (c, _, _) in rows[3]]
    D, E = st.DE(rows)
    lG = lambda dg, l: ell[0] * st.G[dg][0][2, l] + ell[1] * st.G[dg][0][3, l]      # noqa: E731
    M0 = [acb_mat([[lG(0, 2), lG(0, 3)], [A0 * CW, A0 * CV]]), acb_mat([[lG(1, 2), lG(1, 3)], [acb(0), acb(0)]])]
    rhs = [acb_mat([[-lG(0, 0)], [(to_acb(base) - A0) * S0]]), acb_mat([[-lG(1, 0)], [to_acb(S0)]])]
    wv = []
    for k in range(dm):
        b = rhs[k] if k < 2 else acb_mat(2, 1)
        if k >= 1:
            b = b - M0[1] * wv[k - 1]
        wv.append(M0[0].solve(b))
    q = [[acb_mat([[acb(int(k == 0))], [acb(0)], [wv[k][0, 0]], [wv[k][1, 0]]]) for k in range(dm)]]
    for n in range(1, K + 1):
        R = st.level_values(q, n, rows, n - 1, dm)
        M = [D * n + E[0], E[1]]
        qn = []
        for k in range(dm):
            b = -R[k]
            if k >= 1:
                b = b - M[1] * qn[k - 1]
            qn.append(M[0].solve(b))
        q.append(qn)
    return q, D, E


@dataclass
class LinSonicExpansion4:
    kappa_c: acb
    width: float
    m: int
    K: int
    L: LinSystem
    polys: tuple
    bg: object
    rows: list
    point: dict
    box: dict = None
    balls: list = None                 # list over n of 4 acb: (A_p, N_p, W_p, V_p)_n over the box
    cert: object = None
    info: dict = field(default_factory=dict)

    @property
    def kappa(self):
        return kappa_box(self.kappa_c, self.width) if self.width else self.kappa_c

    def certify(self, nu=None, prec=256):
        with precision(prec):
            self.cert = certify4(self, nu=nu)
        return self.cert

    def eval(self, x, with_tail=True):
        """Balls (A_p, N_p, W_p, V_p)(x) over the box (tail included if certified, |x| <= nu)."""
        x = to_arb(x)
        out = []
        for i in range(4):
            acc = acb(0)
            for n in reversed(range(self.K + 1)):
                acc = acc * x + self.balls[n][i]
            out.append(acc)
        if with_tail:
            tb = self.cert.tail_bound(abs(x).abs_upper())
            out = [v + acb(arb(0, tb), arb(0, tb)) for v in out]
        return out

    def delta_model(self, x, r):
        """p(x; kappa_c + delta) = sum_{k<=m} c_k delta^k + R,  |R_i| <= rem_i for |delta| <= r <= width
        (remainder: (m+1)-th delta-coefficient over the box, plus the box tail bound at |x|)."""
        x = to_arb(x)
        if not (self.width and r <= self.width and self.cert is not None and self.cert.ok):
            raise ValueError("delta_model needs a certified box expansion with r <= width")
        coefs = []
        for k in range(self.m + 1):
            ck = []
            for i in range(4):
                acc = acb(0)
                for n in reversed(range(self.K + 1)):
                    acc = acc * x + self.point["coefs"][n][k][i, 0]
                ck.append(acc)
            coefs.append(ck)
        ax, rp = abs(x).abs_upper(), to_arb(r) ** (self.m + 1)
        tb = self.cert.tail_bound(ax)
        rem = []
        for i in range(4):
            top = sum((abs_up(self.box["coefs"][n][self.m + 1][i, 0]) * ax**n for n in range(self.K + 1)), arb(0))
            rem.append((top * rp + tb).abs_upper())
        return coefs, rem


def linear_sonic_expansion4(bg, kappa_c, width=0.0, m=5, K=40, prec=256):
    """4D sonic series for kappa in kappa_c + [-w, w] + i[-w, w] (w = 0: point); bg: certified A1
    expansion with bg.K >= K + 1 (``linmatch.box_background`` for the V0 box)."""
    if bg.K < K + 1:
        raise ValueError("background expansion must have K_bg >= K + 1")
    with precision(prec):
        L = LinSystem()
        polys = polys4(L)
        kc = to_acb(kappa_c)
        Pser, Gser = series_parts4(polys, bg.series())
        ell, row = bg.point["ell"], bg.point["row"]
        rows = [[(arb(1), 0, -1)], [(arb(1), 1, -1)], [(arb(1), row, -1)], [(ell[0], 2, 0), (ell[1], 3, 0)]]
        u0, dm = bg.balls[0], m + 2
        st_pt = _Structure4(Pser, Gser, kc)
        q_pt, D_pt, E_pt = _run4(L, st_pt, u0, rows, kc, K, dm)
        point = dict(coefs=q_pt, D=D_pt, E=E_pt, structure=st_pt)
        box = None
        if width:
            kb = kappa_box(kc, width)
            st_bx = _Structure4(Pser, Gser, kb)
            q_bx, D_bx, E_bx = _run4(L, st_bx, u0, rows, kb, K, dm)
            box = dict(coefs=q_bx, D=D_bx, E=E_bx, structure=st_bx)
            dw = acb(arb(0, width), arb(0, width))
            wpow = abs_up(dw) ** (m + 1)
            balls = []
            for n in range(K + 1):
                bn = []
                for i in range(4):
                    acc = acb(0)
                    for k in reversed(range(m + 1)):
                        acc = acc * dw + q_pt[n][k][i, 0]
                    top = abs_up(q_bx[n][m + 1][i, 0]) * wpow
                    bn.append(acc + acb(arb(0, top), arb(0, top)))
                balls.append(bn)
        else:
            balls = [[q_pt[n][0][i, 0] for i in range(4)] for n in range(K + 1)]
        return LinSonicExpansion4(kc, float(width), m, K, L, polys, bg, rows, point, box, balls)


def sigma(D, E0):
    """The Frobenius exponent sigma(kappa) = -tr(D^{-1} E(kappa)) (acb; box if E0 is over a box)."""
    X = D.inv() * E0
    return -sum((X[i, i] for i in range(4)), acb(0))


def certify4(ex, nu=None):
    """Tail certificate of ``ex.balls`` (n > K) valid for every kappa in the box; the argument of
    Stage 1 (notes s2-theorem-b.md section 1.3) with 4 rows and the rank-1 bound on ||M_n^{-1}||."""
    L, K, rows = ex.L, ex.K, ex.rows
    P, Ga, Gb = ex.polys
    run = ex.box or ex.point
    st, D, E0 = run["structure"], run["D"], run["E"][0]
    cert_u = ex.bg.cert
    if cert_u is None or not cert_u.ok:
        raise ValueError("background expansion must be certified first")
    if not all(E0[i, l].contains(acb(0)) and E0[i, l].rad() == 0 for i in range(3) for l in range(4)):
        raise ValueError("E is expected to have a single nonzero row")
    Dinv = D.inv()
    g, Dn = norm_inf(Dinv * E0), norm_inf(Dinv)
    sg = sigma(D, E0)
    gap = K + 1 - _upper(sg.real)
    if not gap > 0:
        return Certificate(K, arb(0), arb(0), arb(0), arb(0), arb(0), arb(0), False, dict(reason="Re sigma >= K+1"))
    c = Dn * (1 + g / gap)
    ubar = ex.bg.series()
    th_ubar = [s.deriv() for s in ubar]
    P0, G0 = st.P, st.G[0]
    kabs = abs_up(ex.kappa)
    q1 = [[_col(b)] for b in ex.balls]
    jmax = st.lp + K + 1
    F = [st.residual(q1, j, K, 1)[0] for j in range(jmax + 1)]

    def En(n):
        out = [arb(0)] * 4
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
            Bk = [[arb(0)] * 4 for _ in range(4)]
            for i, row in enumerate(rows):
                for cc, r, s in row:
                    for l in range(4):
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
        incP = [[abs_eval_increment(P[r][q], args, eps_v) for q in range(4)] for r in range(4)]
        incG = [[arb(0)] * 4 for _ in range(4)]
        for r in range(4):
            for q in range(4):
                term = abs_eval_increment(Ga[0][r][q], args, eps_v) + kabs * abs_eval_increment(Ga[1][r][q], args, eps_v)
                for i in range(4):
                    if not Gb[r][q][i].is_zero():
                        term += abs_eval_increment(Gb[r][q][i], args, eps_v) * thu[i]
                        term += abs_eval(Gb[r][q][i], args_e) * eps_th
                incG[r][q] = term
        Zm = [[arb(0)] * 4 for _ in range(4)]
        Y2 = arb(0)
        nq = [sum((abs_up(ex.balls[n][l]) * nu_**n for n in range(K + 1)), arb(0)) for l in range(4)]
        nth = [sum((abs_up(ex.balls[n][l]) * n * nu_ ** (n - 1) for n in range(1, K + 1)), arb(0)) for l in range(4)]
        for i, row in enumerate(rows):
            for cc, r, s in row:
                ac = abs_up(cc)
                for l in range(4):
                    Zm[i][l] += ac * (nu_ ** (-(s + 1)) * incP[r][l] + nu_ ** (-s) * incG[r][l] / (K + 2))
                    Y2 += ac * nu_ ** (-s) * (incP[r][l] * nth[l] + incG[r][l] * nq[l])
        Z2 = c * _arb_norm_inf(Zm)
        Y2 = c * Y2 / (K + 1)
        Z, Y = Z1 + Z2, Y1 + Y2
        det = dict(g=g, Dinv_norm=Dn, sigma=sg, eps_v=eps_v, eps_th=eps_th, Y1=Y1, Y2=Y2)
        if Z < 1:
            return Certificate(K, nu_, (Y / (1 - Z)).abs_upper(), Y, Z1, Z2, c, True, det)
        last = Certificate(K, nu_, arb(0), Y, Z1, Z2, c, False, dict(reason="Z >= 1", **det))
    return last


def constraint_residual(ex, upto=None):
    """Coefficients (n <= upto) of  (kappa - A) S A_p - A (C~_N N_p + C~_W W_p + C~_V V_p)  along the
    truncated background: balls that must contain 0 (constraint-surface invariance, ball identity)."""
    L, K = ex.L, ex.K
    upto = K if upto is None else upto
    ubar = [s.with_cap(upto + 1) for s in ex.bg.series()]
    ev = _PolyEvaluator([Series.var(upto + 1)] + ubar)
    S, A = ev(L.S), ubar[0]
    Cq = [ev(L.dC[l + 1]) for l in range(3)]
    kap = ex.kappa
    p = [Series([b[i] for b in ex.balls], cap=upto + 1) for i in range(4)] if False else None
    out = []
    for n in range(upto + 1):
        acc = acb(0)
        for j in range(n + 1):
            SA = sum((S[j - i] * A[i] for i in range(j + 1)), arb(0))
            acc += kap * S[j] * ex.balls[n - j][0] - SA * ex.balls[n - j][0]
            for l in range(3):
                ACl = sum((A[i] * Cq[l][j - i] for i in range(j + 1)), arb(0))
                acc -= ACl * ex.balls[n - j][l + 1]
        out.append(acc)
    return out
