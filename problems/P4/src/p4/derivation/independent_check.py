"""Independent pointwise verification of the CSS reduction in exact rational arithmetic.

Independence from ``einstein_euler`` (deliverable 1): no sympy, no symbolic differentiation,
no symbolic substitution, no shared code.  Everything is evaluated at random rational points
of the *jet space* -- values and first/second partials of (alpha, a, rho, V) at a point
(t0, r0, theta0) -- with forward-mode AD over Q (``qjet``); an identity between rational
functions that holds at such a point holds exactly there, so failures cannot hide behind
simplification or rounding.  The Einstein tensor is recomputed from the metric by the
coordinate formulas *on jets*; the KHA rows are transcribed here again from KHA95 eq. 18 as
displayed in notes/literature-digest.md (with S1's row-3 correction) in Fraction arithmetic;
rho is carried as rhob = 4 pi rho so that no pi occurs (G = 2 * 4 pi T).
Checks at each point (``check_point``):
  (a) E_tt = N^2 e^{2s} (A_x/A - F_A);  E_tr = e^{2s-x} ((A_s+A_x)/A - G);
      E_rr = e^{2s-2x} (2 (N_x/N - F_N) + (A_x/A - F_A))          [hand-derived factors];
      since the second-order jet coefficients are random, this also shows they cancel;
  (b) alpha_t does not enter the Euler equations;
  (c) with a_r, a_t, alpha_r solved from E_tt = E_tr = E_rr = 0 (each is affine in its
      variable), the Euler pair equals M (row3, row4) with M reconstructed from the affine
      structure in (rho_t, V_t) and verified on the (rho_r, V_r) columns and the sources;
      det M != 0;  row1 = row2 = rowM = 0 there;
  (d) Misner-Sharp mass m = r (1 - a^{-2})/2 obeys m_t = r^2 T_t^r, m_r = -r^2 T^t_t
      (Hayward's covariant form of HM01 eqs. 3-4) at such points.
"""
from __future__ import annotations

import random
from fractions import Fraction as Fr

from .qjet import QJet

NV, ORD = 4, 2
T, R, TH, PH = range(4)
SIN0, COS0 = Fr(3, 5), Fr(4, 5)            # theta0 with rational sin, cos
E_T, E_R = (1, 0, 0, 0), (0, 1, 0, 0)


def _set(jet, e, v):
    c = dict(jet.c)
    c[e] = Fr(v)
    return QJet(jet.n, jet.order, c)


class Point:
    """Coordinates (t0, r0, theta0) and the four field jets (functions of t, r only)."""

    def __init__(self, t0, r0, alpha, a, rhob, V):
        self.t0, self.r0 = Fr(t0), Fr(r0)
        self.alpha, self.a, self.rhob, self.V = alpha, a, rhob, V
        self.t = QJet.var(NV, ORD, T, self.t0)
        self.r = QJet.var(NV, ORD, R, self.r0)
        self.sin = QJet(NV, ORD, {(0,) * 4: SIN0, (0, 0, 1, 0): COS0, (0, 0, 2, 0): -SIN0 / 2})

    def copy(self):
        return Point(self.t0, self.r0, self.alpha, self.a, self.rhob, self.V)

    def with_coeff(self, attr, e, v):
        q = self.copy()
        setattr(q, attr, _set(getattr(q, attr), e, v))
        return q


def random_point(rng):
    def pos():
        return Fr(rng.randint(1, 30), rng.randint(1, 10))

    def field(v0):
        c = {(0,) * 4: v0}
        for e in [E_T, E_R, (2, 0, 0, 0), (1, 1, 0, 0), (0, 2, 0, 0)]:
            c[e] = Fr(rng.randint(-40, 40), rng.randint(1, 12))
        return QJet(NV, ORD, c)

    return Point(-pos(), pos(), field(pos()), field(pos()), field(pos()),
                 field(Fr(rng.randint(-9, 9), 10)))


# --- geometry and fluid on jets ------------------------------------------------------------
def metric(p):
    g = [-p.alpha**2, p.a**2, p.r**2, p.r**2 * p.sin**2]
    return g, [1 / x for x in g]


def christoffel(g, ginv):
    """Gamma^l_{mn} = g^{ll} (d_m g_{ln} + d_n g_{lm} - d_l g_{mn}) / 2 (diagonal g)."""
    Gam = [[[None] * 4 for _ in range(4)] for _ in range(4)]
    for l in range(4):
        for m in range(4):
            for n in range(4):
                tot = QJet.const(NV, ORD - 1, 0)
                if l == n:
                    tot = tot + g[l].d(m)
                if l == m:
                    tot = tot + g[l].d(n)
                if m == n:
                    tot = tot - g[m].d(l)
                Gam[l][m][n] = ginv[l] * tot * Fr(1, 2)
    return Gam


def ricci(Gam):
    Rc = [[Fr(0)] * 4 for _ in range(4)]
    for m in range(4):
        for n in range(4):
            tot = Fr(0)
            for l in range(4):
                tot += Gam[l][m][n].d(l).val - Gam[l][m][l].d(n).val
                for s in range(4):
                    tot += Gam[l][l][s].val * Gam[s][m][n].val - Gam[l][n][s].val * Gam[s][m][l].val
            Rc[m][n] = tot
    return Rc


def einstein(g, ginv, Rc):
    Rs = sum(ginv[i].val * Rc[i][i] for i in range(4))
    return [[Rc[m][n] - (g[m].val * Rs / 2 if m == n else 0) for n in range(4)] for m in range(4)]


def stress_up(p, ginv):
    """Tbar^{mn} = 4 pi T^{mn}, p = rho/3, u^t = 1/(alpha sqrt(1-V^2)), u^r = V/(a sqrt(1-V^2))."""
    pb = p.rhob * Fr(1, 3)
    W2 = 1 / (1 - p.V**2)
    Tup = [[QJet.const(NV, ORD, 0) for _ in range(4)] for _ in range(4)]
    Tup[T][T] = (p.rhob + pb) * W2 / p.alpha**2 + pb * ginv[T]
    Tup[T][R] = Tup[R][T] = (p.rhob + pb) * p.V * W2 / (p.alpha * p.a)
    Tup[R][R] = (p.rhob + pb) * p.V**2 * W2 / p.a**2 + pb * ginv[R]
    Tup[TH][TH] = pb * ginv[TH]
    Tup[PH][PH] = pb * ginv[PH]
    return Tup


def divergence(Tup, Gam):
    U = []
    for n in range(4):
        tot = Fr(0)
        for m in range(4):
            tot += Tup[m][n].d(m).val
            for l in range(4):
                tot += Gam[m][m][l].val * Tup[l][n].val + Gam[n][m][l].val * Tup[m][l].val
        U.append(tot)
    return U


def equations(p):
    """E_{mn} = G_{mn} - 2 Tbar_{mn} for (tt, tr, rr) and U^n = nabla_m Tbar^{mn}."""
    g, ginv = metric(p)
    Gam = christoffel(g, ginv)
    G = einstein(g, ginv, ricci(Gam))
    Tup = stress_up(p, ginv)
    E = {k: G[k[0]][k[1]] - 2 * (g[k[0]] * g[k[1]] * Tup[k[0]][k[1]]).val
         for k in ((T, T), (T, R), (R, R))}
    return E, divergence(Tup, Gam), Tup, g


# --- KHA side ---------------------------------------------------------------------------
def kha_jets(p):
    """(A, N, W, V) and their (s, x)-derivatives at the point: d_x = r d_r, d_s = -t d_t - r d_r."""
    A, N, W = p.a**2, p.alpha * (-p.t) / (p.a * p.r), p.r**2 * p.a**2 * p.rhob
    out = {}
    for name, f in (("A", A), ("N", N), ("W", W), ("V", p.V)):
        fx = p.r0 * f.d(R).val
        out[name], out[name + "x"], out[name + "s"] = f.val, fx, -p.t0 * f.d(T).val - fx
    out["es"], out["ex"] = 1 / (-p.t0), p.r0 / (-p.t0)
    return out


def kha_rows(k):
    """KHA95 eq. 18 (digest 1.1, row 3 with S1's 2NV(...) correction) and KHA99 eq. 211."""
    A, N, W, V = k["A"], k["N"], k["W"], k["V"]
    S = 1 - V * V
    FA = 1 - A + 2 * W * (1 + V * V / 3) / S
    FN = -2 + A - 2 * W / 3
    G = -Fr(8, 3) * N * V * W / S
    return {
        "row1": k["Ax"] / A - FA, "row2": k["Nx"] / N - FN, "rowM": (k["As"] + k["Ax"]) / A - G,
        "row3": (k["Ws"] + (1 + N * V) * k["Wx"]) / W + 4 * (V * k["Vs"] + (N + V) * k["Vx"]) / (3 * S)
        - N * V * FA / 3 + 4 * V * N * FN / 3 + 2 * N * V * (1 + 4 * W / (9 * S)),
        "row4": (4 * V * k["Ws"] + (4 * V + N + 3 * N * V * V) * k["Wx"]) / W
        + 4 * ((1 + V * V) * k["Vs"] + (1 + V * V + 2 * N * V) * k["Vx"]) / S
        + N * S * FA + 4 * (1 + V * V) * N * FN + 2 * N * (1 + 3 * V * V)}


# --- checks ---------------------------------------------------------------------------------
def _root(p, attr, e, key):
    """Zero in the coefficient (attr, e) of the affine map v -> E[key](p with coeff v)."""
    v0, v1 = (equations(p.with_coeff(attr, e, v))[0][key] for v in (0, 1))
    return -v0 / (v1 - v0)


def constrain(p):
    """Solve E_tt, E_tr, E_rr = 0 for a_r, a_t, alpha_r (each affine in its own variable)."""
    q = p.with_coeff("a", E_R, _root(p, "a", E_R, (T, T)))
    q = q.with_coeff("a", E_T, _root(q, "a", E_T, (T, R)))
    q = q.with_coeff("alpha", E_R, _root(q, "alpha", E_R, (R, R)))
    E = equations(q)[0]
    assert all(v == 0 for v in E.values()), E
    return q


def check_factors(p):
    """(a): closed-form factors between E_tt, E_tr, E_rr and row1, rowM, row2 (unconstrained p)."""
    E = equations(p)[0]
    k = kha_jets(p)
    rw = kha_rows(k)
    es, ex = k["es"], k["ex"]
    return {"E_tt": E[(T, T)] == k["N"] ** 2 * es**2 * rw["row1"],
            "E_tr": E[(T, R)] == es**2 / ex * rw["rowM"],
            "E_rr": E[(R, R)] == es**2 / ex**2 * (2 * rw["row2"] + rw["row1"])}


def check_lapse_independence(p):
    """(b): the Euler equations do not depend on alpha_t."""
    U0 = equations(p)[1]
    U1 = equations(p.with_coeff("alpha", E_T, p.alpha.coeff(E_T) + 1))[1]
    return U0[T] == U1[T] and U0[R] == U1[R]


def euler_vs_rows(q):
    """(c): at a constrained point, (U^t, U^r) = M (row3, row4) with M from the (rho_t, V_t)
    columns of the affine structure; returns (ok, M, constraint rows)."""
    def ev(pt):
        _, U, _, _ = equations(pt)
        rw = kha_rows(kha_jets(pt))
        return (U[T], U[R]), (rw["row3"], rw["row4"]), (rw["row1"], rw["row2"], rw["rowM"])

    U0, R0, C0 = ev(q)
    cols = []
    for attr, e in (("rhob", E_T), ("rhob", E_R), ("V", E_T), ("V", E_R)):
        base = getattr(q, attr).coeff(e)
        U1, R1, _ = ev(q.with_coeff(attr, e, base + 1))
        U2, R2, _ = ev(q.with_coeff(attr, e, base + 2))
        dU, dR = tuple(U1[i] - U0[i] for i in range(2)), tuple(R1[i] - R0[i] for i in range(2))
        affine = all(U2[i] - U0[i] == 2 * dU[i] and R2[i] - R0[i] == 2 * dR[i] for i in range(2))
        cols.append((dU, dR, affine))
    LU = [[cols[j][0][i] for j in (0, 2)] for i in range(2)]
    LR = [[cols[j][1][i] for j in (0, 2)] for i in range(2)]
    det = LR[0][0] * LR[1][1] - LR[0][1] * LR[1][0]
    assert det != 0
    LRi = [[LR[1][1] / det, -LR[0][1] / det], [-LR[1][0] / det, LR[0][0] / det]]
    M = [[sum(LU[i][k] * LRi[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    ok = all(c[2] for c in cols)
    ok &= all(cols[j][0][i] == sum(M[i][k] * cols[j][1][k] for k in range(2))
              for j in (1, 3) for i in range(2))
    ok &= all(U0[i] == sum(M[i][k] * R0[k] for k in range(2)) for i in range(2))
    ok &= (M[0][0] * M[1][1] - M[0][1] * M[1][0]) != 0
    ok &= all(c == 0 for c in C0)
    return ok, M, C0


def check_misner_sharp(q):
    """(d): m = r(1 - a^{-2})/2 has m_t = r^2 Tbar_t^r, m_r = -r^2 Tbar^t_t at constrained points
    (4 pi absorbed in Tbar)."""
    _, _, Tup, g = equations(q)
    m = q.r * (1 - 1 / q.a**2) * Fr(1, 2)
    return (m.d(T).val == q.r0**2 * (g[T] * Tup[T][R]).val
            and m.d(R).val == -q.r0**2 * (g[T] * Tup[T][T]).val)


def check_point(p):
    q = constrain(p)
    ok_c, M, _ = euler_vs_rows(q)
    rep = dict(check_factors(p), lapse=check_lapse_independence(p), euler_rows=ok_c,
               misner_sharp=check_misner_sharp(q))
    rep["M"] = M
    return rep


def run(n_points=8, seed=20260829):
    """Run all checks at n_points random rational points; raise on any failure."""
    rng = random.Random(seed)
    reports = []
    for _ in range(n_points):
        rep = check_point(random_point(rng))
        bad = [k for k, v in rep.items() if k != "M" and v is not True]
        if bad:
            raise AssertionError(f"independent check failed: {bad}")
        reports.append(rep)
    return reports


def cross_check_css(n_points=8, seed=1, tol=1e-9):
    """Float consistency of this file's transcription with ``p4.css.fluid_residuals``
    (not part of the exact proof; the exact link css <-> rows is in einstein_euler)."""
    from p4 import css
    rng = random.Random(seed)
    worst = 0.0
    for _ in range(n_points):
        k = kha_jets(random_point(rng))
        rw = kha_rows(k)
        f = {n: float(v) for n, v in k.items()}
        r3, r4 = css.fluid_residuals(f["A"], f["N"], f["W"], f["V"], f["Wx"], f["Vx"], f["Ws"], f["Vs"])
        for mine, theirs in ((rw["row3"], r3), (rw["row4"], r4)):
            worst = max(worst, abs(float(mine) - theirs) / (1 + abs(theirs)))
    assert worst < tol, worst
    return worst
