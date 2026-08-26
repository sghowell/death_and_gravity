"""Rigorous verification of dual certificates for the kappa-model bound problems.

Everything here is computed in ball arithmetic (python-flint / Arb) from the frozen inputs;
no floating-point rounding is trusted on the verification path. A certificate is a dual vector z
for the conic program

    min q^T y   s.t.  A y + s = b,  s in K = {0}^n_eq x R_+^n_in x SOC(n_soc),

and weak duality with residual absorption gives, for every primal-feasible y,

    q^T y  >=  -b^T z  -  sum_k |rho_k| * Y_k,      rho := A^T z + q,   |y_k| <= Y_k,

provided z_in >= 0 and z_soc in SOC. The constraint matrices are re-derived here independently
of the solver code (socp2.py) from the frozen data and the class definition, so the check does
not depend on the correctness of the floating-point model, only on the relaxation being valid
(FORMULATION.md §3.2) — the solver merely proposes z.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from flint import arb, arb_mat, ctx

from .model import Frozen
from .socp2 import NodeBounds

ctx.prec = 160  # bits


def _a(x) -> arb:
    """Exact conversion of a float (a rational) to a ball."""
    return arb(float(x))


def _log10(x: arb) -> arb:
    return x.log() / arb(10).log()


def _max(a: arb, b: arb) -> arb:
    """Enclosure of max(a, b)."""
    return (a + b + abs(a - b)) / 2


def _up(x: arb) -> float:
    """Rigorous upper endpoint as a float (rounded up)."""
    return float(x.upper().str(30, radius=False).split(" +/- ")[0]) if False else _endpoint(x, +1)


def _endpoint(x: arb, sign: int) -> float:
    m, r = x.mid(), x.rad()
    e = m + r if sign > 0 else m - r
    # convert to float conservatively
    f = float(e.str(25, radius=False))
    # nudge outward by a few ulps to cover decimal-conversion rounding
    return np.nextafter(np.nextafter(f, np.inf), np.inf) if sign > 0 else np.nextafter(np.nextafter(f, -np.inf), -np.inf)


@dataclass
class ProblemDef:
    """Rigorous (ball-arithmetic) definition of one relaxed problem instance."""
    x: list            # grid nodes as arb (exact rationals)
    L: arb
    u_box: tuple       # (ulo, uhi) arb
    Mp_box: tuple
    T: arb             # rigorous upper bound on chi2_ref + Delta
    lo: list; hi: list # node brackets for D_i (i>=1) and u_0 (i=0), arb (rational)
    tp: list           # tangent points (rational)
    e_seg: list        # rigorous interpolation slacks (upper endpoints), arb


class Verifier:
    def __init__(self, fr: Frozen, nb: NodeBounds, T_float: float, tangent_u: np.ndarray,
                 T_rigorous: arb | None = None, Mp_box=(0.0, 40.0)):
        self.fr = fr
        sp = fr.spec
        self.N = sp.n_seg
        self.n = len(fr.sn.m)
        self.nb_ = len(fr.bao.value)
        self.x = [_a(v) for v in sp.x]
        self.L = _a(sp.L)
        self.ulo, self.uhi = _a(sp.u_box[0]), _a(sp.u_box[1])
        self.Mp_lo, self.Mp_hi = _a(Mp_box[0]), _a(Mp_box[1])
        self.T = T_rigorous if T_rigorous is not None else _a(T_float)
        # brackets (rational) and tangent points
        self.nb = nb
        self.lo = [_a(nb.u_lo[0])] + [_a(v) for v in nb.dm_lo[1:]]
        self.hi = [_a(nb.u_hi[0])] + [_a(v) for v in nb.dm_hi[1:]]
        from .socp2 import tangent_points
        lo_f = np.concatenate([[nb.u_lo[0]], nb.dm_lo[1:]]); hi_f = np.concatenate([[nb.u_hi[0]], nb.dm_hi[1:]])
        self.tps = [[_a(v) for v in tp] for tp in tangent_points(lo_f, hi_f, tangent_u, fr)]
        # --- rigorous geometry: exp(x_k) = 1 + z_k is NOT used (x_k are the exact rationals); use arb exp
        self.ex = [xi.exp() for xi in self.x]
        self.hs = [self.x[k + 1] - self.x[k] for k in range(self.N)]
        # A_nodes rows (D_i as linear map of u), rigorous
        self.A_nodes = self._dm_rows(self.x[1:])
        # SN geometry
        self.xj = [(_a(1.0) + _a(z)).log() for z in fr.sn.zHD]
        self.A_sn = self._dm_rows(self.xj)
        self.k_sn = [self._segment(xj) for xj in self.xj]
        self.t_sn = [(self.xj[j] - self.x[k]) / self.hs[k] for j, k in enumerate(self.k_sn)]
        self.c_sn = [_log10(xj.exp() - 1) for xj in self.xj]
        self.off = [5 * _log10(_a(1.0) + _a(z)) for z in fr.sn.zHEL]
        self.m = [_a(v) for v in fr.sn.m]
        # BAO geometry rows
        self.xb = [(_a(1.0) + _a(z)).log() for z in fr.bao.z]
        Pm = self._dm_rows(self.xb)
        Ph = self._dh_rows(self.xb)
        self.P = [Pm[r] if kind == "DM_over_rs" else Ph[r] for r, kind in enumerate(fr.bao.kind)]
        self.d = [_a(v) for v in fr.bao.value]
        # whiteners as exact rationals (their floats define C~)
        self.Wb = arb_mat([[_a(v) for v in row] for row in fr.Wb])
        self.Wsn = arb_mat([[_a(v) for v in row] for row in fr.Wsn])
        # interpolation slack (rigorous upper endpoints)
        self.e_seg = [self._slack(k) for k in range(self.N)]
        self.c_node = [arb(0)] + [_log10(self.ex[i] - 1) for i in range(1, self.N + 1)]

    # ---------------- geometry in ball arithmetic ----------------
    def _segment(self, xv: arb) -> int:
        # rigorous: find k with x_k <= xv <= x_{k+1}; ties resolved either way are valid (continuity)
        mid = float(xv.mid().str(20, radius=False))
        xs = self.fr.spec.x
        k = int(np.searchsorted(xs, mid, side="right") - 1)
        return int(np.clip(k, 0, self.N - 1))

    def _dm_rows(self, xs):
        rows = []
        for xv in xs:
            k = self._segment(xv)
            row = [arb(0)] * (self.N + 1)
            for i in range(k):
                h = self.hs[i]; ea, eb = self.ex[i], self.ex[i + 1]
                row[i] += (eb - (h + 1) * ea) / h
                row[i + 1] += ((h - 1) * eb + ea) / h
            h = self.hs[k]; xa, xb = self.x[k], self.x[k + 1]; ea = self.ex[k]; exv = xv.exp()
            row[k] += ((xb - xv + 1) * exv - (h + 1) * ea) / h
            row[k + 1] += ((xv - xa - 1) * exv + ea) / h
            rows.append(row)
        return rows

    def _dh_rows(self, xs):
        rows = []
        for xv in xs:
            k = self._segment(xv)
            t = (xv - self.x[k]) / self.hs[k]
            row = [arb(0)] * (self.N + 1)
            row[k] += 1 - t; row[k + 1] += t
            rows.append(row)
        return rows

    def _theta_bounds(self, xv: arb):
        L = self.L
        em1 = xv.exp() - 1

        def I(sign):
            a = 1 - sign * L
            if abs(float(a.mid().str(20, radius=False))) < 1e-9:
                return (sign * L * xv).exp() * xv
            return (sign * L * xv).exp() * ((a * xv).exp() - 1) / a

        return em1 / I(+1), em1 / I(-1)

    def _kappa_diff_bounds(self):
        """Rigorous (outward-rounded) class-only bounds on kappa_{i+1} - kappa_i; mirrors geometry.kappa_difference_bounds."""
        L = self.L
        N = self.N
        c = [arb(0)] + [_log10(self.ex[i] - 1) for i in range(1, N + 1)]

        def J(xa, xb, sign):
            a = 1 + sign * L
            if abs(float(a.mid().str(20, radius=False))) < 1e-9:
                return xa.exp() * (xb - xa)
            return xa.exp() * ((a * (xb - xa)).exp() - 1) / a

        def I(xa, sign):
            a = 1 - sign * L
            if abs(float(a.mid().str(20, radius=False))) < 1e-9:
                return (sign * L * xa).exp() * xa
            return (sign * L * xa).exp() * ((a * xa).exp() - 1) / a

        lo, hi = [], []
        for i in range(N):
            xa, xb = self.x[i], self.x[i + 1]
            Jm, Jp = J(xa, xb, -1), J(xa, xb, +1)
            if i == 0:
                r_lo, r_hi = Jm, Jp
            else:
                r_lo, r_hi = 1 + Jm / I(xa, +1), 1 + Jp / I(xa, -1)
            lo.append(_a(_endpoint(_log10(r_lo) - (c[i + 1] - c[i]), -1)))
            hi.append(_a(_endpoint(_log10(r_hi) - (c[i + 1] - c[i]), +1)))
        return lo, hi

    def _bracket_diff_bounds(self):
        """Rigorous bracket-aware bounds on kappa_{i+1} - kappa_i; mirrors socp2.bracket_kappa_difference_bounds."""
        N = self.N
        nb = self.nb
        c = self.c_node
        lo, hi = [], []
        for i in range(N):
            h = self.hs[i]; ea, eb = self.ex[i], self.ex[i + 1]
            a = (eb - (h + 1) * ea) / h
            b = ((h - 1) * eb + ea) / h
            ulo_i, uhi_i = _a(nb.u_lo[i]), _a(nb.u_hi[i])
            ulo_n, uhi_n = _a(nb.u_lo[i + 1]), _a(nb.u_hi[i + 1])
            if i == 0:
                r_lo = a + b * ulo_n / uhi_i
                r_hi = a + b * uhi_n / ulo_i
            else:
                inc_lo = a * ulo_i + b * ulo_n; inc_hi = a * uhi_i + b * uhi_n
                r_lo = 1 + inc_lo / _a(nb.dm_hi[i]); r_hi = 1 + inc_hi / _a(nb.dm_lo[i])
            lo.append(_a(_endpoint(_log10(r_lo) - (c[i + 1] - c[i]), -1)))
            hi.append(_a(_endpoint(_log10(r_hi) - (c[i + 1] - c[i]), +1)))
        return lo, hi

    def _slack(self, k: int) -> arb:
        L = self.L
        xa, xb, h = self.x[k], self.x[k + 1], self.hs[k]
        if k == 0:
            B = L / 6 + L * L * (arb(1) / 2 + xb / 6) ** 2
        else:
            q = xa.exp() / (xa.exp() - 1)
            th_lo, th_hi = self._theta_bounds(xb)
            dth = _max(abs(th_lo - 1), abs(th_hi - 1))
            B = q * (L * th_hi + dth * (q * (1 + th_hi) - 1))
        e = h * h / 8 * B / arb(10).log()
        return _a(_endpoint(e, +1))   # rigorous upper endpoint as a rational

    # ---------------- certificate check ----------------
    def certify(self, z_eq: np.ndarray, z_in: np.ndarray, z_soc: np.ndarray, q_u: np.ndarray,
                verbose: bool = True, q_k: np.ndarray | None = None) -> float:
        """Return a rigorous lower bound on q_u . u over the relaxed feasible set (hence over F).

        Row order must match socp2.KappaModel: eq = [BAO (nb), SN (n)]; ineq = [box lo (N+1), box hi (N+1),
        slopes (4N) if finite L, SN interp upper (n), SN interp lower (n), bracket lo (N+1), bracket hi (N+1),
        secant (N+1), tangent (N+1), Mp box lo (1), Mp box hi (1)]; soc = [t, w_b (nb), w_s (n)].
        """
        N, n, nb_ = self.N, self.n, self.nb_
        # ---- dual cone membership (exact on the float components) ----
        if np.any(z_in < 0):
            raise ValueError("dual inequality multipliers must be nonnegative")
        z_soc = z_soc.copy()
        nrm = arb(0)
        for v in z_soc[1:]:
            nrm += _a(v) * _a(v)
        nrm = nrm.sqrt()
        if not (_a(z_soc[0]) >= nrm):
            z_soc[0] = _endpoint(nrm, +1)   # push into the cone (valid: only loosens)
        # ---- accumulate rho = A^T z + q and val = -b^T z, blockwise ----
        rho_u = [arb(0)] * (N + 1); rho_k = [arb(0)] * (N + 1); rho_M = arb(0)
        rho_l = [arb(0)] * n; rho_wb = [arb(0)] * nb_; rho_ws = [arb(0)] * n
        val = arb(0)
        # equalities: BAO: w_b + Wb P u = Wb d
        zb = [_a(v) for v in z_eq[:nb_]]
        WbP = [[sum((self.Wb[r, s] * self.P[s][i] for s in range(nb_)), arb(0)) for i in range(N + 1)] for r in range(nb_)]
        Wbd = [sum((self.Wb[r, s] * self.d[s] for s in range(nb_)), arb(0)) for r in range(nb_)]
        for r in range(nb_):
            for i in range(N + 1):
                rho_u[i] += WbP[r][i] * zb[r]
            rho_wb[r] += zb[r]
            val -= Wbd[r] * zb[r]
        # SN: w_s + 5 Wsn ell + (Wsn 1) Mp = Wsn (m - off)
        zs = arb_mat([[_a(v)] for v in z_eq[nb_:nb_ + n]])
        WsnT_zs = self.Wsn.transpose() * zs                     # n x 1
        ones = arb_mat([[arb(1)]] * n)
        Wsn1 = self.Wsn * ones
        rhs = self.Wsn * arb_mat([[self.m[j] - self.off[j]] for j in range(n)])
        for j in range(n):
            rho_l[j] += 5 * WsnT_zs[j, 0]
            rho_ws[j] += zs[j, 0]
            rho_M += Wsn1[j, 0] * zs[j, 0]
            val -= rhs[j, 0] * zs[j, 0]
        # inequalities G y <= h : contribute G^T z_in to rho and -h^T z_in to val
        p = 0

        def take(cnt):
            nonlocal p
            block = z_in[p:p + cnt]; p += cnt
            return block

        # box lo: -u <= -ulo ; box hi: u <= uhi
        zl = take(N + 1); zh = take(N + 1)
        for i in range(N + 1):
            rho_u[i] += -_a(zl[i]) + _a(zh[i])
            val -= (-self.ulo) * _a(zl[i]) + self.uhi * _a(zh[i])
        if np.isfinite(self.fr.spec.L):
            # rows: d - Lh u_k <= 0 ; d - Lh u_{k+1} <= 0 ; -d - Lh u_k <= 0 ; -d - Lh u_{k+1} <= 0, d = u_{k+1}-u_k
            for variant in range(4):
                zz = take(N)
                for k in range(N):
                    zk = _a(zz[k]); Lh = self.L * self.hs[k]
                    sgn = 1 if variant < 2 else -1
                    rho_u[k + 1] += sgn * zk
                    rho_u[k] += -sgn * zk
                    if variant in (0, 2):
                        rho_u[k] += -Lh * zk
                    else:
                        rho_u[k + 1] += -Lh * zk
        # SN interpolation upper: ell_j - (1-t)k_k - t k_{k+1} <= c_j + e ; lower: -ell_j + ... <= -c_j + e
        zu = take(n); zlo = take(n)
        for j in range(n):
            k = self.k_sn[j]; t = self.t_sn[j]; e = self.e_seg[k]
            a1, a2 = _a(zu[j]), _a(zlo[j])
            rho_l[j] += a1 - a2
            rho_k[k] += -(1 - t) * a1 + (1 - t) * a2
            rho_k[k + 1] += -t * a1 + t * a2
            val -= (self.c_sn[j] + e) * a1 + (-self.c_sn[j] + e) * a2
        # bracket rows: -D_i <= -lo_i ; D_i <= hi_i  (D_0 = u_0, D_i = A_nodes[i-1] u)
        def Drow(i):
            if i == 0:
                r = [arb(0)] * (N + 1); r[0] = arb(1); return r
            return self.A_nodes[i - 1]

        zbl = take(N + 1); zbh = take(N + 1)
        for i in range(N + 1):
            r = Drow(i); a1, a2 = _a(zbl[i]), _a(zbh[i])
            for s in range(N + 1):
                rho_u[s] += (-a1 + a2) * r[s]
            val -= (-self.lo[i]) * a1 + self.hi[i] * a2
        # secant: -kappa_i + slope_i D_i <= -log10(lo_i) + slope_i lo_i + c_i
        zsec = take(N + 1)
        for i in range(N + 1):
            lo, hi = self.lo[i], self.hi[i]
            slope = (_log10(hi) - _log10(lo)) / (hi - lo)
            a1 = _a(zsec[i]); r = Drow(i)
            rho_k[i] += -a1
            for s in range(N + 1):
                rho_u[s] += slope * r[s] * a1
            val -= (-_log10(lo) + slope * lo + self.c_node[i]) * a1
        # tangents: kappa_i - D_i/(tp ln10) <= log10(tp) - 1/ln10 - c_i   (one block per tangent point)
        ln10 = arb(10).log()
        for tps in self.tps:
            ztan = take(N + 1)
            for i in range(N + 1):
                tp = tps[i]; a1 = _a(ztan[i]); r = Drow(i)
                rho_k[i] += a1
                for s in range(N + 1):
                    rho_u[s] += -(r[s] / (tp * ln10)) * a1
                val -= (_log10(tp) - 1 / ln10 - self.c_node[i]) * a1
        # kappa differences: kappa_{i+1} - kappa_i <= dhi_i ; -(kappa_{i+1} - kappa_i) <= -dlo_i
        diff_blocks = []
        if np.isfinite(self.fr.spec.L):
            diff_blocks.append(self._kappa_diff_bounds())
        diff_blocks.append(self._bracket_diff_bounds())
        for dlo, dhi in diff_blocks:
            zdh = take(N); zdl = take(N)
            for i in range(N):
                a1, a2 = _a(zdh[i]), _a(zdl[i])
                rho_k[i + 1] += a1 - a2
                rho_k[i] += -a1 + a2
                val -= dhi[i] * a1 + (-dlo[i]) * a2
        # Mp box: -Mp <= -Mp_lo ; Mp <= Mp_hi
        zml = take(1); zmh = take(1)
        rho_M += -_a(zml[0]) + _a(zmh[0])
        val -= (-self.Mp_lo) * _a(zml[0]) + self.Mp_hi * _a(zmh[0])
        assert p == len(z_in), (p, len(z_in))
        # SOC rows: s = (sqrt(T), w_b, w_s) = b - A y with A rows: 0 ; -e_wb ; -e_ws
        z0 = _a(z_soc[0])
        val -= self.T.sqrt() * z0
        for r in range(nb_):
            rho_wb[r] += -_a(z_soc[1 + r])
        for j in range(n):
            rho_ws[j] += -_a(z_soc[1 + nb_ + j])
        # objective q on u (and optionally on kappa)
        for i in range(N + 1):
            rho_u[i] += _a(q_u[i])
            if q_k is not None:
                rho_k[i] += _a(q_k[i])
        self._rho = dict(u=rho_u, k=rho_k, M=rho_M, l=rho_l, wb=rho_wb, ws=rho_ws)
        # ---- residual absorption with variable boxes ----
        Tsq = self.T.sqrt()
        loss = arb(0)
        Y_u = [_max(abs(self.ulo), abs(self.uhi)) for _ in range(N + 1)]
        Y_k = []
        for i in range(N + 1):
            klo = _log10(self.lo[i]) - self.c_node[i]; khi = _log10(self.hi[i]) - self.c_node[i]
            Y_k.append(_max(abs(klo), abs(khi)))
        for i in range(N + 1):
            loss += abs(rho_u[i]) * Y_u[i] + abs(rho_k[i]) * Y_k[i]
        loss += abs(rho_M) * (abs(self.Mp_lo) + abs(self.Mp_hi))
        for j in range(n):
            k = self.k_sn[j]
            Yl = Y_k[k] + Y_k[k + 1] + abs(self.c_sn[j]) + self.e_seg[k]
            loss += abs(rho_l[j]) * Yl + abs(rho_ws[j]) * Tsq
        for r in range(nb_):
            loss += abs(rho_wb[r]) * Tsq
        bound = val - loss
        lb = _endpoint(bound, -1)
        self._val, self._loss = val, loss
        if verbose:
            print(f"    certificate: -b^T z = {val.str(12)}  residual loss = {loss.str(5)}  "
                  f"=> rigorous lower bound {lb:.8f}")
        return lb


def rigorous_chi2(fr: Frozen, u: np.ndarray, Mp: float) -> arb:
    """Enclosure of chi2(u, Mp) with exact log10 (for the reference point / T)."""
    v = Verifier.__new__(Verifier)  # lightweight geometry only
    # build the pieces we need without the full constructor
    sp = fr.spec
    x = [_a(t) for t in sp.x]
    ex = [t.exp() for t in x]
    hs = [x[k + 1] - x[k] for k in range(sp.n_seg)]
    ua = [_a(t) for t in u]

    def dm(xv):
        mid = float(xv.mid().str(20, radius=False))
        k = int(np.clip(np.searchsorted(sp.x, mid, side="right") - 1, 0, sp.n_seg - 1))
        s = arb(0)
        for i in range(k):
            h = hs[i]; ea, eb = ex[i], ex[i + 1]
            s += ua[i] * (eb - (h + 1) * ea) / h + ua[i + 1] * ((h - 1) * eb + ea) / h
        h = hs[k]; xa, xb = x[k], x[k + 1]; ea = ex[k]; exv = xv.exp()
        s += ua[k] * ((xb - xv + 1) * exv - (h + 1) * ea) / h + ua[k + 1] * ((xv - xa - 1) * exv + ea) / h
        return s

    def dh(xv):
        mid = float(xv.mid().str(20, radius=False))
        k = int(np.clip(np.searchsorted(sp.x, mid, side="right") - 1, 0, sp.n_seg - 1))
        t = (xv - x[k]) / hs[k]
        return ua[k] * (1 - t) + ua[k + 1] * t

    # BAO
    rb = []
    for z, val, kind in zip(fr.bao.z, fr.bao.value, fr.bao.kind):
        xv = (_a(1.0) + _a(z)).log()
        pred = dm(xv) if kind == "DM_over_rs" else dh(xv)
        rb.append(_a(val) - pred)
    Wb = arb_mat([[_a(t) for t in row] for row in fr.Wb])
    wb = Wb * arb_mat([[r] for r in rb])
    chi2 = sum((wb[i, 0] * wb[i, 0] for i in range(len(rb))), arb(0))
    # SN
    rs = []
    for j in range(len(fr.sn.m)):
        xv = (_a(1.0) + _a(fr.sn.zHD[j])).log()
        D = dm(xv)
        ell = _log10(D) + 5 * _log10(_a(1.0) + _a(fr.sn.zHEL[j])) / 5  # log10 D + log10(1+zHEL)
        rs.append(_a(fr.sn.m[j]) - 5 * ell - _a(Mp))
    Wsn = arb_mat([[_a(t) for t in row] for row in fr.Wsn])
    ws = Wsn * arb_mat([[r] for r in rs])
    chi2 += sum((ws[j, 0] * ws[j, 0] for j in range(len(rs))), arb(0))
    return chi2
