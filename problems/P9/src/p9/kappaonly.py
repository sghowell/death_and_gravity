"""Kappa-only relaxation (FORMULATION v3 §3): no u variables at all.

Variables y = [kappa (N+1) | Mp | ell (n) | yb (n_bao) | Pb (n_bao) | E (2 n_bao) | w_b (n_bao) | w_s (n)].

- kappa_i = log10[ D_M(x_i)/(e^{x_i}-1) ] (kappa_0 = log10 u_0). Objective: min kappa_0.
- SN: ell_j = interp_j(kappa) + c_j +- e_k                     (exact convex, class-only slack)
- class: first differences  dlo_i <= kappa_{i+1}-kappa_i <= dhi_i   (class-only, geometry.kappa_difference_bounds)
         second divided differences |(k_{i+1}-k_i)/h_i - (k_i-k_{i-1})/h_{i-1}| <= (h_{i-1}+h_i) B_i/(2 ln10)
         kappa box from the class box on u.
- BAO D_M rows:  yb = interp_b(kappa) + c_b +- e ;  Pb in [tangents(yb), chord(yb)]  (10^y convex; chord needs a
  bracket [ylo, yhi] on yb, tightened by OBBT on yb itself — 6 quantities).
- BAO D_H rows:  with E_i = 10^{kappa_i + c_i} sandwiched at the two nodes around z_b, and u piecewise linear:
  D_H(z_b) = u(x_b) satisfies  (E_{k+1}-E_k)/(a_k+b_k)/(1+L h_k) <= D_H <= (E_{k+1}-E_k)/(a_k+b_k)*(1+L h_k).
- chi2 = ||Wb (d - P)||^2 + ||Wsn (m - off - 5 ell - Mp)||^2 <= T   (SOC).

Every step encloses F, so min kappa_0 over this set is a valid lower bound on log10 u_0 over F.
"""

from __future__ import annotations

from dataclasses import dataclass

import clarabel
import numpy as np
import scipy.sparse as sp

from . import C_KM_S
from .geometry import (kappa_difference_bounds, kappa_interp_slack, kappa_second_derivative_bound,
                       segment_index)
from .model import Frozen
from .socp2 import MP_BOX, full_segment_coeffs

LN10 = np.log(10.0)
N_TANGENTS = 4


@dataclass
class YBrackets:
    """Brackets on log10 D_M/r_d at the BAO D_M redshifts (yb) and on log10 D_i at the nodes used by D_H rows."""
    yb_lo: np.ndarray
    yb_hi: np.ndarray
    ynode_lo: dict   # node index -> lower bound on log10 D_i
    ynode_hi: dict


def kappa_box(fr: Frozen) -> tuple[np.ndarray, np.ndarray]:
    """Class-box implied bounds on kappa_i: D_i in [A u_lo, A u_hi] componentwise (A >= 0)."""
    spc = fr.spec
    ulo, uhi = spc.u_box
    N = spc.n_seg
    x = spc.x
    c = np.concatenate([[0.0], np.log10(np.expm1(x[1:]))])
    Dlo = np.concatenate([[ulo], (fr.A_nodes @ np.full(N + 1, ulo))[1:]])
    Dhi = np.concatenate([[uhi], (fr.A_nodes @ np.full(N + 1, uhi))[1:]])
    return np.log10(Dlo) - c, np.log10(Dhi) - c


def curvature_bounds(x: np.ndarray, L: float) -> np.ndarray:
    """For interior nodes i=1..N-1: bound on |(k_{i+1}-k_i)/h_i - (k_i-k_{i-1})/h_{i-1}| (log10 units)."""
    N = len(x) - 1
    hs = np.diff(x)
    B = np.array([kappa_second_derivative_bound(x[k], x[k + 1], L) for k in range(N)])
    out = np.empty(N - 1)
    for i in range(1, N):
        out[i - 1] = (hs[i - 1] + hs[i]) * max(B[i - 1], B[i]) / (2.0 * LN10)
    return out


class KappaOnlyModel:
    def __init__(self, fr: Frozen, yb: YBrackets, T: float | None, tol: float = 1e-8):
        self.fr, self.yb, self.T = fr, yb, T
        spc = fr.spec
        x = spc.x; hs = spc.hs
        N, n = spc.n_seg, len(fr.sn.m)
        self.N, self.n = N, n
        L = spc.L
        if not np.isfinite(L):
            raise ValueError("kappa-only model needs finite L")
        # BAO rows
        kinds = fr.bao.kind
        idx_dm = [r for r, k in enumerate(kinds) if k == "DM_over_rs"]
        idx_dh = [r for r, k in enumerate(kinds) if k == "DH_over_rs"]
        nbao = len(kinds)
        # nodes needed for D_H rows: segment k containing x_b -> nodes k, k+1
        xb = np.log1p(fr.bao.z)
        kb = segment_index(x, xb)
        enodes = sorted({int(kb[r]) for r in idx_dh} | {int(kb[r]) + 1 for r in idx_dh})
        self.enodes = enodes
        epos = {node: p for p, node in enumerate(enodes)}
        nE = len(enodes)
        # variable layout
        ik = np.arange(N + 1); iM = N + 1; il = N + 2 + np.arange(n)
        iyb = N + 2 + n + np.arange(len(idx_dm)); iP = iyb[-1] + 1 + np.arange(nbao) if len(idx_dm) else N + 2 + n + np.arange(nbao)
        iE = iP[-1] + 1 + np.arange(nE)
        iwb = iE[-1] + 1 + np.arange(nbao) if nE else iP[-1] + 1 + np.arange(nbao)
        iws = iwb[-1] + 1 + np.arange(n)
        nvar = iws[-1] + 1
        self.idx = dict(kappa=ik, Mp=iM, ell=il, yb=iyb, P=iP, E=iE, wb=iwb, ws=iws)
        self.nvar = nvar
        self.idx_dm, self.idx_dh = idx_dm, idx_dh

        def E_(rows, trip):
            r, c_, v = zip(*trip) if trip else ((), (), ())
            return sp.csr_matrix((np.array(v, float), (np.array(r), np.array(c_))), shape=(rows, nvar))

        c_node = np.concatenate([[0.0], np.log10(np.expm1(x[1:]))])
        A_blocks, b_blocks = [], []
        # ---- equalities ----
        # w_b + Wb P = Wb d
        Eb = sp.lil_matrix((nbao, nvar)); Eb[:, iP] = fr.Wb; Eb[:, iwb] = np.eye(nbao)
        A_blocks.append(Eb.tocsr()); b_blocks.append(fr.Wb @ fr.bao.value)
        # w_s + 5 Wsn ell + (Wsn 1) Mp = Wsn (m - off)
        Es = sp.lil_matrix((n, nvar)); Es[:, il] = 5.0 * fr.Wsn
        Es[:, iM] = (fr.Wsn @ np.ones(n)).reshape(-1, 1); Es[:, iws] = np.eye(n)
        A_blocks.append(Es.tocsr()); b_blocks.append(fr.Wsn @ (fr.sn.m - fr.sn_offset))
        n_eq = nbao + n
        # ---- inequalities G y <= h ----
        G, h = [], []
        klo, khi = kappa_box(fr)
        Ik = E_(N + 1, [(i, ik[i], 1.0) for i in range(N + 1)])
        G += [-Ik, Ik]; h += [-klo, khi]
        # first differences
        dlo, dhi = kappa_difference_bounds(x, L)
        Dk = E_(N, [(i, ik[i + 1], 1.0) for i in range(N)] + [(i, ik[i], -1.0) for i in range(N)])
        G += [Dk, -Dk]; h += [dhi, -dlo]
        # second divided differences
        cb = curvature_bounds(x, L)
        trip = []
        for i in range(1, N):
            trip += [(i - 1, ik[i + 1], 1.0 / hs[i]), (i - 1, ik[i], -1.0 / hs[i] - 1.0 / hs[i - 1]), (i - 1, ik[i - 1], 1.0 / hs[i - 1])]
        Ck = E_(N - 1, trip)
        G += [Ck, -Ck]; h += [cb, cb]
        # SN interpolation
        xj = np.log1p(fr.sn.zHD); kj = segment_index(x, xj); tj = (xj - x[kj]) / hs[kj]
        e_seg = kappa_interp_slack(x, L); self.e_seg = e_seg
        cj = np.log10(np.expm1(xj)); ej = e_seg[kj]
        Il = E_(n, [(j, il[j], 1.0) for j in range(n)])
        Ikj = E_(n, [(j, ik[kj[j]], 1.0 - tj[j]) for j in range(n)] + [(j, ik[kj[j] + 1], tj[j]) for j in range(n)])
        G += [Il - Ikj, -Il + Ikj]; h += [cj + ej, -cj + ej]
        # BAO D_M rows: yb - interp_b(kappa) in [c_b - e, c_b + e]; P in [tangents(yb), chord(yb)]
        self.tangents_y = []
        for p, r in enumerate(idx_dm):
            k = int(kb[r]); t = (xb[r] - x[k]) / hs[k]; cbr = np.log10(np.expm1(xb[r])); e = e_seg[k]
            row_y = E_(1, [(0, iyb[p], 1.0), (0, ik[k], -(1.0 - t)), (0, ik[k + 1], -t)])
            G += [row_y, -row_y]; h += [np.array([cbr + e]), np.array([-cbr + e])]
            ylo, yhi = yb.yb_lo[p], yb.yb_hi[p]
            # bracket rows on yb
            Iy = E_(1, [(0, iyb[p], 1.0)])
            G += [-Iy, Iy]; h += [np.array([-ylo]), np.array([yhi])]
            # tangents: P >= 10^{yt}(1 + (y - yt) ln10)  ->  -P + 10^{yt} ln10 * y <= 10^{yt}(yt ln10 - 1)
            yts = np.linspace(ylo, yhi, N_TANGENTS + 2)[1:-1]
            self.tangents_y.append(yts)
            for yt in yts:
                G.append(E_(1, [(0, iP[r], -1.0), (0, iyb[p], 10 ** yt * LN10)])); h.append(np.array([10 ** yt * (yt * LN10 - 1.0)]))
            # chord: P <= 10^{ylo} + s (y - ylo), s = (10^{yhi} - 10^{ylo})/(yhi - ylo)
            s = (10 ** yhi - 10 ** ylo) / (yhi - ylo)
            G.append(E_(1, [(0, iP[r], 1.0), (0, iyb[p], -s)])); h.append(np.array([10 ** ylo - s * ylo]))
        # E nodes: E_i = 10^{kappa_i + c_i} sandwiched on [ynode_lo, ynode_hi] (log10 D_i brackets)
        self.tangents_E = {}
        for node in enodes:
            p = epos[node]
            ylo, yhi = yb.ynode_lo[node], yb.ynode_hi[node]
            # log10 D_i = kappa_i + c_i must lie in [ylo, yhi]
            Iki = E_(1, [(0, ik[node], 1.0)])
            G += [-Iki, Iki]; h += [np.array([-(ylo - c_node[node])]), np.array([yhi - c_node[node]])]
            yts = np.linspace(ylo, yhi, N_TANGENTS + 2)[1:-1]
            self.tangents_E[node] = yts
            for yt in yts:   # E >= 10^{yt}(1 + (kappa + c - yt) ln10)
                G.append(E_(1, [(0, iE[p], -1.0), (0, ik[node], 10 ** yt * LN10)]))
                h.append(np.array([10 ** yt * (yt * LN10 - 1.0) - 10 ** yt * LN10 * c_node[node]]))
            s = (10 ** yhi - 10 ** ylo) / (yhi - ylo)
            G.append(E_(1, [(0, iE[p], 1.0), (0, ik[node], -s)])); h.append(np.array([10 ** ylo - s * ylo + s * c_node[node]]))
        # D_H rows: P_r (a+b)/(1+Lh) <= E_{k+1} - E_k ;  E_{k+1} - E_k <= P_r (a+b)(1+Lh)
        a, b = full_segment_coeffs(x)
        for r in idx_dh:
            k = int(kb[r]); ab = a[k] + b[k]; f = 1.0 + L * hs[k]
            pk, pk1 = epos[k], epos[k + 1]
            G.append(E_(1, [(0, iP[r], ab / f), (0, iE[pk1], -1.0), (0, iE[pk], 1.0)])); h.append(np.zeros(1))
            G.append(E_(1, [(0, iP[r], -ab * f), (0, iE[pk1], 1.0), (0, iE[pk], -1.0)])); h.append(np.zeros(1))
        # Mp box
        IM = E_(1, [(0, iM, 1.0)])
        G += [-IM, IM]; h += [np.array([-MP_BOX[0]]), np.array([MP_BOX[1]])]
        Gm = sp.vstack(G).tocsr(); hv = np.concatenate(h)
        A_blocks.append(Gm); b_blocks.append(hv)
        cones = [clarabel.ZeroConeT(n_eq), clarabel.NonnegativeConeT(Gm.shape[0])]
        self._n_eq, self._n_in = n_eq, Gm.shape[0]
        if T is not None:
            nsoc = 1 + nbao + n
            S = sp.lil_matrix((nsoc, nvar)); S[1:1 + nbao, iwb] = -np.eye(nbao); S[1 + nbao:, iws] = -np.eye(n)
            A_blocks.append(S.tocsr()); b_blocks.append(np.concatenate([[np.sqrt(T)], np.zeros(nbao + n)]))
            cones.append(clarabel.SecondOrderConeT(nsoc))
        self.A = sp.vstack(A_blocks).tocsc(); self.b = np.concatenate(b_blocks); self.cones = cones
        self.settings = clarabel.DefaultSettings(); self.settings.verbose = False
        self.settings.tol_gap_abs = tol; self.settings.tol_gap_rel = tol; self.settings.tol_feas = tol
        self.settings.max_iter = 300
        self.P0 = sp.csc_matrix((nvar, nvar))
        dvec = np.zeros(nvar); dvec[iwb] = 1.0; dvec[iws] = 1.0
        self.Pw = sp.diags(dvec).tocsc()
        self.c_node = c_node

    def _solve(self, P, q):
        solver = clarabel.DefaultSolver(P, q, self.A, self.b, self.cones, self.settings)
        sol = solver.solve()
        if "Solved" not in str(sol.status):
            raise RuntimeError(f"clarabel status {sol.status}")
        return sol

    def min_chi2(self):
        sol = self._solve(self.Pw, np.zeros(self.nvar))
        y = np.asarray(sol.x)
        return 2.0 * sol.obj_val, y[self.idx["kappa"]].copy(), float(y[self.idx["Mp"]])

    def extremize(self, q: np.ndarray):
        return self._solve(self.P0, q)

    def min_kappa0(self) -> float:
        q = np.zeros(self.nvar); q[self.idx["kappa"][0]] = 1.0
        return float(self.extremize(q).obj_val)

    def tighten_brackets(self) -> YBrackets:
        """OBBT on yb (BAO D_M redshifts) and on log10 D at the E nodes; intersect with current."""
        yb = self.yb
        yb_lo = yb.yb_lo.copy(); yb_hi = yb.yb_hi.copy()
        for p in range(len(self.idx_dm)):
            q = np.zeros(self.nvar); q[self.idx["yb"][p]] = 1.0
            yb_lo[p] = max(yb_lo[p], self.extremize(q).obj_val)
            yb_hi[p] = min(yb_hi[p], -self.extremize(-q).obj_val)
        ynode_lo = dict(yb.ynode_lo); ynode_hi = dict(yb.ynode_hi)
        for node in self.enodes:
            q = np.zeros(self.nvar); q[self.idx["kappa"][node]] = 1.0
            ynode_lo[node] = max(ynode_lo[node], self.extremize(q).obj_val + self.c_node[node])
            ynode_hi[node] = min(ynode_hi[node], -self.extremize(-q).obj_val + self.c_node[node])
        return YBrackets(yb_lo, yb_hi, ynode_lo, ynode_hi)


def initial_ybrackets(fr: Frozen) -> YBrackets:
    """From the class box only."""
    klo, khi = kappa_box(fr)
    spc = fr.spec; x = spc.x
    c_node = np.concatenate([[0.0], np.log10(np.expm1(x[1:]))])
    kinds = fr.bao.kind
    xb = np.log1p(fr.bao.z); kb = segment_index(x, xb)
    idx_dm = [r for r, k in enumerate(kinds) if k == "DM_over_rs"]
    idx_dh = [r for r, k in enumerate(kinds) if k == "DH_over_rs"]
    ulo, uhi = spc.u_box
    yb_lo = np.array([np.log10(fr.P[r] @ np.full(spc.n_seg + 1, ulo)) for r in idx_dm])
    yb_hi = np.array([np.log10(fr.P[r] @ np.full(spc.n_seg + 1, uhi)) for r in idx_dm])
    enodes = sorted({int(kb[r]) for r in idx_dh} | {int(kb[r]) + 1 for r in idx_dh})
    ynode_lo = {i: klo[i] + c_node[i] for i in enodes}
    ynode_hi = {i: khi[i] + c_node[i] for i in enodes}
    return YBrackets(yb_lo, yb_hi, ynode_lo, ynode_hi)


def kappa_only_bound(fr: Frozen, T: float, n_passes: int = 30, tol: float = 1e-5, verbose: bool = True):
    """Iterate OBBT on the 6+12 bracketed quantities; returns (kappa0_min, H0_max, brackets, history)."""
    yb = initial_ybrackets(fr)
    hist = []
    last = None
    for it in range(n_passes):
        m = KappaOnlyModel(fr, yb, T)
        k0 = m.min_kappa0()
        H0 = C_KM_S / (fr.spec.r_lo * 10 ** k0)
        w = float(np.max(yb.yb_hi - yb.yb_lo))
        hist.append(dict(pass_=it, kappa0_min=k0, H0_max=H0, yb_width=w))
        if verbose:
            print(f"  pass {it:2d}: kappa0_min={k0:.6f} -> H0_max={H0:.4f}  max yb width={w:.4f}", flush=True)
        if last is not None and abs(k0 - last) < tol:
            break
        last = k0
        yb = m.tighten_brackets()
    return k0, H0, yb, hist
