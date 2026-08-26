"""lambda-kappa-rho relaxation (FORMULATION v3 §3): exact class in lambda = log10 u, exact-convex SN in
kappa = log10[D_M/(e^x-1)], and the segment identities in scale-free exponential form.

For segment i (nodes i, i+1), with D_i = D_M(x_i)/r_d and increment coefficients a_i, b_i >= 0:
    D_{i+1} - D_i = a_i u_i + b_i u_{i+1}
    i >= 1:  10^{delta_i} - 1 = a_i 10^{rho_i} + b_i 10^{rho_i + dlam_i}
    i = 0 :  10^{delta_0}     = a_0            + b_0 10^{dlam_0}                       (D_0 = 0, u_0 = 10^{lambda_0})
where delta_i = kappa_{i+1} - kappa_i + c_{i+1} - c_i  (= log10 D_{i+1}/D_i), rho_i = lambda_i - kappa_i - c_i
(= log10 u_i/D_i), dlam_i = lambda_{i+1} - lambda_i, c_i = log10(e^{x_i} - 1), c_0 = 0, kappa_0 = lambda_0.
Each exponential is sandwiched (tangents below, chord above) on a bracket: delta from the class-only ratio
bounds, dlam from the class, rho from the class theta-range and then OBBT.

Variables y = [lambda (N+1) | kappa (N+1) | Mp | ell (n) | Ed (N) | Er (N) | Es (N) | yb (n_dm) | P (n_bao) | U (n_E) | w_b | w_s].
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import clarabel
import numpy as np
import scipy.sparse as sp

from . import C_KM_S
from .geometry import kappa_difference_bounds, kappa_interp_slack, segment_index, theta_bounds
from .model import Frozen
from .socp2 import MP_BOX, full_segment_coeffs

LN10 = np.log(10.0)
N_TANGENTS = 4


@dataclass
class Brackets3:
    rho_lo: np.ndarray; rho_hi: np.ndarray          # log10 u_i/D_i, i = 1..N (index 0 unused)
    yb_lo: np.ndarray; yb_hi: np.ndarray            # log10 D_M/r_d at BAO D_M redshifts
    lam_lo: np.ndarray; lam_hi: np.ndarray          # log10 u_i at all nodes (used for U sandwiches at E nodes)

    def width(self):
        return dict(rho=float(np.max(self.rho_hi[1:] - self.rho_lo[1:])), yb=float(np.max(self.yb_hi - self.yb_lo)),
                    lam=float(np.max(self.lam_hi - self.lam_lo)))


def initial_brackets3(fr: Frozen) -> Brackets3:
    spc = fr.spec; x = spc.x; N = spc.n_seg; L = spc.L
    ulo, uhi = spc.u_box
    rho_lo = np.zeros(N + 1); rho_hi = np.zeros(N + 1)
    for i in range(1, N + 1):
        tl, th = theta_bounds(x[i], L)
        em1 = np.expm1(x[i])
        rho_lo[i] = np.log10(tl / em1); rho_hi[i] = np.log10(th / em1)
    kinds = fr.bao.kind
    idx_dm = [r for r, k in enumerate(kinds) if k == "DM_over_rs"]
    yb_lo = np.array([np.log10(fr.P[r] @ np.full(N + 1, ulo)) for r in idx_dm])
    yb_hi = np.array([np.log10(fr.P[r] @ np.full(N + 1, uhi)) for r in idx_dm])
    return Brackets3(rho_lo, rho_hi, yb_lo, yb_hi, np.full(N + 1, np.log10(ulo)), np.full(N + 1, np.log10(uhi)))


def _sandwich_rows(E_, nvar, iE, lin_coeffs, lo, hi, K=N_TANGENTS):
    """Rows for  E >= 10^{y}  (tangents) and  E <= chord(y)  where y = sum lin_coeffs[var]*y_var is linear.
    Returns (G_rows, h_rows). y in [lo, hi] must be enforced separately (bracket rows)."""
    G, h = [], []
    yts = np.linspace(lo, hi, K + 2)[1:-1]
    for yt in yts:
        # E >= 10^yt (1 + (y - yt) ln10)  ->  -E + 10^yt ln10 * y <= 10^yt (yt ln10 - 1)
        trip = [(0, iE, -1.0)] + [(0, v, 10 ** yt * LN10 * cf) for v, cf in lin_coeffs]
        G.append(E_(1, trip)); h.append(np.array([10 ** yt * (yt * LN10 - 1.0)]))
    s = (10 ** hi - 10 ** lo) / (hi - lo)
    trip = [(0, iE, 1.0)] + [(0, v, -s * cf) for v, cf in lin_coeffs]
    G.append(E_(1, trip)); h.append(np.array([10 ** lo - s * lo]))
    return G, h


class LKRModel:
    def __init__(self, fr: Frozen, br: Brackets3, T: float | None, tol: float = 1e-8):
        self.fr, self.br, self.T = fr, br, T
        spc = fr.spec; x = spc.x; hs = spc.hs; L = spc.L
        N, n = spc.n_seg, len(fr.sn.m)
        self.N, self.n = N, n
        kinds = fr.bao.kind; nbao = len(kinds)
        idx_dm = [r for r, k in enumerate(kinds) if k == "DM_over_rs"]
        idx_dh = [r for r, k in enumerate(kinds) if k == "DH_over_rs"]
        xb = np.log1p(fr.bao.z); kb = segment_index(x, xb); tb = (xb - x[kb]) / hs[kb]
        enodes = sorted({int(kb[r]) for r in idx_dh} | {int(kb[r]) + 1 for r in idx_dh})
        epos = {node: p for p, node in enumerate(enodes)}
        c_node = np.concatenate([[0.0], np.log10(np.expm1(x[1:]))])
        a, b = full_segment_coeffs(x)
        # layout
        pos = 0
        def block(m):
            nonlocal pos
            r = np.arange(pos, pos + m); pos += m; return r
        il_ = block(N + 1); ik = block(N + 1); iM = block(1)[0]; il = block(n)
        iEd = block(N); iEr = block(N); iEs = block(N); iyb = block(len(idx_dm)); iP = block(nbao)
        iU = block(len(enodes)); iwb = block(nbao); iws = block(n)
        nvar = pos
        self.idx = dict(lam=il_, kappa=ik, Mp=iM, ell=il, Ed=iEd, Er=iEr, Es=iEs, yb=iyb, P=iP, U=iU, wb=iwb, ws=iws)
        self.nvar = nvar; self.enodes = enodes; self.idx_dm, self.idx_dh = idx_dm, idx_dh

        def E_(rows, trip):
            r, c_, v = zip(*trip) if trip else ((), (), ())
            return sp.csr_matrix((np.array(v, float), (np.array(r), np.array(c_))), shape=(rows, nvar))

        A_blocks, b_blocks = [], []
        # ---------------- equalities ----------------
        eq_rows, eq_b = [], []
        Eb = sp.lil_matrix((nbao, nvar)); Eb[:, iP] = fr.Wb; Eb[:, iwb] = np.eye(nbao)
        eq_rows.append(Eb.tocsr()); eq_b.append(fr.Wb @ fr.bao.value)
        Es_ = sp.lil_matrix((n, nvar)); Es_[:, il] = 5.0 * fr.Wsn
        Es_[:, iM] = (fr.Wsn @ np.ones(n)).reshape(-1, 1); Es_[:, iws] = np.eye(n)
        eq_rows.append(Es_.tocsr()); eq_b.append(fr.Wsn @ (fr.sn.m - fr.sn_offset))
        # kappa_0 = lambda_0
        eq_rows.append(E_(1, [(0, ik[0], 1.0), (0, il_[0], -1.0)])); eq_b.append(np.zeros(1))
        # segment identities: i=0: Ed_0 - b_0 Es_0 = a_0 ; i>=1: Ed_i - a_i Er_i - b_i Es_i = 1
        trip = [(0, iEd[0], 1.0), (0, iEs[0], -b[0])]
        eq_rows.append(E_(1, trip)); eq_b.append(np.array([a[0]]))
        trip = []
        for i in range(1, N):
            trip += [(i - 1, iEd[i], 1.0), (i - 1, iEr[i], -a[i]), (i - 1, iEs[i], -b[i])]
        eq_rows.append(E_(N - 1, trip)); eq_b.append(np.ones(N - 1))
        # D_H rows: P_r = (1-t) U_k + t U_{k+1}
        trip = []
        for q_, r in enumerate(idx_dh):
            k = int(kb[r]); t = tb[r]
            trip += [(q_, iP[r], 1.0), (q_, iU[epos[k]], -(1.0 - t)), (q_, iU[epos[k + 1]], -t)]
        if idx_dh:
            eq_rows.append(E_(len(idx_dh), trip)); eq_b.append(np.zeros(len(idx_dh)))
        Aeq = sp.vstack(eq_rows).tocsr(); beq = np.concatenate(eq_b)
        A_blocks.append(Aeq); b_blocks.append(beq); n_eq = Aeq.shape[0]
        # ---------------- inequalities ----------------
        G, h = [], []
        ulo, uhi = spc.u_box
        Il = E_(N + 1, [(i, il_[i], 1.0) for i in range(N + 1)])
        G += [-Il, Il]; h += [-np.maximum(np.full(N + 1, np.log10(ulo)), br.lam_lo), np.minimum(np.full(N + 1, np.log10(uhi)), br.lam_hi)]
        # class: |lambda_{i+1} - lambda_i| <= log10(1 + L h_i)
        Dl = E_(N, [(i, il_[i + 1], 1.0) for i in range(N)] + [(i, il_[i], -1.0) for i in range(N)])
        lh = np.log10(1.0 + L * hs)
        G += [Dl, -Dl]; h += [lh, lh]
        # SN interpolation
        xj = np.log1p(fr.sn.zHD); kj = segment_index(x, xj); tj = (xj - x[kj]) / hs[kj]
        e_seg = kappa_interp_slack(x, L); self.e_seg = e_seg
        cj = np.log10(np.expm1(xj)); ej = e_seg[kj]
        Iell = E_(n, [(j, il[j], 1.0) for j in range(n)])
        Ikj = E_(n, [(j, ik[kj[j]], 1.0 - tj[j]) for j in range(n)] + [(j, ik[kj[j] + 1], tj[j]) for j in range(n)])
        G += [Iell - Ikj, -Iell + Ikj]; h += [cj + ej, -cj + ej]
        # brackets and sandwiches for delta_i, rho_i, s_i
        dlo, dhi = kappa_difference_bounds(x, L)              # bounds on kappa_{i+1}-kappa_i
        dc = np.diff(c_node)
        self.delta_br = (dlo + dc, dhi + dc)
        self.s_br = (np.empty(N), np.empty(N))
        for i in range(N):
            # delta_i = kappa_{i+1} - kappa_i + dc_i ; y-lin coeffs on kappa
            lin_d = [(ik[i + 1], 1.0), (ik[i], -1.0)]
            lo_d, hi_d = dlo[i] + dc[i], dhi[i] + dc[i]
            # bracket rows: lo <= kappa_{i+1}-kappa_i + dc <= hi
            G.append(E_(1, [(0, ik[i + 1], -1.0), (0, ik[i], 1.0)])); h.append(np.array([-(lo_d - dc[i])]))
            G.append(E_(1, [(0, ik[i + 1], 1.0), (0, ik[i], -1.0)])); h.append(np.array([hi_d - dc[i]]))
            # Ed_i = 10^{delta_i}: y = (kappa_{i+1}-kappa_i) + dc -> shift handled by using y' = kappa_{i+1}-kappa_i with bracket shifted
            Gd, hd = _sandwich_rows(E_, nvar, iEd[i], lin_d, lo_d - dc[i], hi_d - dc[i])
            # _sandwich_rows models E >= 10^{y'}; we need E = 10^{y' + dc} = 10^{dc} 10^{y'}: scale rows
            for Gr, hr in zip(Gd, hd):
                Gr = Gr.tolil(); Gr[0, iEd[i]] = Gr[0, iEd[i]] / (10 ** dc[i]); G.append(Gr.tocsr()); h.append(hr)
            if i >= 1:
                # rho_i = lambda_i - kappa_i - c_i
                lin_r = [(il_[i], 1.0), (ik[i], -1.0)]
                lo_r, hi_r = br.rho_lo[i], br.rho_hi[i]
                G.append(E_(1, [(0, il_[i], -1.0), (0, ik[i], 1.0)])); h.append(np.array([-(lo_r + c_node[i])]))
                G.append(E_(1, [(0, il_[i], 1.0), (0, ik[i], -1.0)])); h.append(np.array([hi_r + c_node[i]]))
                Gr_, hr_ = _sandwich_rows(E_, nvar, iEr[i], lin_r, lo_r + c_node[i], hi_r + c_node[i])
                for Gr, hr in zip(Gr_, hr_):
                    Gr = Gr.tolil(); Gr[0, iEr[i]] = Gr[0, iEr[i]] * (10 ** c_node[i]); G.append(Gr.tocsr()); h.append(hr)
                # s_i = rho_i + dlam_i = lambda_{i+1} - kappa_i - c_i
                lin_s = [(il_[i + 1], 1.0), (ik[i], -1.0)]
                lo_s, hi_s = lo_r - lh[i], hi_r + lh[i]
            else:
                # s_0 = dlam_0 = lambda_1 - lambda_0 (rho_0 = 0)
                lin_s = [(il_[1], 1.0), (il_[0], -1.0)]
                lo_s, hi_s = -lh[0], lh[0]
            self.s_br[0][i], self.s_br[1][i] = lo_s, hi_s
            shift = c_node[i] if i >= 1 else 0.0
            G.append(E_(1, [(0, v, -cf) for v, cf in lin_s])); h.append(np.array([-(lo_s + shift)]))
            G.append(E_(1, [(0, v, cf) for v, cf in lin_s])); h.append(np.array([hi_s + shift]))
            Gs, hs_ = _sandwich_rows(E_, nvar, iEs[i], lin_s, lo_s + shift, hi_s + shift)
            for Gr, hr in zip(Gs, hs_):
                Gr = Gr.tolil(); Gr[0, iEs[i]] = Gr[0, iEs[i]] * (10 ** shift); G.append(Gr.tocsr()); h.append(hr)
        # BAO D_M rows
        for p, r in enumerate(idx_dm):
            k = int(kb[r]); t = tb[r]; cbr = np.log10(np.expm1(xb[r])); e = e_seg[k]
            row_y = E_(1, [(0, iyb[p], 1.0), (0, ik[k], -(1.0 - t)), (0, ik[k + 1], -t)])
            G += [row_y, -row_y]; h += [np.array([cbr + e]), np.array([-cbr + e])]
            ylo, yhi = br.yb_lo[p], br.yb_hi[p]
            G.append(E_(1, [(0, iyb[p], -1.0)])); h.append(np.array([-ylo]))
            G.append(E_(1, [(0, iyb[p], 1.0)])); h.append(np.array([yhi]))
            Gy, hy = _sandwich_rows(E_, nvar, iP[r], [(iyb[p], 1.0)], ylo, yhi)
            G += Gy; h += hy
        # U nodes: U_i = 10^{lambda_i} on lambda brackets
        for node in enodes:
            lo_l, hi_l = max(br.lam_lo[node], np.log10(ulo)), min(br.lam_hi[node], np.log10(uhi))
            Gu, hu = _sandwich_rows(E_, nvar, iU[epos[node]], [(il_[node], 1.0)], lo_l, hi_l)
            G += Gu; h += hu
        # Mp box
        G.append(E_(1, [(0, iM, -1.0)])); h.append(np.array([-MP_BOX[0]]))
        G.append(E_(1, [(0, iM, 1.0)])); h.append(np.array([MP_BOX[1]]))
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
        self.settings.max_iter = 400
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
        return 2.0 * sol.obj_val, y

    def extremize(self, q):
        return float(self._solve(self.P0, q).obj_val)

    def min_lambda0(self) -> float:
        q = np.zeros(self.nvar); q[self.idx["lam"][0]] = 1.0
        return self.extremize(q)

    def tighten(self, verbose=False, lam_nodes=None, rho_nodes=None) -> Brackets3:
        """OBBT on rho_i (all nodes by default), on lambda at lam_nodes (default: E nodes and node 0), and yb."""
        br = self.br; N = self.N
        rho_lo = br.rho_lo.copy(); rho_hi = br.rho_hi.copy()
        lam_lo = br.lam_lo.copy(); lam_hi = br.lam_hi.copy()
        yb_lo = br.yb_lo.copy(); yb_hi = br.yb_hi.copy()
        if rho_nodes is None:
            rho_nodes = range(1, N + 1)
        if lam_nodes is None:
            lam_nodes = sorted(set(self.enodes) | {0})
        t0 = time.time()
        for cnt, i in enumerate(rho_nodes):
            q = np.zeros(self.nvar); q[self.idx["lam"][i]] = 1.0; q[self.idx["kappa"][i]] = -1.0
            rho_lo[i] = max(rho_lo[i], self.extremize(q) - self.c_node[i])
            rho_hi[i] = min(rho_hi[i], -self.extremize(-q) - self.c_node[i])
            if verbose and cnt % 20 == 0:
                print(f"      rho node {i}: [{rho_lo[i]:.4f}, {rho_hi[i]:.4f}] [{time.time()-t0:.0f}s]", flush=True)
        for i in lam_nodes:
            q = np.zeros(self.nvar); q[self.idx["lam"][i]] = 1.0
            lam_lo[i] = max(lam_lo[i], self.extremize(q)); lam_hi[i] = min(lam_hi[i], -self.extremize(-q))
        for p in range(len(self.idx_dm)):
            q = np.zeros(self.nvar); q[self.idx["yb"][p]] = 1.0
            yb_lo[p] = max(yb_lo[p], self.extremize(q)); yb_hi[p] = min(yb_hi[p], -self.extremize(-q))
        return Brackets3(rho_lo, rho_hi, yb_lo, yb_hi, lam_lo, lam_hi)
