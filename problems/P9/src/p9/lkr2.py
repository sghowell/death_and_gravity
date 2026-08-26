"""LKR model assembled from the shared row builder (lkr_rows.build with FloatArith).

Same public interface as lkr.LKRModel (idx, nvar, A, b, cones, P0, Pw, _n_eq, _n_in, enodes, idx_dm,
c_node, min_chi2, extremize, min_lambda0, tighten), so parallel.tighten_parallel works unchanged.
Row order: eq = [BAO (nbao) | SN dense block (n) | kappa0=lambda0 | segment identities (N) | D_H rows]
           le = as produced by lkr_rows.build; soc = (sqrt T, w_b, w_s).
"""

from __future__ import annotations

import time

import clarabel
import numpy as np
import scipy.sparse as sp

from .lkr import Brackets3
from .lkr_rows import FloatArith, build
from .model import Frozen


def _rows_to_csr(rows, nvar):
    r_idx, c_idx, vals, rhs = [], [], [], []
    for r, (coeffs, b) in enumerate(rows):
        for v, cf in coeffs:
            r_idx.append(r); c_idx.append(int(v)); vals.append(float(cf))
        rhs.append(float(b))
    A = sp.csr_matrix((vals, (r_idx, c_idx)), shape=(len(rows), nvar))
    return A, np.array(rhs)


class LKRModel2:
    def __init__(self, fr: Frozen, br: Brackets3, T: float | None, tol: float = 1e-10):
        self.fr, self.br, self.T = fr, br, T
        ar = FloatArith()
        B = build(fr, br, T, ar)
        self.build = B
        lay = B.layout
        self.lay = lay
        self.N, self.n = lay.N, lay.n
        self.nvar = lay.nvar
        self.idx = dict(lam=lay.lam, kappa=lay.kappa, Mp=lay.Mp, ell=lay.ell, Ed=lay.Ed, Er=lay.Er, Es=lay.Es,
                        yb=lay.yb, P=lay.P, U=lay.U, wb=lay.wb, ws=lay.ws)
        self.enodes, self.idx_dm, self.idx_dh = lay.enodes, lay.idx_dm, lay.idx_dh
        self.c_node = np.array([float(v) for v in B.c_node])
        n, nbao = lay.n, lay.nbao
        # equalities: BAO rows (first nbao of eq_rows), then SN dense block, then the rest
        A_bao, b_bao = _rows_to_csr(B.eq_rows[:nbao], self.nvar)
        Es = sp.lil_matrix((n, self.nvar)); Es[:, lay.ell] = 5.0 * fr.Wsn
        Es[:, lay.Mp] = (fr.Wsn @ np.ones(n)).reshape(-1, 1); Es[:, lay.ws] = np.eye(n)
        b_sn = fr.Wsn @ np.array([float(v) for v in B.sn_rhs])
        A_rest, b_rest = _rows_to_csr(B.eq_rows[nbao:], self.nvar)
        A_eq = sp.vstack([A_bao, Es.tocsr(), A_rest]).tocsr(); b_eq = np.concatenate([b_bao, b_sn, b_rest])
        A_le, b_le = _rows_to_csr(B.le_rows, self.nvar)
        blocks = [A_eq, A_le]; rhs = [b_eq, b_le]
        cones = [clarabel.ZeroConeT(A_eq.shape[0]), clarabel.NonnegativeConeT(A_le.shape[0])]
        self._n_eq, self._n_in = A_eq.shape[0], A_le.shape[0]
        if T is not None:
            nsoc = 1 + nbao + n
            S = sp.lil_matrix((nsoc, self.nvar)); S[1:1 + nbao, lay.wb] = -np.eye(nbao); S[1 + nbao:, lay.ws] = -np.eye(n)
            blocks.append(S.tocsr()); rhs.append(np.concatenate([[float(B.soc_T)], np.zeros(nbao + n)]))
            cones.append(clarabel.SecondOrderConeT(nsoc))
        self.A = sp.vstack(blocks).tocsc(); self.b = np.concatenate(rhs); self.cones = cones
        self.settings = clarabel.DefaultSettings(); self.settings.verbose = False
        self.settings.tol_gap_abs = tol; self.settings.tol_gap_rel = tol; self.settings.tol_feas = tol
        self.settings.max_iter = 400
        self.P0 = sp.csc_matrix((self.nvar, self.nvar))
        dvec = np.zeros(self.nvar); dvec[lay.wb] = 1.0; dvec[lay.ws] = 1.0
        self.Pw = sp.diags(dvec).tocsc()

    def _solve(self, P, q):
        solver = clarabel.DefaultSolver(P, q, self.A, self.b, self.cones, self.settings)
        sol = solver.solve()
        if "Solved" not in str(sol.status):
            raise RuntimeError(f"clarabel status {sol.status}")
        return sol

    def solve_dual(self, q):
        """Returns (primal value, x, z) for min q.y."""
        sol = self._solve(self.P0, q)
        return float(sol.obj_val), np.asarray(sol.x), np.asarray(sol.z)

    def min_chi2(self):
        sol = self._solve(self.Pw, np.zeros(self.nvar)); y = np.asarray(sol.x)
        return 2.0 * sol.obj_val, y

    def extremize(self, q):
        return float(self._solve(self.P0, q).obj_val)

    def min_lambda0(self) -> float:
        q = np.zeros(self.nvar); q[self.idx["lam"][0]] = 1.0
        return self.extremize(q)

    def tighten(self, verbose=False, lam_nodes=None, rho_nodes=None) -> Brackets3:
        br = self.br; N = self.N
        rho_lo = br.rho_lo.copy(); rho_hi = br.rho_hi.copy()
        lam_lo = br.lam_lo.copy(); lam_hi = br.lam_hi.copy()
        yb_lo = br.yb_lo.copy(); yb_hi = br.yb_hi.copy()
        if rho_nodes is None:
            rho_nodes = range(1, N + 1)
        if lam_nodes is None:
            lam_nodes = sorted(set(self.enodes) | {0})
        for i in rho_nodes:
            q = np.zeros(self.nvar); q[self.idx["lam"][i]] = 1.0; q[self.idx["kappa"][i]] = -1.0
            rho_lo[i] = max(rho_lo[i], self.extremize(q) - self.c_node[i])
            rho_hi[i] = min(rho_hi[i], -self.extremize(-q) - self.c_node[i])
        for i in lam_nodes:
            q = np.zeros(self.nvar); q[self.idx["lam"][i]] = 1.0
            lam_lo[i] = max(lam_lo[i], self.extremize(q)); lam_hi[i] = min(lam_hi[i], -self.extremize(-q))
        for p in range(len(self.idx_dm)):
            q = np.zeros(self.nvar); q[self.idx["yb"][p]] = 1.0
            yb_lo[p] = max(yb_lo[p], self.extremize(q)); yb_hi[p] = min(yb_hi[p], -self.extremize(-q))
        return Brackets3(rho_lo, rho_hi, yb_lo, yb_hi, lam_lo, lam_hi)
