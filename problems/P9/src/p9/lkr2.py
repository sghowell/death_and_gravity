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
from .lkr_rows import FloatArith, build, obbt_objectives
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
        """Solve; on numerical trouble retry at looser tolerances (the verifier certifies any returned
        dual, so a looser solve only weakens the bound slightly)."""
        last = None
        for tol in (self.settings.tol_gap_abs, 1e-8, 1e-7, 1e-6):
            st = clarabel.DefaultSettings(); st.verbose = False
            st.tol_gap_abs = tol; st.tol_gap_rel = tol; st.tol_feas = tol; st.max_iter = 500
            solver = clarabel.DefaultSolver(P, q, self.A, self.b, self.cones, st)
            sol = solver.solve()
            last = str(sol.status)
            if "Solved" in last:
                return sol
        raise RuntimeError(f"clarabel status {last}")

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
        """OBBT over lkr_rows.obbt_objectives (rho at rho_nodes, lambda at lam_nodes, all y_b, all y_V)."""
        br = self.br.copy()
        for kind, i, side, qd in obbt_objectives(self.lay, lam_nodes):
            if kind == "rho" and rho_nodes is not None and i not in rho_nodes:
                continue
            q = np.zeros(self.nvar)
            for v, cf in qd.items():
                q[v] = cf
            br.apply_bound(kind, i, side, self.extremize(q), float(self.c_node[i]) if kind == "rho" else 0.0)
        return br
