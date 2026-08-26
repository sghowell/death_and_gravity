"""Direct Clarabel formulation of the relaxed problem (FORMULATION.md §3.3).

Same mathematics as model.Model (cvxpy) but assembled as sparse conic data, so that
hundreds of objective changes cost one factorization each. The cvxpy path remains the
independent implementation for cross-checks.

Variable layout y = [u (N+1) | Mp (1) | ell (n) | w_b (n_b) | w_s (n)].
Cones (in row order): zero (equalities), nonnegative (inequalities), second-order (chi2 <= T).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import clarabel

from . import C_KM_S
from .model import Brackets, Frozen

LN10 = np.log(10.0)


class DirectModel:
    def __init__(self, fr: Frozen, br: Brackets, T: float | None, tol: float = 1e-9):
        self.fr, self.br, self.T = fr, br, T
        spc = fr.spec
        N, n, nb = spc.n_seg, len(fr.sn.m), len(fr.bao.value)
        self.N, self.n, self.nb = N, n, nb
        iu = np.arange(N + 1)
        iM = N + 1
        il = N + 2 + np.arange(n)
        iwb = N + 2 + n + np.arange(nb)
        iws = N + 2 + n + nb + np.arange(n)
        nvar = N + 2 + 2 * n + nb
        self.idx = dict(u=iu, Mp=iM, ell=il, wb=iwb, ws=iws)
        self.nvar = nvar

        rows_A, rows_b = [], []

        # ---- equalities (zero cone) ----
        # w_b + Wb P u = Wb d
        Eb = sp.lil_matrix((nb, nvar))
        Eb[:, iu] = fr.Wb @ fr.P
        Eb[:, iwb] = np.eye(nb)
        rows_A.append(Eb.tocsr()); rows_b.append(fr.Wb @ fr.bao.value)
        # w_s + 5 Wsn ell + Wsn 1 Mp = Wsn (m - off)
        Es = sp.lil_matrix((n, nvar))
        Es[:, il] = 5.0 * fr.Wsn
        Es[:, iM] = (fr.Wsn @ np.ones(n)).reshape(-1, 1)
        Es[:, iws] = np.eye(n)
        rows_A.append(Es.tocsr()); rows_b.append(fr.Wsn @ (fr.sn.m - fr.sn_offset))
        n_eq = nb + n

        # ---- inequalities (nonnegative cone): G y <= h ----
        G_rows, h_rows = [], []
        ulo, uhi = spc.u_box
        I_u = sp.lil_matrix((N + 1, nvar)); I_u[:, iu] = np.eye(N + 1); I_u = I_u.tocsr()
        G_rows += [-I_u, I_u]; h_rows += [-np.full(N + 1, ulo), np.full(N + 1, uhi)]
        if np.isfinite(spc.L):
            Lh = spc.L * spc.hs
            Dm = sp.lil_matrix((N, nvar))                       # d_i = u_{i+1} - u_i
            for i in range(N):
                Dm[i, iu[i + 1]] = 1.0; Dm[i, iu[i]] = -1.0
            Dm = Dm.tocsr()
            Sl = sp.lil_matrix((N, nvar)); Sr = sp.lil_matrix((N, nvar))
            for i in range(N):
                Sl[i, iu[i]] = Lh[i]; Sr[i, iu[i + 1]] = Lh[i]
            Sl, Sr = Sl.tocsr(), Sr.tocsr()
            G_rows += [Dm - Sl, Dm - Sr, -Dm - Sl, -Dm - Sr]
            h_rows += [np.zeros(N)] * 4
        # brackets lo <= D <= hi with D = A_sn u
        Asn = sp.lil_matrix((n, nvar)); Asn[:, iu] = fr.A_sn; Asn = Asn.tocsr()
        G_rows += [-Asn, Asn]; h_rows += [-br.lo, br.hi]
        # sandwich
        lo, hi, t = br.lo, br.hi, br.t
        slope = (np.log10(hi) - np.log10(lo)) / (hi - lo)
        I_l = sp.lil_matrix((n, nvar)); I_l[:, il] = np.eye(n); I_l = I_l.tocsr()
        Sd = sp.diags(slope) @ Asn
        G_rows.append(-I_l + Sd); h_rows.append(-np.log10(lo) + slope * lo)        # secant
        Td = sp.diags(1.0 / (t * LN10)) @ Asn
        G_rows.append(I_l - Td); h_rows.append(np.log10(t) - 1.0 / LN10)           # tangent
        G = sp.vstack(G_rows).tocsr(); h = np.concatenate(h_rows)
        n_in = G.shape[0]
        rows_A.append(G); rows_b.append(h)

        cones = [clarabel.ZeroConeT(n_eq), clarabel.NonnegativeConeT(n_in)]
        # ---- SOC: (sqrt(T), w_b, w_s) ----
        if T is not None:
            nsoc = 1 + nb + n
            S = sp.lil_matrix((nsoc, nvar))
            S[1:1 + nb, iwb] = -np.eye(nb)
            S[1 + nb:, iws] = -np.eye(n)
            rows_A.append(S.tocsr()); rows_b.append(np.concatenate([[np.sqrt(T)], np.zeros(nb + n)]))
            cones.append(clarabel.SecondOrderConeT(nsoc))
        self.A = sp.vstack(rows_A).tocsc()
        self.b = np.concatenate(rows_b)
        self.cones = cones
        self.settings = clarabel.DefaultSettings()
        self.settings.verbose = False
        self.settings.tol_gap_abs = tol; self.settings.tol_gap_rel = tol
        self.settings.tol_feas = tol
        self.settings.max_iter = 200
        self.P0 = sp.csc_matrix((nvar, nvar))
        dvec = np.zeros(nvar); dvec[iwb] = 1.0; dvec[iws] = 1.0
        self.Pw = sp.diags(dvec).tocsc()

    def _solve(self, P, q):
        solver = clarabel.DefaultSolver(P, q, self.A, self.b, self.cones, self.settings)
        sol = solver.solve()
        st = str(sol.status)
        if "Solved" not in st:
            raise RuntimeError(f"clarabel status {st}")
        return sol

    def min_chi2(self):
        """min chi2 (relaxed) over the constraints (SOC ignored if present). Returns (chi2, u, Mp)."""
        sol = self._solve(self.Pw, np.zeros(self.nvar))
        y = np.asarray(sol.x)
        return 2.0 * sol.obj_val, y[self.idx["u"]].copy(), float(y[self.idx["Mp"]])

    def extremize(self, c_u: np.ndarray) -> tuple[float, np.ndarray]:
        """min c_u @ u over the relaxed feasible set; returns (value, y)."""
        q = np.zeros(self.nvar); q[self.idx["u"]] = c_u
        sol = self._solve(self.P0, q)
        return float(sol.obj_val), np.asarray(sol.x)

    def node_brackets(self, verbose=False):
        N = self.N
        u_lo = np.empty(N + 1); u_hi = np.empty(N + 1); dm_lo = np.empty(N + 1); dm_hi = np.empty(N + 1)
        A = self.fr.A_nodes
        for i in range(N + 1):
            e = np.zeros(N + 1); e[i] = 1.0
            u_lo[i] = self.extremize(e)[0]
            u_hi[i] = -self.extremize(-e)[0]
            if i == 0:
                dm_lo[i] = dm_hi[i] = 0.0
            else:
                dm_lo[i] = self.extremize(A[i])[0]
                dm_hi[i] = -self.extremize(-A[i])[0]
            if verbose and i % 10 == 0:
                print(f"    node {i}: u in [{u_lo[i]:.4f},{u_hi[i]:.4f}]  DM in [{dm_lo[i]:.4f},{dm_hi[i]:.4f}]")
        return u_lo, u_hi, dm_lo, dm_hi


def a_priori_brackets(fr: Frozen) -> Brackets:
    """Brackets from the class box alone (always valid; very loose)."""
    spc = fr.spec
    N = spc.n_seg
    ulo, uhi = spc.u_box
    u_lo = np.full(N + 1, ulo); u_hi = np.full(N + 1, uhi)
    dm_lo = fr.A_nodes @ u_lo; dm_hi = fr.A_nodes @ u_hi
    return Brackets.from_nodes(fr, u_lo, u_hi, dm_lo, dm_hi)
