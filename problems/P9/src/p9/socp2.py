"""Direct Clarabel formulation, kappa (log-distance-node) relaxation — FORMULATION v2 §3.

Variables y = [u (N+1) | kappa (N+1) | Mp | ell (n) | w_b (n_b) | w_s (n)].

- u: node values of c/(r_d H) (the class lives here; BAO predictions are exact linear maps of u).
- kappa_i = log10[ D_M(x_i)/(e^{x_i}-1) ] for i >= 1, kappa_0 = log10 u_0.
  Linked to u by a secant/tangent sandwich on certified brackets of D_i = (A_nodes u)_i (and u_0).
- ell_j (per SN) = kappa(x_j) + log10(e^{x_j}-1) with kappa(x_j) the linear interpolation of the
  node kappas, up to a rigorous class-only interpolation error e_k (geometry.kappa_interp_slack).
- chi2_BAO(u) + chi2_SN(ell, Mp) <= T as a second-order cone.

Every relaxation step encloses the true feasible set, so min u_0 over this set is a valid lower
bound on u_0 over F, i.e. an upper bound on H0.
"""

from __future__ import annotations

from dataclasses import dataclass

import clarabel
import numpy as np
import scipy.sparse as sp

from .geometry import kappa_difference_bounds, kappa_interp_slack, segment_index
from .model import Frozen

LN10 = np.log(10.0)
MP_BOX = (0.0, 40.0)   # box on the SN nuisance M' (part of the definition of F)


@dataclass
class NodeBounds:
    """Certified bounds on u_i (i=0..N) and D_i = D_M(x_i)/r_d (i=1..N; D_0 = 0)."""
    u_lo: np.ndarray
    u_hi: np.ndarray
    dm_lo: np.ndarray
    dm_hi: np.ndarray

    def width(self) -> float:
        w_u = np.max(np.log10(self.u_hi / self.u_lo))
        w_d = np.max(np.log10(self.dm_hi[1:] / self.dm_lo[1:]))
        return float(max(w_u, w_d))

    @staticmethod
    def a_priori(fr: Frozen) -> "NodeBounds":
        ulo, uhi = fr.spec.u_box
        N = fr.spec.n_seg
        u_lo = np.full(N + 1, ulo); u_hi = np.full(N + 1, uhi)
        return NodeBounds(u_lo, u_hi, fr.A_nodes @ u_lo, fr.A_nodes @ u_hi)


N_TANGENTS = 3


def tangent_points(lo: np.ndarray, hi: np.ndarray, tangent_u, fr: Frozen) -> list[np.ndarray]:
    """K tangent points per node inside [lo, hi]: geometric quantiles, plus the reference point if given."""
    K = N_TANGENTS
    pts = [lo ** (1 - s) * hi ** s for s in np.linspace(0.15, 0.85, K)]
    if tangent_u is not None:
        tp = np.concatenate([[tangent_u[0]], (fr.A_nodes @ tangent_u)[1:]])
        pts[K // 2] = np.clip(tp, lo, hi)
    return pts


def full_segment_coeffs(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """inc_k = alpha_k u_k + beta_k u_{k+1} for the full segment k (both coefficients >= 0)."""
    hs = np.diff(x); ex = np.exp(x)
    alpha = (ex[1:] - (hs + 1.0) * ex[:-1]) / hs
    beta = ((hs - 1.0) * ex[1:] + ex[:-1]) / hs
    return alpha, beta


def bracket_kappa_difference_bounds(fr: Frozen, nb: NodeBounds) -> tuple[np.ndarray, np.ndarray]:
    """Bounds on kappa_{i+1} - kappa_i from the current node bounds (valid on F given nb valid).

    i >= 1: D_{i+1}/D_i = 1 + inc_i/D_i, inc_i in [a u_lo,i + b u_lo,i+1, a u_hi,i + b u_hi,i+1], D_i in [dm_lo, dm_hi].
    i = 0:  D_1/u_0 = a_0 + b_0 u_1/u_0 in [a_0 + b_0 u_lo,1/u_hi,0, a_0 + b_0 u_hi,1/u_lo,0].
    """
    x = fr.spec.x
    N = len(x) - 1
    a, b = full_segment_coeffs(x)
    c = np.concatenate([[0.0], np.log10(np.expm1(x[1:]))])
    lo = np.empty(N); hi = np.empty(N)
    for i in range(N):
        inc_lo = a[i] * nb.u_lo[i] + b[i] * nb.u_lo[i + 1]
        inc_hi = a[i] * nb.u_hi[i] + b[i] * nb.u_hi[i + 1]
        if i == 0:
            r_lo = a[0] + b[0] * nb.u_lo[1] / nb.u_hi[0]
            r_hi = a[0] + b[0] * nb.u_hi[1] / nb.u_lo[0]
        else:
            r_lo = 1.0 + inc_lo / nb.dm_hi[i]
            r_hi = 1.0 + inc_hi / nb.dm_lo[i]
        lo[i] = np.log10(r_lo) - (c[i + 1] - c[i])
        hi[i] = np.log10(r_hi) - (c[i + 1] - c[i])
    return lo, hi


class KappaModel:
    MP_BOX = MP_BOX

    def __init__(self, fr: Frozen, nb: NodeBounds, T: float | None, tangent_u: np.ndarray | None = None,
                 tol: float = 1e-8):
        self.fr, self.nb, self.T = fr, nb, T
        spc = fr.spec
        assert not spc.use_dv, "this model has no D_V row (use the LKR relaxation)"
        x = spc.x
        N, n, nbao = spc.n_seg, len(fr.sn.m), len(fr.bao.value)
        self.N, self.n = N, n
        iu = np.arange(N + 1)
        ik = N + 1 + np.arange(N + 1)
        iM = 2 * N + 2
        il = 2 * N + 3 + np.arange(n)
        iwb = 2 * N + 3 + n + np.arange(nbao)
        iws = 2 * N + 3 + n + nbao + np.arange(n)
        nvar = 2 * N + 3 + 2 * n + nbao
        self.idx = dict(u=iu, kappa=ik, Mp=iM, ell=il, wb=iwb, ws=iws)
        self.nvar = nvar

        def E(rows, cols_vals):
            """helper: build a sparse block from (row, col, val) triplets"""
            r, c, v = zip(*cols_vals) if cols_vals else ((), (), ())
            return sp.csr_matrix((np.array(v, float), (np.array(r), np.array(c))), shape=(rows, nvar))

        A_blocks, b_blocks = [], []
        # ---------------- equalities ----------------
        Eb = sp.lil_matrix((nbao, nvar)); Eb[:, iu] = fr.Wb @ fr.P; Eb[:, iwb] = np.eye(nbao)
        A_blocks.append(Eb.tocsr()); b_blocks.append(fr.Wb @ fr.bao.value)
        Es = sp.lil_matrix((n, nvar)); Es[:, il] = 5.0 * fr.Wsn
        Es[:, iM] = (fr.Wsn @ np.ones(n)).reshape(-1, 1); Es[:, iws] = np.eye(n)
        A_blocks.append(Es.tocsr()); b_blocks.append(fr.Wsn @ (fr.sn.m - fr.sn_offset))
        n_eq = nbao + n

        # ---------------- inequalities G y <= h ----------------
        G, h = [], []
        ulo, uhi = spc.u_box
        Iu = sp.lil_matrix((N + 1, nvar)); Iu[:, iu] = np.eye(N + 1); Iu = Iu.tocsr()
        G += [-Iu, Iu]; h += [-np.full(N + 1, ulo), np.full(N + 1, uhi)]
        if np.isfinite(spc.L):
            hs = spc.hs
            Dm = E(N, [(i, iu[i + 1], 1.0) for i in range(N)] + [(i, iu[i], -1.0) for i in range(N)])
            Sl = E(N, [(i, iu[i], spc.L * hs[i]) for i in range(N)])
            Sr = E(N, [(i, iu[i + 1], spc.L * hs[i]) for i in range(N)])
            G += [Dm - Sl, Dm - Sr, -Dm - Sl, -Dm - Sr]; h += [np.zeros(N)] * 4
        # SN interpolation: ell_j - (1-t) kappa_k - t kappa_{k+1} in [c_j - e_k, c_j + e_k]
        xj = np.log1p(fr.sn.zHD)
        k = segment_index(x, xj)
        t = (xj - x[k]) / spc.hs[k]
        e_seg = kappa_interp_slack(x, spc.L) if np.isfinite(spc.L) else np.full(N, np.inf)
        self.e_seg = e_seg
        cj = np.log10(np.expm1(xj))
        ej = e_seg[k]
        if not np.all(np.isfinite(ej)):
            raise ValueError("infinite interpolation slack (L = inf not supported by the kappa model)")
        Il = E(n, [(j, il[j], 1.0) for j in range(n)])
        Ik = E(n, [(j, ik[k[j]], 1.0 - t[j]) for j in range(n)] + [(j, ik[k[j] + 1], t[j]) for j in range(n)])
        G += [Il - Ik, -Il + Ik]; h += [cj + ej, -cj + ej]
        # node links: y_i := kappa_i + c_i = log10 D_i, sandwiched on bracket [lo_i, hi_i] of D_i = A_i u
        # i = 0: D_0 := u_0 (c_0 = 0); i >= 1: D_i = A_nodes[i] u, c_i = log10(e^{x_i} - 1)
        lo = np.concatenate([[nb.u_lo[0]], nb.dm_lo[1:]])
        hi = np.concatenate([[nb.u_hi[0]], nb.dm_hi[1:]])
        cnode = np.concatenate([[0.0], np.log10(np.expm1(x[1:]))])
        tps = tangent_points(lo, hi, tangent_u, fr)
        self.tps = tps
        Dlin = sp.lil_matrix((N + 1, nvar))
        Dlin[0, iu[0]] = 1.0
        Dlin[1:, iu] = fr.A_nodes[1:]
        Dlin = Dlin.tocsr()
        Ikn = E(N + 1, [(i, ik[i], 1.0) for i in range(N + 1)])
        slope = (np.log10(hi) - np.log10(lo)) / (hi - lo)
        G += [-Dlin, Dlin]; h += [-lo, hi]
        G.append(-Ikn + sp.diags(slope) @ Dlin); h.append(-np.log10(lo) + slope * lo + cnode)     # secant
        for tp in tps:                                                                            # tangents
            G.append(Ikn - sp.diags(1.0 / (tp * LN10)) @ Dlin); h.append(np.log10(tp) - 1.0 / LN10 - cnode)
        Dk = E(N, [(i, ik[i + 1], 1.0) for i in range(N)] + [(i, ik[i], -1.0) for i in range(N)])
        if np.isfinite(spc.L):
            # class-implied smoothness of kappa: dlo_i <= kappa_{i+1} - kappa_i <= dhi_i (bracket-free)
            dlo, dhi = kappa_difference_bounds(x, spc.L)
            G += [Dk, -Dk]; h += [dhi, -dlo]
            self.kappa_diff_bounds = (dlo, dhi)
        # bracket-aware smoothness of kappa (valid given the current node bounds)
        blo, bhi = bracket_kappa_difference_bounds(fr, nb)
        G += [Dk, -Dk]; h += [bhi, -blo]
        self.bracket_diff_bounds = (blo, bhi)
        # Mp box (part of the definition of F)
        IM = E(1, [(0, iM, 1.0)])
        G += [-IM, IM]; h += [np.array([-self.MP_BOX[0]]), np.array([self.MP_BOX[1]])]
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
        self.Dlin = Dlin

    def _solve(self, P, q):
        solver = clarabel.DefaultSolver(P, q, self.A, self.b, self.cones, self.settings)
        sol = solver.solve()
        st = str(sol.status)
        if "Solved" not in st:
            raise RuntimeError(f"clarabel status {st}")
        return sol

    def min_chi2(self):
        sol = self._solve(self.Pw, np.zeros(self.nvar))
        y = np.asarray(sol.x)
        return 2.0 * sol.obj_val, y[self.idx["u"]].copy(), float(y[self.idx["Mp"]])

    def extremize(self, c_u: np.ndarray, c_k: np.ndarray | None = None) -> float:
        q = np.zeros(self.nvar); q[self.idx["u"]] = c_u
        if c_k is not None:
            q[self.idx["kappa"]] = c_k
        return float(self._solve(self.P0, q).obj_val)

    def node_bounds(self, verbose: bool = False) -> NodeBounds:
        """Min/max of u_i, D_i (via u) and kappa_i over the relaxed set; D brackets are the
        intersection of the u-side bounds with 10^(kappa bounds + c_i)."""
        N = self.N
        A = self.fr.A_nodes
        x = self.fr.spec.x
        c = np.concatenate([[0.0], np.log10(np.expm1(x[1:]))])
        z = np.zeros(N + 1)
        u_lo = np.empty(N + 1); u_hi = np.empty(N + 1); dm_lo = np.zeros(N + 1); dm_hi = np.zeros(N + 1)
        for i in range(N + 1):
            e = np.zeros(N + 1); e[i] = 1.0
            u_lo[i] = self.extremize(e); u_hi[i] = -self.extremize(-e)
            k_lo = self.extremize(z, e); k_hi = -self.extremize(z, -e)
            if i > 0:
                dm_lo[i] = max(self.extremize(A[i]), 10 ** (k_lo + c[i]))
                dm_hi[i] = min(-self.extremize(-A[i]), 10 ** (k_hi + c[i]))
            else:
                u_lo[0] = max(u_lo[0], 10 ** k_lo); u_hi[0] = min(u_hi[0], 10 ** k_hi)
            if verbose and i % 15 == 0:
                print(f"    node {i:3d} x={x[i]:.4f}: u in [{u_lo[i]:.4f},{u_hi[i]:.4f}]"
                      f"  D in [{dm_lo[i]:.5f},{dm_hi[i]:.5f}]", flush=True)
        # intersect with previous bounds (both valid)
        return NodeBounds(np.maximum(u_lo, self.nb.u_lo), np.minimum(u_hi, self.nb.u_hi),
                          np.maximum(dm_lo, self.nb.dm_lo), np.minimum(dm_hi, self.nb.dm_hi))
