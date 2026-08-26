"""Convex model for the certified H0 bound (FORMULATION.md §3).

Variables: u (N+1 node values of c/(r_d H)), Mp (SN nuisance), ell (per-SN log10 D_M).
Relaxation of ell_j = log10 D_j by a secant/tangent sandwich on certified brackets.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import cvxpy as cp
import numpy as np
from scipy.linalg import cholesky, solve_triangular

from . import C_KM_S
from .data import BAO, SN
from .geometry import dh_matrix, dm_matrix, grid, segment_increment_factor

LN10 = np.log(10.0)


@dataclass(frozen=True)
class ClassSpec:
    N: int = 50               # segments (uniform grid) — ignored when grid_kind == "geometric"
    z_max: float = 2.5
    L: float = 1.0            # slope bound; np.inf allowed (uniform grid / v1 model only)
    H_min: float = 20.0       # a priori box on H(z) [km/s/Mpc], part of the class definition
    H_max: float = 2000.0
    r_lo: float = 146.57      # Mpc, r_d box
    r_hi: float = 147.61
    grid_kind: str = "uniform"   # "uniform" | "geometric"
    z_first: float = 0.005
    ratio: float = 1.1
    h_max: float = 0.02

    @property
    def x(self) -> np.ndarray:
        if self.grid_kind == "uniform":
            return grid(self.N, self.z_max)
        from .geometry import geometric_grid
        return geometric_grid(self.z_max, self.z_first, self.ratio, self.h_max)

    @property
    def n_seg(self) -> int:
        return len(self.x) - 1

    @property
    def hs(self) -> np.ndarray:
        return np.diff(self.x)

    @property
    def h(self) -> float:
        hs = self.hs
        if not np.allclose(hs, hs[0]):
            raise ValueError("non-uniform grid: use hs")
        return float(hs[0])

    @property
    def u_box(self) -> tuple[float, float]:
        return C_KM_S / (self.H_max * self.r_hi), C_KM_S / (self.H_min * self.r_lo)


def whitener(C: np.ndarray) -> np.ndarray:
    """W with W @ C @ W.T = I (W = L^{-1}, C = L L^T)."""
    Lc = cholesky(C, lower=True)
    return solve_triangular(Lc, np.eye(C.shape[0]), lower=True)


@dataclass
class Frozen:
    bao: BAO
    sn: SN
    spec: ClassSpec
    P: np.ndarray = field(init=False)      # BAO prediction matrix (rows in data order)
    A_sn: np.ndarray = field(init=False)   # D_M/r_d at SN zHD
    A_nodes: np.ndarray = field(init=False)  # D_M/r_d at grid nodes
    Wb: np.ndarray = field(init=False)
    Wsn: np.ndarray = field(init=False)
    sn_offset: np.ndarray = field(init=False)  # 5 log10(1+zHEL)

    def __post_init__(self):
        x = self.spec.x
        AM = dm_matrix(x, self.bao.z)
        AH = dh_matrix(x, self.bao.z)
        rows = []
        for r, kind in enumerate(self.bao.kind):
            if kind == "DM_over_rs":
                rows.append(AM[r])
            elif kind == "DH_over_rs":
                rows.append(AH[r])
            else:
                raise ValueError(f"unsupported BAO kind {kind} (drop D_V rows)")
        self.P = np.array(rows)
        self.A_sn = dm_matrix(x, self.sn.zHD)
        self.A_nodes = dm_matrix(x, np.expm1(x))
        self.Wb = whitener(self.bao.cov)
        self.Wsn = whitener(self.sn.cov)
        self.sn_offset = 5.0 * np.log10(1.0 + self.sn.zHEL)

    # ---- exact statistic at a point (true logarithm, no relaxation) ----
    def chi2(self, u: np.ndarray, Mp: float) -> float:
        rb = self.bao.value - self.P @ u
        D = self.A_sn @ u
        rs = self.sn.m - 5.0 * np.log10(D) - self.sn_offset - Mp
        return float(np.sum((self.Wb @ rb) ** 2) + np.sum((self.Wsn @ rs) ** 2))

    def best_Mp(self, u: np.ndarray) -> float:
        """Analytic minimizer of the SN term over Mp for fixed u."""
        D = self.A_sn @ u
        y = self.sn.m - 5.0 * np.log10(D) - self.sn_offset
        Cinv1 = self.Wsn.T @ (self.Wsn @ np.ones_like(y))
        return float((Cinv1 @ y) / (Cinv1 @ np.ones_like(y)))


@dataclass
class Brackets:
    """Certified per-SN brackets D_j in [lo, hi] and tangent points t."""
    lo: np.ndarray
    hi: np.ndarray
    t: np.ndarray

    @staticmethod
    def from_nodes(fr: Frozen, u_lo, u_hi, dm_lo, dm_hi, t=None) -> "Brackets":
        k, fac = segment_increment_factor(fr.spec.x, fr.sn.zHD)
        lo = dm_lo[k] + np.minimum(u_lo[k], u_lo[k + 1]) * fac
        hi = dm_hi[k] + np.maximum(u_hi[k], u_hi[k + 1]) * fac
        if t is None:
            t = np.sqrt(lo * hi)
        return Brackets(lo=lo, hi=hi, t=np.clip(t, lo, hi))


class Model:
    """One cvxpy problem with parameterized linear objective (DPP), reused for all solves."""

    def __init__(self, fr: Frozen, br: Brackets, T: float | None):
        self.fr, self.br, self.T = fr, br, T
        sp = fr.spec
        N = sp.n_seg
        self.u = cp.Variable(N + 1, name="u")
        self.Mp = cp.Variable(name="Mp")
        self.ell = cp.Variable(len(fr.sn.m), name="ell")
        self.cobj = cp.Parameter(N + 1, name="cobj")
        cons = []
        ulo, uhi = sp.u_box
        cons += [self.u >= ulo, self.u <= uhi]
        if np.isfinite(sp.L):
            Lh = sp.L * sp.hs
            d = self.u[1:] - self.u[:-1]
            cons += [d <= cp.multiply(Lh, self.u[:-1]), d <= cp.multiply(Lh, self.u[1:]),
                     -d <= cp.multiply(Lh, self.u[:-1]), -d <= cp.multiply(Lh, self.u[1:])]
        D = fr.A_sn @ self.u
        lo, hi, t = br.lo, br.hi, br.t
        slope = (np.log10(hi) - np.log10(lo)) / (hi - lo)
        cons += [D >= lo, D <= hi,
                 self.ell >= np.log10(lo) + cp.multiply(slope, D - lo),        # secant (needs bracket)
                 self.ell <= np.log10(t) + cp.multiply(1.0 / (t * LN10), D - t)]  # tangent (global)
        rb = fr.bao.value - fr.P @ self.u
        rs = fr.sn.m - 5.0 * self.ell - fr.sn_offset - self.Mp
        self.chi2 = cp.sum_squares(fr.Wb @ rb) + cp.sum_squares(fr.Wsn @ rs)
        if T is not None:
            cons += [self.chi2 <= T]
        self.cons = cons
        self.prob_lin = cp.Problem(cp.Minimize(self.cobj @ self.u), cons)
        self.prob_chi2 = cp.Problem(cp.Minimize(self.chi2), cons)

    def _solve(self, prob, **kw):
        prob.solve(solver=cp.CLARABEL, **kw)
        if prob.status not in ("optimal", "optimal_inaccurate"):
            raise RuntimeError(f"solver status {prob.status}")
        return prob.value

    def min_chi2(self):
        v = self._solve(self.prob_chi2)
        return v, self.u.value.copy(), float(self.Mp.value)

    def extremize(self, c: np.ndarray) -> float:
        """min c @ u over the current (relaxed) feasible set."""
        self.cobj.value = c
        return self._solve(self.prob_lin)

    def node_brackets(self):
        """Certified min/max of u_i and D_M(x_i)/r_d over the current relaxed set."""
        N = self.fr.spec.n_seg
        u_lo = np.empty(N + 1); u_hi = np.empty(N + 1)
        dm_lo = np.empty(N + 1); dm_hi = np.empty(N + 1)
        A = self.fr.A_nodes
        for i in range(N + 1):
            e = np.zeros(N + 1); e[i] = 1.0
            u_lo[i] = self.extremize(e)
            u_hi[i] = -self.extremize(-e)
            if i == 0:
                dm_lo[i] = dm_hi[i] = 0.0
            else:
                dm_lo[i] = self.extremize(A[i])
                dm_hi[i] = -self.extremize(-A[i])
        return u_lo, u_hi, dm_lo, dm_hi


def initial_brackets(fr: Frozen, T: float) -> Brackets:
    """Brackets from the outer set {class box, slope, BAO chi2 <= T} (SN dropped)."""
    u_lo, u_hi, dm_lo, dm_hi = bao_only_bounds(fr, T)
    return Brackets.from_nodes(fr, u_lo, u_hi, dm_lo, dm_hi)


def bao_only_bounds(fr: Frozen, T: float):
    """Node bounds over the outer set {class box, slope, BAO chi2 <= T} (SN dropped: valid, loose)."""
    sp = fr.spec
    N = sp.n_seg
    u = cp.Variable(N + 1)
    c = cp.Parameter(N + 1)
    ulo, uhi = sp.u_box
    cons = [u >= ulo, u <= uhi]
    if np.isfinite(sp.L):
        Lh = sp.L * sp.hs
        d = u[1:] - u[:-1]
        cons += [d <= cp.multiply(Lh, u[:-1]), d <= cp.multiply(Lh, u[1:]),
                 -d <= cp.multiply(Lh, u[:-1]), -d <= cp.multiply(Lh, u[1:])]
    cons += [cp.sum_squares(fr.Wb @ (fr.bao.value - fr.P @ u)) <= T]
    prob = cp.Problem(cp.Minimize(c @ u), cons)

    def ext(v):
        c.value = v
        prob.solve(solver=cp.CLARABEL)
        if prob.status not in ("optimal", "optimal_inaccurate"):
            raise RuntimeError(prob.status)
        return prob.value

    u_lo = np.empty(N + 1); u_hi = np.empty(N + 1); dm_lo = np.zeros(N + 1); dm_hi = np.zeros(N + 1)
    for i in range(N + 1):
        e = np.zeros(N + 1); e[i] = 1.0
        u_lo[i] = ext(e); u_hi[i] = -ext(-e)
        if i > 0:
            dm_lo[i] = ext(fr.A_nodes[i]); dm_hi[i] = -ext(-fr.A_nodes[i])
    return u_lo, u_hi, dm_lo, dm_hi


@dataclass
class BoundResult:
    L: float
    Delta: float
    N: int
    T: float
    chi2_ref: float
    u0_min: float
    H0_max: float          # c/(r_lo u0_min)
    H0_ref: float          # c/(147.09 u_ref[0]) for context
    u_ref: np.ndarray
    Mp_ref: float
    iterations: int
    history: list


def certified_bound(fr: Frozen, Delta: float, u_ref: np.ndarray, Mp_ref: float | None = None,
                    n_iter: int = 3, tol: float = 5e-3, verbose: bool = True) -> BoundResult:
    """Compute H0_max for the class in fr.spec (FORMULATION §1) by iterated bracket tightening.

    u_ref must lie in the class; it seeds T = chi2(u_ref) + Delta (any feasible point is valid).
    """
    sp = fr.spec
    if Mp_ref is None:
        Mp_ref = fr.best_Mp(u_ref)
    chi2_ref = fr.chi2(u_ref, Mp_ref)
    T = chi2_ref + Delta
    br = initial_brackets(fr, T)
    history = []
    u0_min = None
    for it in range(n_iter):
        m = Model(fr, br, T)
        # (a) improve the reference point: relaxed chi2 minimizer, re-evaluated exactly
        _, u_c, _ = m.min_chi2()
        Mp_c = fr.best_Mp(u_c)
        chi2_c = fr.chi2(u_c, Mp_c)
        if chi2_c < chi2_ref:
            u_ref, Mp_ref, chi2_ref = u_c, Mp_c, chi2_c
            T = chi2_ref + Delta
            m = Model(fr, br, T)
        # (b) bound
        e0 = np.zeros(sp.N + 1); e0[0] = 1.0
        new_u0 = m.extremize(e0)
        H = C_KM_S / (sp.r_lo * new_u0)
        history.append(dict(iter=it, T=T, chi2_ref=chi2_ref, u0_min=new_u0, H0_max=H,
                            bracket_width_max=float(np.max(np.log10(br.hi / br.lo)))))
        if verbose:
            print(f"  it={it} chi2_ref={chi2_ref:.3f} T={T:.3f} u0_min={new_u0:.5f} "
                  f"H0_max={H:.3f} maxlog10(hi/lo)={history[-1]['bracket_width_max']:.4f}")
        if u0_min is not None and abs(new_u0 - u0_min) / u0_min < tol / 10:
            u0_min = new_u0
            break
        u0_min = new_u0
        # (c) tighten brackets over the current relaxed set (still a superset of F)
        u_lo, u_hi, dm_lo, dm_hi = m.node_brackets()
        # tangent points at the current chi2-minimizer's distances
        br = Brackets.from_nodes(fr, u_lo, u_hi, dm_lo, dm_hi, t=fr.A_sn @ u_ref)
    return BoundResult(L=sp.L, Delta=Delta, N=sp.N, T=T, chi2_ref=chi2_ref, u0_min=u0_min,
                       H0_max=C_KM_S / (sp.r_lo * u0_min), H0_ref=C_KM_S / (147.09 * u_ref[0]),
                       u_ref=u_ref, Mp_ref=Mp_ref, iterations=len(history), history=history)
