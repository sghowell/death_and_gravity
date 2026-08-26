"""Bound computation with the kappa relaxation and subset-accelerated bound tightening.

Validity of the subset trick: for any subset S of the SNe, the marginal chi2 over S is <= the
full chi2 (Gaussian marginalization), so {class, chi2_S <= T} ⊇ F. Node bounds computed over
the subset problem are therefore valid for F. Reference points and the final bound use the
full sample.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import C_KM_S
from .data import SN
from .classmin import minimize_chi2_over_class
from .model import Frozen, bao_only_bounds
from .socp2 import KappaModel, NodeBounds


def sn_subset(sn: SN, size: int) -> SN:
    """Stratified-in-redshift subset (every k-th SN in redshift order)."""
    order = np.argsort(sn.zHD)
    k = max(1, len(order) // size)
    pick = np.sort(order[::k][:size])
    return SN(zHD=sn.zHD[pick], zHEL=sn.zHEL[pick], m=sn.m[pick],
              cov=sn.cov[np.ix_(pick, pick)], index=sn.index[pick])


@dataclass
class BoundResult:
    L: float
    Delta: float
    T: float
    chi2_ref: float
    u_ref: np.ndarray
    Mp_ref: float
    u0_min: float
    H0_max: float
    H0_ref: float
    node_bounds: NodeBounds
    history: list


def compute_bound(fr_full: Frozen, fr_sub: Frozen, Delta: float, u_ref: np.ndarray,
                  nb: NodeBounds | None = None, n_passes: int = 5, tol_rel: float = 1e-4,
                  verbose: bool = True) -> BoundResult:
    sp = fr_full.spec
    # reference point: true chi2 minimum over the class (local, polished, exactly feasible)
    u_ref, Mp_ref, chi2_ref = minimize_chi2_over_class(fr_full, u_ref, verbose=verbose)
    T = chi2_ref + Delta
    if nb is None:
        nb = NodeBounds(*bao_only_bounds(fr_full, T))
    history = []
    N = sp.n_seg
    e0 = np.zeros(N + 1); e0[0] = 1.0
    last = None
    u0_min = None
    for it in range(n_passes):
        # (a) try to improve the reference point from the relaxed chi2 minimizer (full sample)
        m_full = KappaModel(fr_full, nb, T, tangent_u=u_ref)
        chi2_rel, u_c, _ = m_full.min_chi2()
        u_c, Mp_c, chi2_c = minimize_chi2_over_class(fr_full, u_c)
        if chi2_c < chi2_ref - 1e-9:
            u_ref, Mp_ref, chi2_ref = u_c, Mp_c, chi2_c
            T = chi2_ref + Delta
            m_full = KappaModel(fr_full, nb, T, tangent_u=u_ref)
        # (b) bound with the full sample
        u0_min = m_full.extremize(e0)
        H0 = C_KM_S / (sp.r_lo * u0_min)
        history.append(dict(pass_=it, T=T, chi2_ref=chi2_ref, chi2_relaxed_min=chi2_rel,
                            u0_min=u0_min, H0_max=H0, width=nb.width()))
        if verbose:
            print(f"  pass {it}: chi2_ref={chi2_ref:.3f} (relaxed min {chi2_rel:.3f}) T={T:.3f} "
                  f"u0_min={u0_min:.5f} H0_max={H0:.3f} node-width(log10)={nb.width():.4f}", flush=True)
        if last is not None and abs(u0_min - last) <= tol_rel * last:
            break
        last = u0_min
        # (c) tighten node bounds on the subset model (valid outer set), intersected with previous
        m_sub = KappaModel(fr_sub, nb, T, tangent_u=u_ref)
        nb = m_sub.node_bounds()
    return BoundResult(L=sp.L, Delta=Delta, T=T, chi2_ref=chi2_ref, u_ref=u_ref, Mp_ref=Mp_ref,
                       u0_min=u0_min, H0_max=C_KM_S / (sp.r_lo * u0_min),
                       H0_ref=C_KM_S / (147.09 * u_ref[0]), node_bounds=nb, history=history)
