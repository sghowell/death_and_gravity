"""Dissect a relaxed optimum: where does the relaxation cheat?

usage: PYTHONPATH=problems/P9/src uv run python -m p9.diagnose results/run_L1.5_D4_v3.json
"""

from __future__ import annotations

import json
import sys

import clarabel
import numpy as np

from . import C_KM_S
from .bound import sn_subset
from .data import load_desi, load_pantheon
from .model import ClassSpec, Frozen
from .socp2 import KappaModel, NodeBounds


def main(path: str, L: float = 1.5):
    r = json.load(open(path))
    bao = load_desi(); sn = load_pantheon()
    spec = ClassSpec(L=L, grid_kind="geometric"); x = spec.x; z = np.expm1(x)
    fr = Frozen(bao, sn, spec)
    nb = NodeBounds(np.array(r["u_lo"]), np.array(r["u_hi"]), np.array(r["dm_lo"]), np.array(r["dm_hi"]))
    u_ref = np.array(r["u_ref"]); T = r["T"]
    N = spec.n_seg
    print("node widths (log10 hi/lo) for u and D:")
    for i in range(0, N + 1, 6):
        wu = np.log10(nb.u_hi[i] / nb.u_lo[i]); wd = np.log10(nb.dm_hi[i] / nb.dm_lo[i]) if i > 0 else 0
        print(f"  z={z[i]:.4f}  u: {wu:.4f} ({(10**wu-1)*100:5.1f}%)   D: {wd:.4f} ({(10**wd-1)*100:5.1f}%)")
    m = KappaModel(fr, nb, T, tangent_u=u_ref)
    q = np.zeros(m.nvar); q[m.idx["u"][0]] = 1.0
    solver = clarabel.DefaultSolver(m.P0, q, m.A, m.b, m.cones, m.settings)
    sol = solver.solve(); y = np.asarray(sol.x)
    u = y[m.idx["u"]]; kap = y[m.idx["kappa"]]; Mp = y[m.idx["Mp"]]; ell = y[m.idx["ell"]]
    print(f"\nrelaxed optimum: u0={u[0]:.5f} -> H0={C_KM_S/(spec.r_lo*u[0]):.3f};  true chi2 at (u, Mp*)={fr.chi2(u, fr.best_Mp(u)):.2f} (T={T:.2f})")
    c = np.concatenate([[0.0], np.log10(np.expm1(x[1:]))])
    D = np.concatenate([[u[0]], (fr.A_nodes @ u)[1:]])
    gap_used = kap + c - np.log10(D)          # kappa-side minus u-side (log10)
    print("\nlink gap usage (kappa + c - log10 D) in mag (x5), by node:")
    for i in range(0, N + 1, 6):
        lo = nb.u_lo[0] if i == 0 else nb.dm_lo[i]; hi = nb.u_hi[0] if i == 0 else nb.dm_hi[i]
        maxgap = (np.log(hi / lo)) ** 2 / (8 * np.log(10))
        print(f"  z={z[i]:.4f}: used {5*gap_used[i]:+.4f} mag  (max possible {5*maxgap:.4f})  u/u_ref={u[i]/u_ref[i]:.4f}")
    Dj = fr.A_sn @ u
    dev = ell - np.log10(Dj)                  # per-SN deviation of modelled log-distance from the u-side truth
    print(f"\nper-SN |ell - log10 D_j(u)| in mag: mean {5*np.mean(np.abs(dev)):.4f}, max {5*np.max(np.abs(dev)):.4f}")
    for lo_, hi_ in [(0.01, 0.05), (0.05, 0.2), (0.2, 0.5), (0.5, 1.0), (1.0, 2.3)]:
        s = (sn.zHD > lo_) & (sn.zHD <= hi_)
        print(f"  z in ({lo_},{hi_}]: n={s.sum():4d}  mean dev {5*np.mean(dev[s]):+.4f} mag  rms {5*np.std(dev[s]):.4f}")
    # chi2 decomposition at the relaxed optimum vs true
    rs_rel = sn.m - 5 * ell - fr.sn_offset - Mp
    rs_true = sn.m - 5 * np.log10(Dj) - fr.sn_offset - Mp
    print(f"\nSN chi2: relaxed {np.sum((fr.Wsn @ rs_rel)**2):.2f}  true(same u, Mp) {np.sum((fr.Wsn @ rs_true)**2):.2f}")
    rb = bao.value - fr.P @ u
    print(f"BAO chi2 at relaxed optimum: {np.sum((fr.Wb @ rb)**2):.2f};  BAO residuals/sigma: {np.round(rb/np.sqrt(np.diag(bao.cov)), 2)}")


if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 1.5)
