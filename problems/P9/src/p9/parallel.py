"""Parallel bound tightening for the LKR model: each OBBT objective is an independent conic solve."""

from __future__ import annotations

import multiprocessing as mp
import os
import time

import numpy as np

from .lkr import Brackets3, LKRModel
from .model import Frozen

_MODEL = None


def _init(fr: Frozen, br: Brackets3, T: float):
    global _MODEL
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    _MODEL = LKRModel(fr, br, T)


def _solve_obj(q):
    return _MODEL.extremize(q)


def tighten_parallel(fr: Frozen, br: Brackets3, T: float, n_workers: int = 6, lam_nodes=None,
                     verbose: bool = True) -> Brackets3:
    m = LKRModel(fr, br, T)  # only for indices / c_node
    N = m.N
    if lam_nodes is None:
        lam_nodes = sorted(set(m.enodes) | {0})
    jobs = []   # (kind, node, sign, q)
    for i in range(1, N + 1):
        q = np.zeros(m.nvar); q[m.idx["lam"][i]] = 1.0; q[m.idx["kappa"][i]] = -1.0
        jobs.append(("rho", i, +1, q)); jobs.append(("rho", i, -1, -q))
    for i in lam_nodes:
        q = np.zeros(m.nvar); q[m.idx["lam"][i]] = 1.0
        jobs.append(("lam", i, +1, q)); jobs.append(("lam", i, -1, -q))
    for p in range(len(m.idx_dm)):
        q = np.zeros(m.nvar); q[m.idx["yb"][p]] = 1.0
        jobs.append(("yb", p, +1, q)); jobs.append(("yb", p, -1, -q))
    t0 = time.time()
    # "spawn": fork after BLAS/Rust threads have started deadlocks on macOS
    ctx = mp.get_context("spawn")
    with ctx.Pool(n_workers, initializer=_init, initargs=(fr, br, T)) as pool:
        vals = pool.map(_solve_obj, [j[3] for j in jobs], chunksize=4)
    rho_lo = br.rho_lo.copy(); rho_hi = br.rho_hi.copy()
    lam_lo = br.lam_lo.copy(); lam_hi = br.lam_hi.copy()
    yb_lo = br.yb_lo.copy(); yb_hi = br.yb_hi.copy()
    for (kind, i, sgn, _), v in zip(jobs, vals):
        if kind == "rho":
            if sgn > 0: rho_lo[i] = max(rho_lo[i], v - m.c_node[i])
            else: rho_hi[i] = min(rho_hi[i], -v - m.c_node[i])
        elif kind == "lam":
            if sgn > 0: lam_lo[i] = max(lam_lo[i], v)
            else: lam_hi[i] = min(lam_hi[i], -v)
        else:
            if sgn > 0: yb_lo[i] = max(yb_lo[i], v)
            else: yb_hi[i] = min(yb_hi[i], -v)
    if verbose:
        print(f"    tighten_parallel: {len(jobs)} solves on {n_workers} workers in {time.time()-t0:.0f}s", flush=True)
    return Brackets3(rho_lo, rho_hi, yb_lo, yb_hi, lam_lo, lam_hi)
