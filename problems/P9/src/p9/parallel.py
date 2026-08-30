"""Parallel bound tightening for the LKR model: each OBBT objective is an independent conic solve."""

from __future__ import annotations

import multiprocessing as mp
import os
import time

import numpy as np

from .lkr import Brackets3, LKRModel
from .lkr_rows import obbt_objectives
from .model import Frozen

_MODEL = None


def _model_cls(fr: Frozen):
    """The independent LKRModel for the baseline; the shared-row LKRModel2 when the D_V row is on."""
    if fr.spec.use_dv:
        from .lkr2 import LKRModel2
        return LKRModel2
    return LKRModel


def _init(fr: Frozen, br: Brackets3, T: float):
    global _MODEL
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    _MODEL = _model_cls(fr)(fr, br, T)


def _solve_obj(q):
    return _MODEL.extremize(q)


def tighten_parallel(fr: Frozen, br: Brackets3, T: float, n_workers: int = 6, lam_nodes=None,
                     verbose: bool = True) -> Brackets3:
    from .lkr_rows import layout_for
    lay = layout_for(fr); c_node = np.concatenate([[0.0], np.log10(np.expm1(fr.spec.x[1:]))])
    jobs = []   # (kind, i, side, q)
    for kind, i, side, qd in obbt_objectives(lay, lam_nodes):
        q = np.zeros(lay.nvar)
        for v, cf in qd.items():
            q[v] = cf
        jobs.append((kind, i, side, q))
    t0 = time.time()
    # "spawn": fork after BLAS/Rust threads have started deadlocks on macOS
    ctx = mp.get_context("spawn")
    with ctx.Pool(n_workers, initializer=_init, initargs=(fr, br, T)) as pool:
        vals = pool.map(_solve_obj, [j[3] for j in jobs], chunksize=4)
    out = br.copy()
    for (kind, i, side, _), v in zip(jobs, vals):
        out.apply_bound(kind, i, side, v, float(c_node[i]) if kind == "rho" else 0.0)
    if verbose:
        print(f"    tighten_parallel: {len(jobs)} solves on {n_workers} workers in {time.time()-t0:.0f}s", flush=True)
    return out
