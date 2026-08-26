"""Certified LKR chain: every bound-tightening solve and the final bound are verified in ball
arithmetic (Verifier3); brackets passed between passes are the rigorous outward-rounded values;
T is a rigorous enclosure of chi2 at the (exactly class-feasible) reference point plus Delta.
Dual vectors are stored under results/certificates/lkr_<tag>/.

usage: PYTHONPATH=problems/P9/src uv run python -m p9.run_lkr_certified --L 1.5 --Delta 4 --refine 2 --workers 8 --passes 8
Resumes from results/certificates/lkr_<tag>/state.json if present.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

from . import C_KM_S
from .certify_feasible import in_class_exact
from .classmin import minimize_chi2_over_class
from .data import load_desi, load_pantheon, verify_manifest
from .geometry import lcdm_u_nodes
from .lcdm import fit_bao_sn
from .lkr import Brackets3, initial_brackets3
from .lkr2 import LKRModel2
from .model import ClassSpec, Frozen
from .verify import _endpoint, rigorous_chi2
from .verify3 import Verifier3

RESULTS = Path(__file__).resolve().parents[2] / "results"

_M = None; _V = None; _DIR = None


def _init(fr, br, T, cert_dir):
    global _M, _V, _DIR
    _M = LKRModel2(fr, br, T); _V = Verifier3(fr, br, T); _DIR = Path(cert_dir)


def _job(args):
    """args = (tag, qdict). Solve, certify, store the dual, return the rigorous bound."""
    tag, qd = args
    q = np.zeros(_M.nvar)
    for v, cf in qd.items():
        q[v] = cf
    try:
        val, x, z = _M.solve_dual(q)
    except RuntimeError as e:      # solver failure: no tightening for this quantity (still valid)
        return tag, None, None
    n_eq, n_in = _M._n_eq, _M._n_in
    lb = _V.certify(z[:n_eq], z[n_eq:n_eq + n_in], z[n_eq + n_in:], qd, verbose=False)
    np.savez_compressed(_DIR / f"{tag}.npz", z=z, q=list(qd.items()), primal=val, rigorous=lb)
    return tag, val, lb


def certified_pass(fr, br, T, cert_dir, pass_id, workers, lam_nodes=None):
    m = LKRModel2(fr, br, T)          # for layout only
    lay = m.lay; N = lay.N
    if lam_nodes is None:
        lam_nodes = sorted(set(lay.enodes) | {0})
    jobs = []
    for i in range(1, N + 1):
        jobs.append((f"p{pass_id}_rho{i}_lo", {int(lay.lam[i]): 1.0, int(lay.kappa[i]): -1.0}))
        jobs.append((f"p{pass_id}_rho{i}_hi", {int(lay.lam[i]): -1.0, int(lay.kappa[i]): 1.0}))
    for i in lam_nodes:
        jobs.append((f"p{pass_id}_lam{i}_lo", {int(lay.lam[i]): 1.0}))
        jobs.append((f"p{pass_id}_lam{i}_hi", {int(lay.lam[i]): -1.0}))
    for p in range(len(lay.idx_dm)):
        jobs.append((f"p{pass_id}_yb{p}_lo", {int(lay.yb[p]): 1.0}))
        jobs.append((f"p{pass_id}_yb{p}_hi", {int(lay.yb[p]): -1.0}))
    t0 = time.time()
    ctx = mp.get_context("spawn")
    with ctx.Pool(workers, initializer=_init, initargs=(fr, br, T, str(cert_dir))) as pool:
        res = pool.map(_job, jobs, chunksize=2)
    c_node = m.c_node
    rho_lo = br.rho_lo.copy(); rho_hi = br.rho_hi.copy(); lam_lo = br.lam_lo.copy(); lam_hi = br.lam_hi.copy()
    yb_lo = br.yb_lo.copy(); yb_hi = br.yb_hi.copy()
    worst_gap = 0.0
    n_fail = sum(1 for _, val, _ in res if val is None)
    for tag, val, lb in res:
        if val is None:
            continue
        worst_gap = max(worst_gap, abs(val - lb))
        kind = tag.split("_")[1]
        i = int(kind[3:]) if kind.startswith("rho") else int(kind[3:]) if kind.startswith("lam") else int(kind[2:])
        side = tag.split("_")[2]
        if kind.startswith("rho"):
            if side == "lo": rho_lo[i] = max(rho_lo[i], lb - c_node[i])
            else: rho_hi[i] = min(rho_hi[i], -lb - c_node[i])
        elif kind.startswith("lam"):
            if side == "lo": lam_lo[i] = max(lam_lo[i], lb)
            else: lam_hi[i] = min(lam_hi[i], -lb)
        else:
            if side == "lo": yb_lo[i] = max(yb_lo[i], lb)
            else: yb_hi[i] = min(yb_hi[i], -lb)
    print(f"    certified pass {pass_id}: {len(jobs)} solves+certificates on {workers} workers in {time.time()-t0:.0f}s; "
          f"max |solver - rigorous| = {worst_gap:.2e}; solver failures (left untightened): {n_fail}", flush=True)
    return Brackets3(rho_lo, rho_hi, yb_lo, yb_hi, lam_lo, lam_hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=float, default=1.5); ap.add_argument("--Delta", type=float, default=4.0)
    ap.add_argument("--refine", type=int, default=2); ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--passes", type=int, default=8); ap.add_argument("--tol", type=float, default=2e-5)
    args = ap.parse_args()
    verify_manifest()
    bao = load_desi(); sn = load_pantheon()
    spec = ClassSpec(L=args.L, grid_kind="geometric", refine=args.refine); fr = Frozen(bao, sn, spec)
    tag = f"L{args.L:g}_D{args.Delta:g}_r{args.refine}"
    cert_dir = RESULTS / "certificates" / f"lkr_{tag}"; cert_dir.mkdir(parents=True, exist_ok=True)
    state = cert_dir / "state.json"
    if state.exists():
        s = json.loads(state.read_text())
        T = s["T"]; br = Brackets3(*(np.array(s[k]) for k in ["rho_lo", "rho_hi", "yb_lo", "yb_hi", "lam_lo", "lam_hi"]))
        hist = s["history"]; ref = s["reference"]
        print(f"resuming {tag}: T={T:.6f} passes done={len(hist)}", flush=True)
    else:
        (om, hrd), _ = fit_bao_sn(bao, sn)
        u = lcdm_u_nodes(spec.x, om, hrd)
        cands = [minimize_chi2_over_class(fr, u), minimize_chi2_over_class(fr, np.full(spec.n_seg + 1, u[0]))]
        u_star, Mp_star, chi2_star = min(cands, key=lambda c: c[2])
        assert in_class_exact(spec, u_star)
        ball = rigorous_chi2(fr, u_star, Mp_star)
        T = _endpoint(ball, +1) + args.Delta
        ref = dict(u=u_star.tolist(), Mp=Mp_star, chi2_float=chi2_star, chi2_enclosure=ball.str(20), T=T)
        br = initial_brackets3(fr); hist = []
        print(f"{tag}: class-min chi2 {ball.str(12)} -> T={T:.6f}; nodes={spec.n_seg+1}", flush=True)
    last = hist[-1]["lambda0_min"] if hist else None
    for it in range(len(hist), len(hist) + args.passes):
        t = time.time()
        br = certified_pass(fr, br, T, cert_dir, it, args.workers)
        # certified bound on lambda_0 with the new brackets (single solve in the main process)
        m = LKRModel2(fr, br, T); ver = Verifier3(fr, br, T)
        qd = {int(m.lay.lam[0]): 1.0}
        val, x, z = m.solve_dual(np.eye(m.nvar)[m.lay.lam[0]])
        lb = ver.certify(z[:m._n_eq], z[m._n_eq:m._n_eq + m._n_in], z[m._n_eq + m._n_in:], qd, verbose=False)
        np.savez_compressed(cert_dir / f"p{it}_final_lambda0.npz", z=z, primal=val, rigorous=lb, T=T)
        H0 = C_KM_S / (spec.r_lo * 10 ** lb)
        hist.append(dict(pass_=it, lambda0_min=lb, solver=val, H0_max=H0, widths=br.width(), seconds=time.time() - t))
        print(f"pass {it}: rigorous lambda0_min={lb:.6f} (solver {val:.6f}) -> H0_max={H0:.4f}  widths {br.width()} [{time.time()-t:.0f}s]", flush=True)
        state.write_text(json.dumps(dict(tag=tag, L=args.L, Delta=args.Delta, refine=args.refine, T=T, reference=ref,
                                         history=hist, H0_max=H0, lambda0_min=lb, r_lo=spec.r_lo,
                                         **{k: getattr(br, k).tolist() for k in ["rho_lo", "rho_hi", "yb_lo", "yb_hi", "lam_lo", "lam_hi"]}), indent=1))
        if last is not None and abs(lb - last) < args.tol:
            break
        last = lb
    print("TOTAL done", flush=True)


if __name__ == "__main__":
    main()
