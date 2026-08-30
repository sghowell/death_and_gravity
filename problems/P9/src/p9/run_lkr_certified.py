"""Certified LKR chain: every bound-tightening solve and the final bound are verified in ball
arithmetic (Verifier3); brackets passed between passes are the rigorous outward-rounded values;
T is a rigorous enclosure of chi2 at the (exactly class-feasible) reference point plus Delta.
Dual vectors are stored under results/certificates/lkr_<tag>/.

usage: PYTHONPATH=problems/P9/src uv run python -m p9.run_lkr_certified --L 1.5 --Delta 4 --refine 2 --workers 8 --passes 8
       [--sn {pantheon,dessn5yr,union3}] [--dv] [--rd_box LO HI | planck | bbn]      (FORMULATION §6.1 variants)
Resumes from results/certificates/lkr_<tag>/state.json if present. The tag is L<L>_D<Delta>_r<refine> plus the
variant suffix (empty for the baseline), so baseline certificate directories are unchanged.
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
from .geometry import lcdm_u_nodes
from .lcdm import fit_bao_sn
from .lkr import Brackets3, initial_brackets3
from .lkr2 import LKRModel2
from .lkr_rows import obbt_objectives
from .variants import Variant, add_variant_args
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


def c_node_balls(fr):
    """Rigorous c_i = log10(e^{x_i} - 1) (c_0 = 0) as Arb balls, for the outward-rounded rho update."""
    from flint import arb
    from .lkr_rows import ArbArith
    ar = ArbArith()
    return [arb(0)] + [ar.log10(ar.expm1(ar.c(v))) for v in fr.spec.x[1:]]


def certified_pass(fr, br, T, cert_dir, pass_id, workers, lam_nodes=None):
    lay = LKRModel2(fr, br, T).lay          # for layout only
    objs = obbt_objectives(lay, lam_nodes)
    jobs = [(f"p{pass_id}_{kind}{i}_{side}", qd) for kind, i, side, qd in objs]
    t0 = time.time()
    ctx = mp.get_context("spawn")
    with ctx.Pool(workers, initializer=_init, initargs=(fr, br, T, str(cert_dir))) as pool:
        res = pool.map(_job, jobs, chunksize=2)
    c_ball = c_node_balls(fr)
    out = br.copy()
    worst_gap = 0.0
    n_fail = sum(1 for _, val, _ in res if val is None)
    for (kind, i, side, _), (tag, val, lb) in zip(objs, res):
        if val is None:
            continue
        worst_gap = max(worst_gap, abs(val - lb))
        out.apply_bound(kind, i, side, lb, c_ball[i] if kind == "rho" else 0.0)
    print(f"    certified pass {pass_id}: {len(jobs)} solves+certificates on {workers} workers in {time.time()-t0:.0f}s; "
          f"max |solver - rigorous| = {worst_gap:.2e}; solver failures (left untightened): {n_fail}", flush=True)
    return out


def reference_point(fr, bao, sn, spec, Delta):
    """Class minimizer (from the LCDM fit and from a flat start), exactly in class, with rigorous T."""
    (om, hrd), _ = fit_bao_sn(bao, sn)
    u = lcdm_u_nodes(spec.x, om, hrd)
    cands = [minimize_chi2_over_class(fr, u), minimize_chi2_over_class(fr, np.full(spec.n_seg + 1, u[0]))]
    u_star, Mp_star, chi2_star = min(cands, key=lambda c: c[2])
    from flint import arb
    from .socp2 import MP_BOX
    assert in_class_exact(spec, u_star) and MP_BOX[0] <= Mp_star <= MP_BOX[1]
    ball = rigorous_chi2(fr, u_star, Mp_star)
    T = _endpoint(ball + arb(Delta), +1)
    ref = dict(u=u_star.tolist(), Mp=Mp_star, chi2_float=chi2_star, chi2_enclosure=ball.str(20), T=T,
               lcdm_fit=dict(omega_m=float(om), h_rd=float(hrd)))
    return ref, T, ball


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=float, default=1.5); ap.add_argument("--Delta", type=float, default=4.0)
    ap.add_argument("--refine", type=int, default=2); ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--passes", type=int, default=8); ap.add_argument("--tol", type=float, default=2e-5)
    add_variant_args(ap)
    args = ap.parse_args()
    var = Variant.from_args(args)
    bao, sn, spec, fr = var.frozen(args.L, args.refine)
    tag = f"L{args.L:g}_D{args.Delta:g}_r{args.refine}" + var.tag
    cert_dir = RESULTS / "certificates" / f"lkr_{tag}"; cert_dir.mkdir(parents=True, exist_ok=True)
    state = cert_dir / "state.json"
    if state.exists():
        s = json.loads(state.read_text())
        assert Variant.from_state(s) == var, f"state.json variant {Variant.from_state(s)} != {var}"
        T = s["T"]; br = Brackets3.from_dict(s); hist = s["history"]; ref = s["reference"]
        print(f"resuming {tag}: T={T:.6f} passes done={len(hist)}", flush=True)
    else:
        ref, T, ball = reference_point(fr, bao, sn, spec, args.Delta)
        br = initial_brackets3(fr); hist = []
        print(f"{tag} [{var.describe()}; n_SN={len(sn.m)}, n_BAO={len(bao.value)}]: class-min chi2 {ball.str(12)} -> "
              f"T={T:.6f}; nodes={spec.n_seg+1}", flush=True)
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
                                         history=hist, H0_max=H0, lambda0_min=lb, **var.as_dict(),
                                         n_sn=len(sn.m), n_bao=len(bao.value), **br.to_dict()), indent=1))
        if last is not None and abs(lb - last) < args.tol:
            break
        last = lb
    print("TOTAL done", flush=True)


if __name__ == "__main__":
    main()
