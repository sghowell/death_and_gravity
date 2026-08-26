"""LKR bound driver with parallel (spawn) tightening.

usage: PYTHONPATH=problems/P9/src uv run python -m p9.run_lkr --L 1.5 --Delta 4 --workers 6 --passes 12
Resumes from results/lkr_L<L>_D<Delta>.json if present.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")

import numpy as np

from . import C_KM_S
from .classmin import minimize_chi2_over_class
from .data import load_desi, load_pantheon
from .geometry import lcdm_u_nodes
from .lcdm import fit_bao_sn
from .lkr import Brackets3, LKRModel, initial_brackets3
from .model import ClassSpec, Frozen
from .parallel import tighten_parallel

RESULTS = Path(__file__).resolve().parents[2] / "results"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=float, default=1.5)
    ap.add_argument("--Delta", type=float, default=4.0)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--passes", type=int, default=12)
    ap.add_argument("--tol", type=float, default=2e-5)
    ap.add_argument("--refine", type=int, default=0)
    args = ap.parse_args()
    bao = load_desi(); sn = load_pantheon()
    spec = ClassSpec(L=args.L, grid_kind="geometric", refine=args.refine); fr = Frozen(bao, sn, spec)
    out = RESULTS / f"lkr_L{args.L:g}_D{args.Delta:g}_r{args.refine}.json"
    if out.exists():
        r = json.loads(out.read_text())
        T = r["T"]; chi2_min = r.get("chi2_class_min")
        br = Brackets3(*(np.array(r[k]) for k in ["rho_lo", "rho_hi", "yb_lo", "yb_hi", "lam_lo", "lam_hi"]))
        hist = r.get("history", [])
        print(f"resuming from {out}: T={T:.4f} widths {br.width()}", flush=True)
    else:
        (om, hrd), _ = fit_bao_sn(bao, sn)
        u = lcdm_u_nodes(spec.x, om, hrd)
        cands = [minimize_chi2_over_class(fr, u), minimize_chi2_over_class(fr, np.full(spec.n_seg + 1, u[0]))]
        u_star, Mp_star, chi2_min = min(cands, key=lambda c: c[2])
        T = chi2_min + args.Delta
        br = initial_brackets3(fr); hist = []
        print(f"class-min chi2={chi2_min:.4f}  T={T:.4f}  initial widths {br.width()}", flush=True)
    last = hist[-1]["lambda0_min"] if hist else None
    for it in range(len(hist), len(hist) + args.passes):
        t = time.time()
        br = tighten_parallel(fr, br, T, n_workers=args.workers)
        m = LKRModel(fr, br, T); c2, _ = m.min_chi2(); l0 = m.min_lambda0(); H0 = C_KM_S / (spec.r_lo * 10 ** l0)
        hist.append(dict(pass_=it, relaxed_min=c2, lambda0_min=l0, H0_max=H0, widths=br.width(), seconds=time.time() - t))
        print(f"pass {it}: relaxed min chi2={c2:.3f}  lambda0_min={l0:.6f} -> H0_max={H0:.4f}  widths {br.width()} [{time.time()-t:.0f}s]", flush=True)
        out.write_text(json.dumps(dict(L=args.L, Delta=args.Delta, T=T, chi2_class_min=chi2_min, H0_max=H0, lambda0_min=l0,
                                       relaxed_min=c2, history=hist,
                                       **{k: getattr(br, k).tolist() for k in ["rho_lo", "rho_hi", "yb_lo", "yb_hi", "lam_lo", "lam_hi"]}), indent=1))
        if last is not None and abs(l0 - last) < args.tol:
            break
        last = l0
    print("TOTAL done", flush=True)


if __name__ == "__main__":
    main()
