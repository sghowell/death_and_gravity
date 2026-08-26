"""Production driver: certified H0_max(L, Delta) with stored certificates, plus verified-feasible
high-H0 class members (lower bounds on the true maximum). Writes results/certified_curve.json.

usage: PYTHONPATH=problems/P9/src uv run python -m p9.run_certified [--L 1.5 3.0] [--Delta 4 1 9] [--subset 300]
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from . import C_KM_S
from .bound import sn_subset
from .certified import certified_bound
from .classmin import minimize_chi2_over_class
from .data import load_desi, load_pantheon, verify_manifest
from .feasible import max_H0_point
from .geometry import lcdm_u_nodes
from .lcdm import fit_bao, fit_bao_sn
from .model import ClassSpec, Frozen

RESULTS = Path(__file__).resolve().parents[2] / "results"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=float, nargs="+", default=[1.0, 1.5, 2.0, 3.0, 5.0, 10.0])
    ap.add_argument("--Delta", type=float, nargs="+", default=[4.0, 1.0, 9.0])
    ap.add_argument("--subset", type=int, default=300)
    ap.add_argument("--passes", type=int, default=6)
    ap.add_argument("--out", default="certified_curve.json")
    args = ap.parse_args()

    verify_manifest()
    bao = load_desi(drop_dv=True)
    sn = load_pantheon()
    sub = sn_subset(sn, args.subset)
    (om_b, hrd_b), chi_b = fit_bao(bao)
    (om_j, hrd_j), chi_j = fit_bao_sn(bao, sn)
    print(f"LCDM BAO-only: Omega_m={om_b:.4f} h*r_d={hrd_b:.2f} chi2={chi_b:.2f}; "
          f"BAO+SN: Omega_m={om_j:.4f} h*r_d={hrd_j:.2f} chi2={chi_j:.2f} H0(147.09)={100*hrd_j/147.09:.2f}", flush=True)
    out = RESULTS / args.out
    rows = json.loads(out.read_text())["rows"] if out.exists() else []
    done = {(r["L"], r["Delta"]) for r in rows}
    for L in args.L:
        spec = ClassSpec(L=L, grid_kind="geometric")
        fr = Frozen(bao, sn, spec); frs = Frozen(bao, sub, spec)
        u_lcdm = lcdm_u_nodes(spec.x, om_j, hrd_j)
        # class minimum (two starts)
        cands = [minimize_chi2_over_class(fr, u_lcdm), minimize_chi2_over_class(fr, np.full(spec.n_seg + 1, u_lcdm[0]))]
        u_star, Mp_star, chi_star = min(cands, key=lambda c: c[2])
        lcdm_in_class = float(fr.chi2(u_lcdm, fr.best_Mp(u_lcdm))) if np.all(np.abs(np.diff(u_lcdm)) <= L * spec.hs * np.minimum(u_lcdm[:-1], u_lcdm[1:])) else None
        for Delta in args.Delta:
            if (L, Delta) in done:
                continue
            tag = f"L{L:g}_D{Delta:g}"
            print(f"\n== {tag}: class-min chi2={chi_star:.4f} (LCDM chi2 in class: {lcdm_in_class}) ==", flush=True)
            t0 = time.time()
            res = certified_bound(fr, frs, Delta, u_star, tag, n_passes=args.passes)
            u_f, Mp_f, chi2_f, H0_f = max_H0_point(fr, res["T"], np.asarray(res["u_ref"]))
            row = dict(L=L, Delta=Delta, chi2_class_min=chi_star, chi2_lcdm=lcdm_in_class, T=res["T"],
                       H0_max_certified=res["H0_max"], u0_min_certified=res["u0_min"],
                       H0_feasible=H0_f, chi2_feasible=chi2_f, u_feasible=u_f.tolist(),
                       H0_ref=res["H0_ref"], seconds=time.time() - t0, passes=len(res["history"]),
                       history=res["history"], certificates=f"certificates/{tag}")
            rows.append(row)
            print(f"-> {tag}: certified H0 <= {res['H0_max']:.3f}; feasible H0 >= {H0_f}; "
                  f"[{time.time()-t0:.0f}s]", flush=True)
            out.write_text(json.dumps(dict(lcdm=dict(bao_only=dict(om=om_b, h_rd=hrd_b, chi2=chi_b),
                                                     bao_sn=dict(om=om_j, h_rd=hrd_j, chi2=chi_j)),
                                           r_lo=spec.r_lo, subset=args.subset, rows=rows), indent=1))
    print("\nL     Delta  chi2_min   H0_max(cert)  H0_feasible")
    for r in sorted(rows, key=lambda r: (r["L"], r["Delta"])):
        print(f"{r['L']:<5g} {r['Delta']:<6g} {r['chi2_class_min']:9.3f}  {r['H0_max_certified']:8.3f}     "
              f"{r['H0_feasible'] if r['H0_feasible'] is None else round(r['H0_feasible'], 3)}")


if __name__ == "__main__":
    main()
