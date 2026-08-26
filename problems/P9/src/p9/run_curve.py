"""S1 driver: known-answer tests, then the curve H0_max(L, Delta). Writes results/curve.json."""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

from . import C_KM_S
from .bound import compute_bound, sn_subset
from .data import load_desi, load_pantheon, verify_manifest
from .geometry import lcdm_u_nodes
from .lcdm import fit_bao, fit_bao_sn
from .model import ClassSpec, Frozen

RESULTS = Path(__file__).resolve().parents[2] / "results"


def lcdm_max_logslope(om: float, z_max: float) -> float:
    z = z_max
    return 1.5 * om * (1 + z) ** 3 / (om * (1 + z) ** 3 + 1 - om)


def main(Ls=(0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0), Deltas=(9.0, 4.0, 1.0), N=None,
         subset_size=300, out_name="curve.json"):
    verify_manifest()
    bao = load_desi(drop_dv=True)
    sn = load_pantheon()
    sub = sn_subset(sn, subset_size)
    print(f"BAO rows: {len(bao.z)}  SN: {len(sn.m)}  bracket subset: {len(sub.m)}", flush=True)

    (om_b, hrd_b), chi_b = fit_bao(bao)
    print(f"LCDM BAO-only: Omega_m={om_b:.4f} h*r_d={hrd_b:.2f} Mpc chi2={chi_b:.2f} "
          f"(DESI DR2 published: 0.2975+-0.0086, 101.54+-0.73)", flush=True)
    (om_j, hrd_j), chi_j = fit_bao_sn(bao, sn)
    print(f"LCDM BAO+SN: Omega_m={om_j:.4f} h*r_d={hrd_j:.2f} Mpc chi2={chi_j:.2f} "
          f"-> H0(r_d=147.09)={100*hrd_j/147.09:.2f}, beta=c/(r_d H0)={C_KM_S/(100*hrd_j):.2f}", flush=True)

    rows = []
    Deltas = sorted(Deltas, reverse=True)   # largest first: its node bounds are valid seeds for smaller Delta
    for L in Ls:
        spec = ClassSpec(L=L, grid_kind="geometric")
        fr = Frozen(bao, sn, spec)
        frs = Frozen(bao, sub, spec)
        u_ref = lcdm_u_nodes(spec.x, om_j, hrd_j)
        if L < lcdm_max_logslope(om_j, spec.z_max):
            u_ref = np.full(spec.n_seg + 1, u_ref[0])       # LCDM not in the class; a constant-u point is
        nb = None
        for Delta in Deltas:
            t0 = time.time()
            print(f"\n== L={L} Delta={Delta} nodes={spec.n_seg + 1} ==", flush=True)
            res = compute_bound(fr, frs, Delta, u_ref, nb=nb)
            nb = res.node_bounds
            u_ref = res.u_ref
            rows.append(dict(L=L, Delta=Delta, nodes=spec.n_seg + 1, T=res.T, chi2_ref=res.chi2_ref,
                             u0_min=res.u0_min, H0_max=res.H0_max, H0_ref=res.H0_ref,
                             seconds=time.time() - t0, history=res.history,
                             u_ref=res.u_ref.tolist(), Mp_ref=res.Mp_ref,
                             node_bounds=dict(u_lo=nb.u_lo.tolist(), u_hi=nb.u_hi.tolist(),
                                              dm_lo=nb.dm_lo.tolist(), dm_hi=nb.dm_hi.tolist())))
            print(f"-> H0_max = {res.H0_max:.3f} km/s/Mpc  (ref point H0={res.H0_ref:.3f}, "
                  f"chi2_ref={res.chi2_ref:.2f}) [{time.time()-t0:.0f}s]", flush=True)
            RESULTS.mkdir(exist_ok=True)
            (RESULTS / out_name).write_text(json.dumps(dict(
                lcdm=dict(bao_only=dict(om=om_b, h_rd=hrd_b, chi2=chi_b), bao_sn=dict(om=om_j, h_rd=hrd_j, chi2=chi_j)),
                subset_size=subset_size, rows=rows), indent=1))
    print("\nL      Delta   H0_max   H0_ref")
    for r in rows:
        print(f"{str(r['L']):6} {r['Delta']:5.1f}  {r['H0_max']:.3f}  {r['H0_ref']:.3f}")


if __name__ == "__main__":
    args = sys.argv[1:]
    if args == ["smoke"]:
        main(Ls=(1.5,), Deltas=(4.0,), out_name="smoke.json")
    else:
        main()
