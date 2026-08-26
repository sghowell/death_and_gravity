"""Rigorous certificate for a feasible point: proves  max_{F} H0 >= H0(u_f)  by exhibiting u_f.

The certificate is the point itself. Verification: (i) u_f in the class (exact rational check of the
slope/box conditions), (ii) an enclosure of chi2(u_f, Mp_f) in ball arithmetic whose upper endpoint is
<= T := upper(chi2(u*, Mp*)) + Delta for the recorded reference point u*. Everything is exact/rigorous.

usage: PYTHONPATH=problems/P9/src uv run python -m p9.certify_feasible --L 1.5 --Delta 4
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

from . import C_KM_S
from .classmin import minimize_chi2_over_class
from .data import load_desi, load_pantheon, verify_manifest
from .feasible import max_H0_point
from .geometry import lcdm_u_nodes
from .lcdm import fit_bao_sn
from .model import ClassSpec, Frozen
from .verify import _endpoint, rigorous_chi2

RESULTS = Path(__file__).resolve().parents[2] / "results"


def in_class_exact(spec: ClassSpec, u: np.ndarray) -> bool:
    """Exact rational check of the class constraints (grid nodes and L as the recorded floats)."""
    x = [Fraction(float(v)) for v in spec.x]
    L = Fraction(float(spec.L))
    uu = [Fraction(float(v)) for v in u]
    ulo = Fraction(C_KM_S) / (Fraction(float(spec.H_max)) * Fraction(float(spec.r_hi)))
    uhi = Fraction(C_KM_S) / (Fraction(float(spec.H_min)) * Fraction(float(spec.r_lo)))
    if any(v < ulo or v > uhi for v in uu):
        return False
    for k in range(len(x) - 1):
        h = x[k + 1] - x[k]
        if abs(uu[k + 1] - uu[k]) > L * h * min(uu[k], uu[k + 1]):
            return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=float, default=1.5)
    ap.add_argument("--Delta", type=float, default=4.0)
    args = ap.parse_args()
    verify_manifest()
    bao = load_desi(); sn = load_pantheon()
    spec = ClassSpec(L=args.L, grid_kind="geometric"); fr = Frozen(bao, sn, spec)
    (om, hrd), _ = fit_bao_sn(bao, sn)
    u = lcdm_u_nodes(spec.x, om, hrd)
    cands = [minimize_chi2_over_class(fr, u), minimize_chi2_over_class(fr, np.full(spec.n_seg + 1, u[0]))]
    u_star, Mp_star, chi2_star = min(cands, key=lambda c: c[2])
    assert in_class_exact(spec, u_star)
    ball_star = rigorous_chi2(fr, u_star, Mp_star)
    T = _endpoint(ball_star, +1) + args.Delta
    u_f, Mp_f, chi2_f, H0_f = max_H0_point(fr, T, u_star, verbose=False)
    ok_class = in_class_exact(spec, u_f)
    ball_f = rigorous_chi2(fr, u_f, Mp_f)
    chi2_f_up = _endpoint(ball_f, +1)
    ok = bool(ok_class) and bool(chi2_f_up <= T)
    H0 = C_KM_S / (spec.r_lo * u_f[0])
    print(f"L={args.L} Delta={args.Delta}: reference chi2 in {ball_star.str(12)} -> T={T:.6f}")
    print(f"  feasible point: in class (exact) = {ok_class}; chi2 in {ball_f.str(12)}; upper {chi2_f_up:.6f} <= T: {chi2_f_up <= T}")
    print(f"  => CERTIFIED lower bound on max H0 over F: {H0:.4f} (at r_lo={spec.r_lo})" if ok else "  => NOT certified")
    out = RESULTS / "certificates" / f"feasible_L{args.L:g}_D{args.Delta:g}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(dict(L=args.L, Delta=args.Delta, grid_x=[float(v) for v in spec.x], r_lo=spec.r_lo,
                                   u_ref=u_star.tolist(), Mp_ref=Mp_star, chi2_ref_enclosure=ball_star.str(20),
                                   T=T, u_feasible=u_f.tolist(), Mp_feasible=Mp_f, chi2_feasible_enclosure=ball_f.str(20),
                                   in_class_exact=bool(ok_class), certified=bool(ok), H0_lower_bound=(float(H0) if ok else None)), indent=1))
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
