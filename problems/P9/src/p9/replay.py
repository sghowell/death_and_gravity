"""Offline replay of a certified LKR chain from the class-box base.

Independent of the solver: only the stored dual vectors are consumed. For each pass, the relaxed
program is rebuilt from the frozen data and the *replayed* brackets (never the stored ones), every
stored dual is re-verified in ball arithmetic (Verifier3), brackets are recomputed with outward
rounding, and finally the last pass's bound certificate is re-verified. The replayed final bound is
reported next to the stored one; the replay never uses a number from state.json except T, the
reference point (re-checked: exact class membership + rigorous chi2) and the pass count.

usage: PYTHONPATH=problems/P9/src uv run python -m p9.replay results/certificates/lkr_L1.5_D4_r2 [--workers 8]
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
from flint import arb

from . import C_KM_S
from .certify_feasible import in_class_exact
from .data import load_desi, load_pantheon, verify_manifest
from .lkr import Brackets3, initial_brackets3
from .lkr2 import LKRModel2
from .lkr_rows import ArbArith
from .model import ClassSpec, Frozen
from .socp2 import MP_BOX
from .verify import _endpoint, rigorous_chi2
from .verify3 import Verifier3

_V = None; _LAY = None


def _init(fr, br, T):
    global _V, _LAY
    _V = Verifier3(fr, br, T)
    _LAY = LKRModel2(fr, br, T).lay


def _verify(args):
    path, expect_q = args
    d = np.load(path, allow_pickle=True)
    z = d["z"]; q_items = d["q"]
    qd = {int(v): float(cf) for v, cf in q_items}
    n_eq = len(_V.B.eq_rows) + _V.lay.n - 0  # eq rows count = bao rows + SN block (n) + rest
    # recompute the split from the layout: eq = BAO(nbao) + SN(n) + rest(len(eq_rows) - nbao)
    nb = _V.lay.nbao; n = _V.lay.n
    n_eq = (len(_V.B.eq_rows) - nb) + nb + n
    n_le = len(_V.B.le_rows)
    lb = _V.certify(z[:n_eq], z[n_eq:n_eq + n_le], z[n_eq + n_le:], qd, verbose=False)
    return os.path.basename(path), lb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cert_dir"); ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    cdir = Path(a.cert_dir)
    st = json.loads((cdir / "state.json").read_text())
    verify_manifest()
    bao = load_desi(); sn = load_pantheon()
    spec = ClassSpec(L=st["L"], grid_kind="geometric", refine=st["refine"]); fr = Frozen(bao, sn, spec)
    # reference point and T, re-checked
    u_ref = np.array(st["reference"]["u"]); Mp_ref = st["reference"]["Mp"]
    assert in_class_exact(spec, u_ref) and MP_BOX[0] <= Mp_ref <= MP_BOX[1], "reference point not in class"
    ball = rigorous_chi2(fr, u_ref, Mp_ref)
    T_replay = _endpoint(ball + arb(st["Delta"]), +1)
    assert T_replay <= st["T"] + 1e-12, (T_replay, st["T"])
    T = st["T"]                          # using the stored T is valid since T_replay <= T (F only grows)
    print(f"replay {cdir.name}: reference chi2 {ball.str(12)} -> T_replay={T_replay:.9f} (stored {T:.9f})", flush=True)
    br = initial_brackets3(fr)
    n_pass = len(st["history"])
    ar = ArbArith()
    c_ball = [arb(0)] + [ar.log10(ar.expm1(ar.c(v))) for v in spec.x[1:]]
    N = spec.n_seg
    ctx = mp.get_context("spawn")
    for p in range(n_pass):
        files = sorted(cdir.glob(f"p{p}_*.npz"))
        files = [f for f in files if "final" not in f.name]
        t0 = time.time()
        with ctx.Pool(a.workers, initializer=_init, initargs=(fr, br, T)) as pool:
            res = pool.map(_verify, [(str(f), None) for f in files], chunksize=4)
        rho_lo = br.rho_lo.copy(); rho_hi = br.rho_hi.copy(); lam_lo = br.lam_lo.copy(); lam_hi = br.lam_hi.copy()
        yb_lo = br.yb_lo.copy(); yb_hi = br.yb_hi.copy()
        for name, lb in res:
            kind, side = name.split("_")[1], name.split("_")[2].split(".")[0]
            i = int(kind[3:]) if kind.startswith(("rho", "lam")) else int(kind[2:])
            if kind.startswith("rho"):
                if side == "lo": rho_lo[i] = max(rho_lo[i], _endpoint(arb(lb) - c_ball[i], -1))
                else: rho_hi[i] = min(rho_hi[i], _endpoint(arb(-lb) - c_ball[i], +1))
            elif kind.startswith("lam"):
                if side == "lo": lam_lo[i] = max(lam_lo[i], lb)
                else: lam_hi[i] = min(lam_hi[i], -lb)
            else:
                if side == "lo": yb_lo[i] = max(yb_lo[i], lb)
                else: yb_hi[i] = min(yb_hi[i], -lb)
        br = Brackets3(rho_lo, rho_hi, yb_lo, yb_hi, lam_lo, lam_hi)
        # final bound of this pass
        fin = cdir / f"p{p}_final_lambda0.npz"
        m = LKRModel2(fr, br, T); ver = Verifier3(fr, br, T)
        d = np.load(fin); z = d["z"]
        qd = {int(m.lay.lam[0]): 1.0}
        lb = ver.certify(z[:m._n_eq], z[m._n_eq:m._n_eq + m._n_in], z[m._n_eq + m._n_in:], qd, verbose=False)
        print(f"  pass {p}: {len(files)} certificates re-verified in {time.time()-t0:.0f}s; replayed bound lambda0 >= {lb:.9f} "
              f"-> H0 <= {C_KM_S/(spec.r_lo*10**lb):.6f} (stored pass value {st['history'][p]['lambda0_min']:.9f})", flush=True)
    print(f"REPLAY RESULT: H0 <= {C_KM_S/(spec.r_lo*10**lb):.6f} (stored {st['H0_max']:.6f})", flush=True)


if __name__ == "__main__":
    main()
