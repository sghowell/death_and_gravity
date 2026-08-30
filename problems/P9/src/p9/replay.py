"""Offline replay of a certified LKR chain from the class-box base.

Independent of the solver: only the stored dual vectors are consumed. For each pass, the relaxed
program is rebuilt from the frozen data and the *replayed* brackets (never the stored ones), every
stored dual is re-verified in ball arithmetic (Verifier3), brackets are recomputed with outward
rounding, and finally the last pass's bound certificate is re-verified. The replayed final bound is
reported next to the stored one; the replay never uses a number from state.json except T, the
reference point (re-checked: exact class membership + rigorous chi2), the pass count and the variant
(SN sample, D_V row, r_d box: these define the frozen inputs and the class, FORMULATION §6.1).

usage: PYTHONPATH=problems/P9/src uv run python -m p9.replay results/certificates/lkr_L1.5_D4_r2 [--workers 8]
       [--sn ... --dv --rd_box ...]   (optional: must agree with the chain's state.json)
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import re
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
from flint import arb

from . import C_KM_S
from .certify_feasible import in_class_exact
from .lkr import initial_brackets3
from .lkr2 import LKRModel2
from .socp2 import MP_BOX
from .variants import Variant, add_variant_args
from .verify import _endpoint, rigorous_chi2
from .verify3 import Verifier3

_V = None
TAG_RE = re.compile(r"p(\d+)_([a-z]+?)(\d+)_(lo|hi)\.npz$")


def _init(fr, br, T):
    global _V
    _V = Verifier3(fr, br, T)


def _verify(args):
    path, _ = args
    d = np.load(path, allow_pickle=True)
    z = d["z"]; qd = {int(v): float(cf) for v, cf in d["q"]}
    nb = _V.lay.nbao; n = _V.lay.n
    n_eq = (len(_V.B.eq_rows) - nb) + nb + n          # eq = BAO(nbao) + SN dense block(n) + rest
    n_le = len(_V.B.le_rows)
    lb = _V.certify(z[:n_eq], z[n_eq:n_eq + n_le], z[n_eq + n_le:], qd, verbose=False)
    return os.path.basename(path), lb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cert_dir"); ap.add_argument("--workers", type=int, default=8)
    add_variant_args(ap)
    a = ap.parse_args()
    cdir = Path(a.cert_dir)
    st = json.loads((cdir / "state.json").read_text())
    var = Variant.from_state(st)
    if a.sn != "pantheon" or a.dv or a.rd_box is not None:
        assert Variant.from_args(a) == var, f"command-line variant {Variant.from_args(a)} != chain's {var}"
    bao, sn, spec, fr = var.frozen(st["L"], st["refine"])
    # reference point and T, re-checked
    u_ref = np.array(st["reference"]["u"]); Mp_ref = st["reference"]["Mp"]
    assert in_class_exact(spec, u_ref) and MP_BOX[0] <= Mp_ref <= MP_BOX[1], "reference point not in class"
    ball = rigorous_chi2(fr, u_ref, Mp_ref)
    T_replay = _endpoint(ball + arb(st["Delta"]), +1)
    assert T_replay <= st["T"] + 1e-12, (T_replay, st["T"])
    T = st["T"]                          # using the stored T is valid since T_replay <= T (F only grows)
    print(f"replay {cdir.name} [{var.describe()}]: reference chi2 {ball.str(12)} -> T_replay={T_replay:.9f} (stored {T:.9f})", flush=True)
    br = initial_brackets3(fr)
    n_pass = len(st["history"])
    from .run_lkr_certified import c_node_balls
    c_ball = c_node_balls(fr)
    ctx = mp.get_context("spawn")
    for p in range(n_pass):
        files = sorted(f for f in cdir.glob(f"p{p}_*.npz") if "final" not in f.name)
        t0 = time.time()
        with ctx.Pool(a.workers, initializer=_init, initargs=(fr, br, T)) as pool:
            res = pool.map(_verify, [(str(f), None) for f in files], chunksize=4)
        out = br.copy()
        for name, lb in res:
            m = TAG_RE.match(name)
            assert m and int(m.group(1)) == p, name
            kind, i, side = m.group(2), int(m.group(3)), m.group(4)
            out.apply_bound(kind, i, side, lb, c_ball[i] if kind == "rho" else 0.0)
        br = out
        # final bound of this pass
        fin = cdir / f"p{p}_final_lambda0.npz"
        m = LKRModel2(fr, br, T); ver = Verifier3(fr, br, T)
        d = np.load(fin); z = d["z"]
        qd = {int(m.lay.lam[0]): 1.0}
        lb = ver.certify(z[:m._n_eq], z[m._n_eq:m._n_eq + m._n_in], z[m._n_eq + m._n_in:], qd, verbose=False)
        print(f"  pass {p}: {len(files)} certificates re-verified in {time.time()-t0:.0f}s; replayed bound lambda0 >= {lb:.9f} "
              f"-> H0 <= {C_KM_S/(spec.r_lo*10**lb):.6f} (stored pass value {st['history'][p]['lambda0_min']:.9f})", flush=True)
    print(f"REPLAY RESULT: H0 <= {C_KM_S/(spec.r_lo*10**lb):.6f} (stored {st['H0_max']:.6f}) at r_lo = {spec.r_lo} Mpc", flush=True)


if __name__ == "__main__":
    main()
