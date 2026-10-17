"""S3 review: independent (floating-point) argument-principle counts and per-segment increments.

Uses only the reviewer's W (indep_w.py).  For each l in the list:
  * every stored segment of certificates/part3/winding_PB_l{l}.json is re-traced with a dense,
    adaptively refined grid (consecutive principal-argument steps < 0.2 rad, and the n / 2n-1 grid
    sums agreeing to 1e-9), and the float increment is compared with the stored ball;
  * the LEFT half of the boundary of R'_l, (0, e) -> (-a, e) -> (-a, -b) -> (0, -b), is traced
    the same way (no symmetry used), so the full-boundary winding is computed without the
    conjugation symmetry;
  * the full-boundary winding (sum over both halves)/(2 pi) is reported.
Not rigorous: a cross-check by an independent code and formulation.
Usage: python indep_winding.py 0-19     (writes indep_winding.json)
"""
from __future__ import annotations

import json
import math
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from indep_w import Geo, W_l  # noqa: E402

HERE = Path(__file__).resolve().parent
CERT = HERE.parent / "part3"
CHUNK = 400


def ball_mid_rad(s: str):
    s = s.strip()[1:-1]
    if "+/-" in s:
        m, r = s.split("+/-")
        m = m.strip()
        return (float(m) if m else 0.0), float(r)
    return float(s), 0.0


def eval_W(g, l, pts: np.ndarray) -> np.ndarray:
    """W at many points, in chunks of similar |w| (keeps the RMS error control honest)."""
    order = np.argsort(np.abs(pts))
    out = np.empty(len(pts), dtype=complex)
    for i in range(0, len(pts), CHUNK):
        idx = order[i:i + CHUNK]
        out[idx] = W_l(g, l, pts[idx], rtol=1e-13)
    return out


def trace(g, l, paths, n0, max_iter=8):
    """paths: list of (z0, z1) complex. Returns per-path increments, max step, points used."""
    ns = [max(8, int(math.ceil(abs(z1 - z0) / n0))) for z0, z1 in paths]
    prev = [None] * len(paths)
    done = [False] * len(paths)
    res = [None] * len(paths)
    for it in range(max_iter):
        todo = [k for k in range(len(paths)) if not done[k]]
        if not todo:
            break
        grids = []
        for k in todo:
            z0, z1 = paths[k]
            grids.append(z0 + (z1 - z0) * np.linspace(0.0, 1.0, ns[k] + 1))
        allp = np.concatenate(grids)
        Wv = eval_W(g, l, allp)
        pos = 0
        for k, gr in zip(todo, grids):
            seg = Wv[pos:pos + len(gr)]
            pos += len(gr)
            steps = np.angle(seg[1:] / seg[:-1])
            inc = float(np.sum(steps))
            mx = float(np.max(np.abs(steps)))
            ok = mx < 0.2 and prev[k] is not None and abs(prev[k] - inc) < 1e-9
            res[k] = {"inc": inc, "max_step": mx, "n": ns[k], "min_absW": float(np.min(np.abs(seg)))}
            if ok:
                done[k] = True
            else:
                prev[k] = inc
                ns[k] = 2 * ns[k]
    for k in range(len(paths)):
        res[k]["converged"] = done[k]
    return res


def run_l(g, l):
    t0 = time.time()
    st = json.loads((CERT / f"winding_PB_l{l}.json").read_text())
    reg = st["region"]
    a, b, e = (float(Fr(reg[k])) for k in ("a", "b", "e"))
    pieces = st["pieces"]
    paths = [(complex(float(Fr(r["p0"][0])), float(Fr(r["p0"][1]))),
              complex(float(Fr(r["p1"][0])), float(Fr(r["p1"][1])))) for r in pieces]
    h0 = g.km / 32
    seg = trace(g, l, paths, h0)
    devs = []
    worst = None
    n_in = 0
    for r, s in zip(pieces, seg):
        mid, rad = ball_mid_rad(r["increment"])
        dv = abs(s["inc"] - mid)
        inside = dv <= rad + 1e-8
        n_in += inside
        devs.append(dv)
        q = dv / rad if rad > 0 else float("inf")
        if worst is None or q > worst[0]:
            worst = (q, r["key"], s["inc"], r["increment"])
    # left half, (0, e) -> (-a, e) -> (-a, -b) -> (0, -b), cut into pieces of length <= kappa_-
    corners = [complex(0, e), complex(-a, e), complex(-a, -b), complex(0, -b)]
    lpaths = []
    for z0, z1 in zip(corners[:-1], corners[1:]):
        n = max(1, int(math.ceil(abs(z1 - z0) / g.km)))
        for j in range(n):
            lpaths.append((z0 + (z1 - z0) * j / n, z0 + (z1 - z0) * (j + 1) / n))
    left = trace(g, l, lpaths, h0)
    sum_right = sum(s["inc"] for s in seg)
    sum_left = sum(s["inc"] for s in left)
    total = sum_right + sum_left
    return {"l": l, "a_over_km": a / g.km, "segments": len(pieces),
            "all_segments_converged": all(s["converged"] for s in seg),
            "left_converged": all(s["converged"] for s in left),
            "segments_inside_stored_ball": n_in,
            "max_abs_dev_from_stored_mid": max(devs),
            "worst_dev_over_radius": {"ratio": worst[0], "key": worst[1], "float_inc": worst[2], "stored": worst[3]},
            "sum_right_over_pi": sum_right / math.pi, "sum_left_over_pi": sum_left / math.pi,
            "winding_full_boundary": total / (2 * math.pi),
            "min_absW_right_float": min(s["min_absW"] for s in seg),
            "min_absW_left_float": min(s["min_absW"] for s in left),
            "max_points_per_segment": max(s["n"] for s in seg),
            "wall_s": round(time.time() - t0, 1)}


def main():
    spec = sys.argv[1] if len(sys.argv) > 1 else "0-19"
    ls = []
    for part in spec.split(","):
        if "-" in part:
            u, v = part.split("-")
            ls += list(range(int(u), int(v) + 1))
        else:
            ls.append(int(part))
    g = Geo()
    outp = HERE / "indep_winding.json"
    out = json.loads(outp.read_text()) if outp.exists() else {"per_l": {}}
    for l in ls:
        r = run_l(g, l)
        out["per_l"][str(l)] = r
        print(json.dumps(r), flush=True)
        outp.write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
