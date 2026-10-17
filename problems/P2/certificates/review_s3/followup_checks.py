"""S3 review follow-up (author's fixes e874b5f8, b19065b6): checks of the new committed outputs.

F1  certificates/part3/replay_PB.json (the C7 record):
    - the sampled segments are exactly the deterministic sample of part3.verify
      (k = max(3, ceil(0.05 n)), indices (j*7919 + l*104729) mod n), 132 in all;
    - every "stored" ball and K equal the certificate's piece;
    - for l = 0, 1, 2, 19 every "recomputed" ball equals, string for string, the reviewer's earlier
      full replay (logs/replay_l*.jsonl, K re-derived at db2a37eb);
    - every recomputed ball contains the reviewer's independent float increment
      (indep_winding: re-traced here for the sampled segments).
F2  The new code path: part3.verify(..., workers=2, out=<this dir>) re-run for l = 0, 1, 5, 10, 19;
    its recomputed balls must equal those of replay_PB.json string for string, and its problems list
    must be empty.  (Writes followup_verify_out.json here, never into certificates/.)
F3  certificates/part3/symmetry_PB_l0.json: the left half of the boundary of R'_0 recomputed with
    part3.left_half_sum (nothing written); the sum must equal the committed string, and its /pi
    must contain 1.
F4  The fast box test's numbers: at the corner (a_19, -b), half-width 1e-4 kappa_-, the relative
    distance of the two tested box corners from the centre value (reviewer's float W), against the
    zsol box ball's relative half-width.
Writes followup_checks.json.
"""
from __future__ import annotations

import json
import math
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
from flint import acb, arb, fmpq  # noqa: E402

from p2 import part3, zsol  # noqa: E402
from p2.arbutil import git_commit, precision  # noqa: E402
from p2.params import Point  # noqa: E402

CERT = HERE.parent / "part3"


def f1(g):
    rp = json.loads((CERT / "replay_PB.json").read_text())
    probs = []
    by_l = {}
    for r in rp["recomputed"]:
        by_l.setdefault(r["l"], []).append(r)
    total = 0
    float_dev = []
    cross = 0
    for l in range(20):
        st = json.loads((CERT / f"winding_PB_l{l}.json").read_text())
        pcs = st["pieces"]
        n = len(pcs)
        k = max(3, math.ceil(0.05 * n))
        want = [pcs[j]["key"] for j in sorted({(j * 7919 + l * 104729) % n for j in range(k)})]
        got = [r["key"] for r in by_l.get(l, [])]
        if got != want:
            probs.append(f"l={l}: sampled keys differ from the deterministic sample")
        piece = {p["key"]: p for p in pcs}
        for r in by_l.get(l, []):
            p = piece[r["key"]]
            if r["stored"] != p["increment"] or r["K"] != p["K"]:
                probs.append(f"l={l}: stored ball or K differs at {r['key']}")
            if not r["overlap"]:
                probs.append(f"l={l}: no overlap at {r['key']}")
        total += len(got)
        # string-for-string against the reviewer's earlier full replay
        lf = HERE / "logs" / f"replay_l{l}.jsonl"
        if lf.exists():
            mine = {json.loads(x)["key"]: json.loads(x) for x in lf.read_text().splitlines()}
            for r in by_l.get(l, []):
                if mine[r["key"]]["increment"] != r["recomputed"]:
                    probs.append(f"l={l}: recomputed ball differs from the reviewer's replay at {r['key']}")
                cross += 1
        # independent float increments for the sampled segments
        from indep_winding import trace
        paths = []
        for r in by_l.get(l, []):
            p = piece[r["key"]]
            paths.append((complex(float(Fr(p["p0"][0])), float(Fr(p["p0"][1]))),
                          complex(float(Fr(p["p1"][0])), float(Fr(p["p1"][1])))))
        res = trace(g, l, paths, g.km / 32)
        with precision(128):
            for r, s in zip(by_l.get(l, []), res):
                b = arb(r["recomputed"])
                dv = abs(float(b.mid()) - s["inc"])
                rad = float(b.rad())
                float_dev.append(dv / rad if rad > 0 else float("inf"))
                if not (dv <= rad + 1e-8):
                    probs.append(f"l={l}: float increment outside the recomputed ball at {r['key']}")
    if total != 132 or rp["segments_recomputed"] != 132:
        probs.append(f"total sampled {total}")
    return {"segments": total, "cross_checked_with_reviewer_replay": cross,
            "max_float_dev_over_radius": max(float_dev), "commit_of_record": rp["commit"],
            "record_pass": rp["pass"], "problems": probs, "ok": not probs}


def f2():
    t0 = time.time()
    outp = HERE / "followup_verify_out.json"
    ls = [0, 1, 5, 10, 19]
    probs = part3.verify("PB", ls, workers=2, out=str(outp))
    mine = json.loads(outp.read_text())
    rp = {r["key"] + f"@{r['l']}": r for r in json.loads((CERT / "replay_PB.json").read_text())["recomputed"]}
    diff = [r["key"] for r in mine["recomputed"] if rp[r["key"] + f"@{r['l']}"]["recomputed"] != r["recomputed"]]
    return {"ls": ls, "problems": probs, "segments": len(mine["recomputed"]),
            "identical_to_replay_PB": not diff, "differing": diff,
            "cpu_segments": mine["timing_s"]["cpu_segments"], "wall_s": round(time.time() - t0, 1),
            "ok": (not probs) and (not diff)}


def f3():
    t0 = time.time()
    st = json.loads((CERT / "symmetry_PB_l0.json").read_text())
    lh = part3.left_half_sum("PB", 0, workers=2)
    with precision(part3.PREC):
        over_pi = arb(lh["sum_left"]) / arb.pi()
        contains1 = bool(over_pi.contains(1))
    return {"committed_sum_left": st["left"]["sum_left"], "recomputed_sum_left": lh["sum_left"],
            "identical": lh["sum_left"] == st["left"]["sum_left"], "segments": lh["segments"],
            "left_over_pi": over_pi.str(10), "contains_1": contains1,
            "wall_s": round(time.time() - t0, 1), "ok": lh["sum_left"] == st["left"]["sum_left"] and contains1}


def f4(g):
    from indep_w import W_l
    pt = Point("PB", 256)
    st = json.loads((CERT / "winding_PB_l19.json").read_text())
    a, b = part3._fq(st["region"]["a"]), part3._fq(st["region"]["b"])
    rho = fmpq(round(float(pt.km.mid()) * 1e-4 * 2 ** 70), 2 ** 70)
    with precision(part3.PREC):
        c = acb(arb(a), arb(-b))
        box = acb(c.real + arb(0, 1) * arb(rho), c.imag + arb(0, 1) * arb(rho))
        _, Wb = zsol.wronskian_z(pt, 19, c, box, theta=part3.THETA_BOX)
        Wc = zsol.wronskian_z(pt, 19, c, theta=part3.THETA_POINT)[0]
        relw = float(max(Wb.real.rad(), Wb.imag.rad()) / abs(Wc).mid())
    cf = complex(float(Fr(str(a))), -float(Fr(str(b))))
    rf = float(Fr(str(rho)))
    vals = W_l(g, 19, np.array([cf, cf + rf * (-1 - 1j), cf + rf * (1 + 1j)]), rtol=1e-13)
    dev = [abs(v - vals[0]) / abs(vals[0]) for v in vals[1:]]
    return {"box_rel_halfwidth": relw, "corner_rel_dev": dev,
            "detects_underestimate_factor_above": relw / max(dev)}


def main():
    from indep_w import Geo
    g = Geo()
    out = {"head": git_commit()}
    out["F1_replay_record"] = f1(g)
    print("F1", json.dumps(out["F1_replay_record"]), flush=True)
    out["F4_box_test_numbers"] = f4(g)
    print("F4", json.dumps(out["F4_box_test_numbers"]), flush=True)
    out["F3_symmetry_l0"] = f3()
    print("F3", json.dumps(out["F3_symmetry_l0"]), flush=True)
    out["F2_verify_new_path"] = f2()
    print("F2", json.dumps(out["F2_verify_new_path"]), flush=True)
    out["ok"] = all(out[k].get("ok", True) for k in ("F1_replay_record", "F2_verify_new_path", "F3_symmetry_l0"))
    (HERE / "followup_checks.json").write_text(json.dumps(out, indent=1) + "\n")
    print("ok:", out["ok"])


if __name__ == "__main__":
    main()
