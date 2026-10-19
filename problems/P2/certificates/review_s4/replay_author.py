"""Exact replays with the author's code at the clean review head (S4 review).

  * the three rows stamped 35201758+dirty in the author's working files: S+ leaf [22/25, 9/10] and the L boxes
    [159/20, 8] x [-1/2, -9/20] and [159/20, 8] x [-11/20, -1/2];
  * the thin S leaves [0.82, 0.84], [0.84, 0.86], [0.86, 0.88], [0.88, 0.9] on both sides and both closing boxes,
    recomputed from their stored designs and targets;
  * every T zeta-box and the A box recomputed in full from the stored partitions (not a 5% sample).
Every stored field that the computation produces is compared string for string.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

W = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(W / "src"))

from p2 import lowfreq as LF  # noqa: E402
from p2 import part1_T as TT  # noqa: E402
from p2 import throat as TH  # noqa: E402
from p2.arbutil import git_commit, precision  # noqa: E402
from p2.params import Point  # noqa: E402

CERT = W / "certificates" / "part1"
WB = Fraction(99, 12500)


def norm(x):
    return json.loads(json.dumps(x))


def cmp_zleaf(new, old):
    new = norm(new)
    keys = ["design", "X_a", "k", "exit", "w1_check", "target", "eps_max"]
    diffs = [k for k in keys if new.get(k) != old.get(k)]
    for a, b in zip(new["eps_boxes"], old["eps_boxes"]):
        for k in ("eps", "mu_cells", "mu_tail", "tail_disc", "worst_cell", "cells", "depths", "n_fails", "pass"):
            if a.get(k) != b.get(k):
                diffs.append(f"eps_boxes.{k}")
    if len(new["eps_boxes"]) != len(old["eps_boxes"]):
        diffs.append("n_eps_boxes")
    if new["pass"] != old["pass"]:
        diffs.append("pass")
    return diffs


def cmp_closing(new, old):
    new = norm(new)
    keys = ["design", "near_horizon", "exit", "w1_check", "target", "mu_cells", "mu_tail", "worst_cell", "cells", "n_fails", "pass", "eps_max"]
    return [k for k in keys if new.get(k) != old.get(k)]


def main():
    pt = Point("PB", 256)
    g = TH.FGeo(pt)
    out = {"commit": git_commit(), "S": [], "L": [], "T": [], "A": None}
    t0 = time.process_time()
    with precision(128):
        for side, tag, picks in ((1, "plus", [("22/25", "9/10"), ("41/50", "21/25"), ("21/25", "43/50"), ("43/50", "22/25")]),
                                 (-1, "c", [("41/50", "21/25"), ("21/25", "43/50"), ("43/50", "22/25"), ("22/25", "9/10")])):
            cert = json.loads((CERT / f"S_PB_{tag}.json").read_text())
            for x in cert["leaves"]:
                if x.get("closing"):
                    r = TH.certify_closing(pt, g, side, WB, design={"design": x["design"], "target": x["target"]})
                    d = cmp_closing(r, x)
                    out["S"].append({"side": tag, "zeta": x["zeta_over_OmegaP"], "closing": True, "diffs": d,
                                     "mu_cells": r["mu_cells"], "stored_mu_cells": x["mu_cells"]})
                    print(tag, "closing", d, r["mu_cells"], flush=True)
                elif tuple(x["zeta_over_OmegaP"]) in picks:
                    zl, zh = (Fraction(v) for v in x["zeta_over_OmegaP"])
                    r = TH.certify_zbox(pt, g, side, zl, zh, WB,
                                        design={"design": x["design"], "X_a": x["X_a"], "k": x["k"], "target": x["target"]})
                    d = cmp_zleaf(r, x)
                    out["S"].append({"side": tag, "zeta": x["zeta_over_OmegaP"], "diffs": d,
                                     "mu_cells": [e["mu_cells"] for e in r["eps_boxes"]],
                                     "stored_mu_cells": [e["mu_cells"] for e in x["eps_boxes"]]})
                    print(tag, x["zeta_over_OmegaP"], d, flush=True)
        lp = json.loads((CERT / "L_PB_plus.json").read_text())["leaves_detail"]
        lc = json.loads((CERT / "L_PB_c.json").read_text())["leaves_detail"]
        ch = LF.LChecker(pt, WB)
        for box in (["159/20", "8", "-1/2", "-9/20"], ["159/20", "8", "-11/20", "-1/2"], ["0", "1/20", "-1/20", "0"],
                    ["0", "1/20", "-11/20", "-1/2"]):
            i = [j for j, y in enumerate(lp) if y["box"] == box][0]
            r = LF.LChecker.box(ch, tuple(Fraction(v) for v in box))
            d = []
            for name, st in (("plus", lp[i]), ("c", lc[i])):
                for k in ("tau_a", "Re_y_bound", "y_sign_ok", "exit_rule_ok", "terms"):
                    if r[name].get(k) != st.get(k):
                        d.append(f"{name}.{k}")
            if r["M"] != lp[i]["M"]:
                d.append("M")
            out["L"].append({"box": box, "diffs": d, "Re_y_plus": r["plus"]["Re_y_bound"], "Re_y_c": r["c"]["Re_y_bound"]})
            print("L", box, d, flush=True)
        _, rays, tails = TT.setup(pt)
        for name, design in (("T", TT.DESIGN_T), ("A", TT.DESIGN_A)):
            cert = json.loads((CERT / f"{name}_PB.json").read_text())
            ch = TT.TChecker(pt, rays, tails, design)
            for bx in cert["boxes"]:
                zb = tuple(Fraction(v) for v in bx["zeta_over_OmegaP"])
                r = norm(TT.certify_box(ch, pt, zb, design, leaves=bx["leaves"]))
                d = [k for k in ("mu_cert", "E0_low_min", "c_min", "tail_left", "tail_right", "worst_cell", "n_fails", "pass",
                                 "jost_decay") if r[k] != bx[k]]
                rec = {"zeta": bx["zeta_over_OmegaP"], "diffs": d, "mu_cert": r["mu_cert"], "stored": bx["mu_cert"]}
                if name == "T":
                    out["T"].append(rec)
                else:
                    out["A"] = rec
                print(name, bx["zeta_over_OmegaP"], d, flush=True)
    out["cpu_s"] = round(time.process_time() - t0, 1)
    out["all_identical"] = all(not x["diffs"] for x in out["S"] + out["L"] + out["T"] + [out["A"]])
    Path(sys.argv[1] if len(sys.argv) > 1 else "replay_author.json").write_text(json.dumps(out, indent=1) + "\n")
    print("all identical:", out["all_identical"], "cpu", out["cpu_s"])


if __name__ == "__main__":
    main()
