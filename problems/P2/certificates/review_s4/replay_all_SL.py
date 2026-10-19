"""Full exact replay of every S zeta-leaf (both sides, from the stored designs and targets) and every L omega-box
with the author's code at the clean review head; every produced field is compared with the stored certificate.
Usage: replay_all_SL.py S|L out.json
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

W = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(W / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from p2 import lowfreq as LF  # noqa: E402
from p2 import throat as TH  # noqa: E402
from p2.arbutil import git_commit, precision  # noqa: E402
from p2.params import Point  # noqa: E402
from replay_author import cmp_closing, cmp_zleaf  # noqa: E402

CERT = W / "certificates" / "part1"
WB = Fraction(99, 12500)


def main():
    which, outp = sys.argv[1], sys.argv[2]
    pt = Point("PB", 256)
    out = {"commit": git_commit(), "rows": []}
    t0 = time.process_time()
    with precision(128):
        if which == "S":
            g = TH.FGeo(pt)
            for side, tag in ((1, "plus"), (-1, "c")):
                cert = json.loads((CERT / f"S_PB_{tag}.json").read_text())
                for x in cert["leaves"]:
                    if x.get("closing"):
                        r = TH.certify_closing(pt, g, side, WB, design={"design": x["design"], "target": x["target"]})
                        d = cmp_closing(r, x)
                    else:
                        zl, zh = (Fraction(v) for v in x["zeta_over_OmegaP"])
                        r = TH.certify_zbox(pt, g, side, zl, zh, WB,
                                            design={"design": x["design"], "X_a": x["X_a"], "k": x["k"], "target": x["target"]})
                        d = cmp_zleaf(r, x)
                    out["rows"].append({"side": tag, "zeta": x["zeta_over_OmegaP"], "diffs": d, "pass": r["pass"]})
                    print(tag, x["zeta_over_OmegaP"], d, r["pass"], flush=True)
        else:
            lp = json.loads((CERT / "L_PB_plus.json").read_text())["leaves_detail"]
            lc = json.loads((CERT / "L_PB_c.json").read_text())["leaves_detail"]
            ch = LF.LChecker(pt, WB)
            for i, (bp, bc) in enumerate(zip(lp, lc)):
                assert bp["box"] == bc["box"]
                r = ch.box(tuple(Fraction(v) for v in bp["box"]))
                d = []
                for name, st in (("plus", bp), ("c", bc)):
                    for k in ("tau_a", "Re_y_bound", "y_sign_ok", "exit_rule_ok", "terms"):
                        if r[name].get(k) != st.get(k):
                            d.append(f"{name}.{k}")
                if r["M"] != bp["M"]:
                    d.append("M")
                out["rows"].append({"box": bp["box"], "diffs": d, "pass": r["pass"]})
                if d or i % 100 == 0:
                    print(i, bp["box"], d, r["pass"], flush=True)
    out["cpu_s"] = round(time.process_time() - t0, 1)
    out["n"] = len(out["rows"])
    out["all_identical"] = all(not x["diffs"] for x in out["rows"])
    out["all_pass"] = all(x["pass"] for x in out["rows"])
    Path(outp).write_text(json.dumps(out, indent=1) + "\n")
    print("n", out["n"], "identical", out["all_identical"], "pass", out["all_pass"], "cpu", out["cpu_s"])


if __name__ == "__main__":
    main()
