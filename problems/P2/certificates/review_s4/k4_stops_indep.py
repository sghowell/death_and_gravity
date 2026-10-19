"""K4's secant stops in H+ (Amendment 6), checked with the reviewer's independent float code (indep_points.D_stat):
D = |y_+ - y_c|/(|y_+| + |y_c|) at the stop; a zero of W_l has D = 0.  Every 6th of the 61 stops, plus the three
named in Amendment 6."""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from indep_points import D_stat, km  # noqa: E402

W = Path(__file__).resolve().parents[2]
K4 = json.loads((W / "certificates/part1/controls_PB.json").read_text())["K4"]
stops = [r for r in K4["rows"] if r["root_in_H_plus"]]
pick = sorted({0, 11, 20} | {r["i"] for r in stops[::6]})
out = []
for r in stops:
    if r["i"] not in pick:
        continue
    w = complex(float(r["secant_root_over_km"][0]), float(r["secant_root_over_km"][1])) * km
    D, yp, yc, _ = D_stat(r["l"], w)
    out.append({"i": r["i"], "l": r["l"], "regime": r["regime"], "root_over_km": r["secant_root_over_km"], "D": D})
    print(out[-1], flush=True)
Path(sys.argv[1] if len(sys.argv) > 1 else "k4_stops_indep.json").write_text(json.dumps(out, indent=1) + "\n")
print("n", len(out), "min D", min(x["D"] for x in out))
