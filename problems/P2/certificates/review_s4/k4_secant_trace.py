"""Trace the float secant of K4 (qnm_float.Solver.find, dps 20, maxsteps 12) at three of the 61 stops in H+, to see
why it stopped (S4 review, Amendment 6).  Uses the author's float solver unchanged; prints each iterate, |W|, arg W."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import mpmath as mp

W = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(W / "src"))
from p2.qnm_float import Geom, Solver  # noqa: E402

rows = {r["i"]: r for r in json.loads((W / "certificates/part1/controls_PB.json").read_text())["K4"]["rows"]}
out = []
g = Geom.frozen("PB")
for i in (0, 11, 23):
    r = rows[i]
    with mp.workdps(20):
        sol = Solver(g, K=6)
        km = g.km
        seed = mp.mpc(mp.mpf(r["w_over_km"][0]), mp.mpf(r["w_over_km"][1])) * km
        f = lambda x: sol.wronskian(x, r["l"])
        x0 = seed
        x1 = x0 * (1 + mp.mpf("1e-6")) + mp.mpc(0, "1e-9")
        f0, f1 = f(x0), f(x1)
        tr = [[mp.nstr(x0 / km, 15), mp.nstr(abs(f0), 10), mp.nstr(mp.arg(f0), 10)],
              [mp.nstr(x1 / km, 15), mp.nstr(abs(f1), 10), mp.nstr(mp.arg(f1), 10)]]
        tol = mp.mpf(10) ** (-mp.mp.dps + 12)
        stop = None
        for k in range(12):
            x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
            x0, f0, x1 = x1, f1, x2
            f1 = f(x1)
            rel = abs(x1 - x0) / abs(x1)
            tr.append([mp.nstr(x1 / km, 15), mp.nstr(abs(f1), 10), mp.nstr(mp.arg(f1), 10), mp.nstr(rel, 3)])
            if rel <= tol:
                stop = k
                break
        # finite-difference phase rate d arg W / d Re w at the seed, at two step sizes
        dw1, dw2 = mp.mpf("1e-6") * abs(seed), mp.mpf("1e-3") * km
        rate = [mp.nstr((mp.arg(f(seed + dw)) - mp.arg(f(seed))) / dw, 6) for dw in (dw1, dw2)]
    rec = {"i": i, "l": r["l"], "regime": r["regime"], "trace": tr, "stopped_at_iteration": stop, "phase_rate_fd": rate}
    print(json.dumps(rec), flush=True)
    out.append(rec)
Path(sys.argv[1] if len(sys.argv) > 1 else "k4_secant_trace.json").write_text(json.dumps(out, indent=1) + "\n")
