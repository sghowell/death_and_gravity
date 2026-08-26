"""Aggregate certified results (chains + feasible points) into results/summary.md, and optionally
certify feasible points for every finished chain at that chain's own T.

usage: PYTHONPATH=problems/P9/src uv run python -m p9.report_certified [--certify-feasible]
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

RESULTS = Path(__file__).resolve().parents[2] / "results"
CERT = RESULTS / "certificates"


def chains():
    out = []
    for d in sorted(CERT.glob("lkr_L*_D*_r*")):
        st = d / "state.json"
        if not st.exists():
            continue
        s = json.loads(st.read_text())
        done = (RESULTS / f"lkr_cert_{s['tag']}.log").exists() and "TOTAL done" in (RESULTS / f"lkr_cert_{s['tag']}.log").read_text() \
            or any("TOTAL done" in p.read_text() for p in RESULTS.glob(f"lkr_cert_{s['tag']}*.log"))
        out.append(dict(dir=d.name, tag=s["tag"], L=s["L"], Delta=s["Delta"], refine=s["refine"], T=s["T"],
                        H0_max=s["H0_max"], passes=len(s["history"]), done=done,
                        chi2_ref=s["reference"]["chi2_enclosure"]))
    return out


def feasible(tag_L, tag_D, refine):
    p = CERT / f"feasible_L{tag_L:g}_D{tag_D:g}_r{refine}.json"
    if not p.exists():
        return None
    f = json.loads(p.read_text())
    return f


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--certify-feasible", action="store_true"); a = ap.parse_args()
    rows = chains()
    if a.certify_feasible:
        for c in rows:
            if not c["done"]:
                continue
            f = feasible(c["L"], c["Delta"], c["refine"])
            if f is None or abs(f.get("T", -1) - c["T"]) > 1e-9:
                print(f"certifying feasible point for {c['tag']} at T={c['T']}", flush=True)
                subprocess.run([sys.executable, "-m", "p9.certify_feasible", "--L", str(c["L"]), "--Delta", str(c["Delta"]),
                                "--refine", str(c["refine"]), "--T", repr(c["T"])], check=False)
    lines = ["# P9(a) certified results", "",
             "Class C(G, L) on the midpoint-refined grid (refine r), M' ∈ [0, 40], frozen DESI DR2 BAO + Pantheon+ (C̃ as recorded),",
             "r_d ≥ 146.57 Mpc; T = upper(χ²(reference) + Δ). Upper bounds: chains of Arb-verified conic dual certificates;",
             "lower bounds: certified feasible class members. Values in km s⁻¹ Mpc⁻¹.", "",
             "| L | Δ | r | T | certified H₀ max (upper) | certified feasible H₀ (lower) | passes | status |", "|---|---|---|---|---|---|---|---|"]
    import math
    for c in sorted(rows, key=lambda r: (r["L"], r["Delta"])):
        f = feasible(c["L"], c["Delta"], c["refine"])
        lo = f["H0_lower_bound"] if (f and f.get("certified") and abs(f.get("T", -1) - c["T"]) < 1e-9) else None
        up = math.ceil(c["H0_max"] * 1e4) / 1e4
        lines.append(f"| {c['L']:g} | {c['Delta']:g} | {c['refine']} | {c['T']:.4f} | **≤ {up:.4f}** | "
                     f"{'≥ %.4f' % lo if lo is not None else '-'} | {c['passes']} | {'done' if c['done'] else 'running'} ({c['dir']}) |")
    (RESULTS / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
