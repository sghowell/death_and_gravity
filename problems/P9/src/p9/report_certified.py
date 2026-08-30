"""Aggregate certified results (chains + feasible points) into results/summary.md, and optionally
certify feasible points for every finished chain at that chain's own T.

usage: PYTHONPATH=problems/P9/src uv run python -m p9.report_certified [--certify-feasible]
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path

from .variants import RD_BOX_PLANCK, Variant

RESULTS = Path(__file__).resolve().parents[2] / "results"
CERT = RESULTS / "certificates"
R_D_REF = 147.09     # Mpc: Planck 2018 central value, the reference point of the rescaled column


def chain_log_done(tag: str) -> bool:
    """A chain is finished when a run log that opens with its tag ends with the 'TOTAL done' marker
    (log names vary between the curve and variant queues, so match on the header line, not the file name)."""
    for p in RESULTS.glob("lkr_cert_*.log"):
        txt = p.read_text()
        head = txt.split("\n", 1)[0]
        if (head.startswith(tag + " ") or head.startswith(tag + ":") or head.startswith(f"resuming {tag}:")) and "TOTAL done" in txt:
            return True
    return False


def chains():
    out = []
    for d in sorted(CERT.glob("lkr_L*_D*_r*")):
        st = d / "state.json"
        if not st.exists():
            continue
        s = json.loads(st.read_text())
        done = bool(s.get("done")) or chain_log_done(s["tag"])
        var = Variant.from_state(s)
        out.append(dict(dir=d.name, tag=s["tag"], L=s["L"], Delta=s["Delta"], refine=s["refine"], T=s["T"],
                        H0_max=s["H0_max"], passes=len(s["history"]), done=done, var=var,
                        chi2_ref=s["reference"]["chi2_enclosure"]))
    return out


def feasible(tag_L, tag_D, refine, var: Variant):
    p = CERT / f"feasible_L{tag_L:g}_D{tag_D:g}_r{refine}{var.tag}.json"
    return json.loads(p.read_text()) if p.exists() else None


def variant_flags(var: Variant) -> list[str]:
    fl = ["--sn", var.sn] + (["--dv"] if var.use_dv else [])
    if (var.r_lo, var.r_hi) != RD_BOX_PLANCK:
        fl += ["--rd_box", repr(var.r_lo), repr(var.r_hi)]
    return fl


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--certify-feasible", action="store_true"); a = ap.parse_args()
    rows = chains()
    if a.certify_feasible:
        for c in rows:
            if not c["done"]:
                continue
            f = feasible(c["L"], c["Delta"], c["refine"], c["var"])
            if f is None or abs(f.get("T", -1) - c["T"]) > 1e-9:
                print(f"certifying feasible point for {c['tag']} at T={c['T']}", flush=True)
                subprocess.run([sys.executable, "-m", "p9.certify_feasible", "--L", str(c["L"]), "--Delta", str(c["Delta"]),
                                "--refine", str(c["refine"]), "--T", repr(c["T"])] + variant_flags(c["var"]), check=False)
    lines = ["# P9(a) certified results", "",
             "Class C(G, L) on the midpoint-refined grid (refine r), M' ∈ [0, 40], frozen DESI DR2 BAO + the SN sample of the row",
             "(C̃ as recorded); T = upper(χ²(reference) + Δ). Upper bounds: chains of Arb-verified conic dual certificates;",
             "lower bounds: certified feasible class members. Values in km s⁻¹ Mpc⁻¹. The certified bound is on ũ₀ = c/(r_d H₀);",
             f"r_d enters only through H₀ = c/(r_d ũ₀), so for every r_d in the row's box H₀ ≤ H₀max · r_lo/r_d exactly, and the",
             f"'at r_d = {R_D_REF}' column is H₀max · r_lo/{R_D_REF} (rescaling to any other r_d is the same one-line identity).",
             "Baseline rows: Pantheon+, BGS D_V row dropped, r_d box = Planck 2018 ±2σ (FORMULATION §1); variants per §6.1.", "",
             f"| L | Δ | r | SN | D_V | r_d box [Mpc] | T | certified H₀ max (upper, at r_lo) | at r_d = {R_D_REF} | "
             "certified feasible H₀ (lower) | passes | status |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for c in sorted(rows, key=lambda r: (r["var"].sn, r["var"].use_dv, r["var"].r_lo, r["L"], r["Delta"])):
        v = c["var"]
        f = feasible(c["L"], c["Delta"], c["refine"], v)
        lo = f["H0_lower_bound"] if (f and f.get("certified") and abs(f.get("T", -1) - c["T"]) < 1e-9) else None
        # upper bounds rounded UP, lower bounds rounded DOWN (never round-to-nearest a certified bound)
        up = math.ceil(c["H0_max"] * 1e4) / 1e4
        up_ref = math.ceil(c["H0_max"] * v.r_lo / R_D_REF * 1e4) / 1e4
        lo_txt = "≥ %.4f" % (math.floor(lo * 1e4) / 1e4) if lo is not None else "-"
        lines.append(f"| {c['L']:g} | {c['Delta']:g} | {c['refine']} | {v.sn} | {'yes' if v.use_dv else 'no'} | "
                     f"[{v.r_lo:g}, {v.r_hi:g}] | {c['T']:.4f} | **≤ {up:.4f}** | ≤ {up_ref:.4f} | "
                     f"{lo_txt} | {c['passes']} | {'done' if c['done'] else 'running'} ({c['dir']}) |")
    (RESULTS / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
