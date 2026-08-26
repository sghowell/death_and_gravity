"""Render results/certified_curve.json as a markdown summary (results/summary.md)."""

from __future__ import annotations

import json
from pathlib import Path

RESULTS = Path(__file__).resolve().parents[2] / "results"


def main(name: str = "certified_curve.json"):
    d = json.loads((RESULTS / name).read_text())
    rows = sorted(d["rows"], key=lambda r: (r["L"], r["Delta"]))
    lines = ["# P9(a) — certified upper bounds on H0 (frozen DESI DR2 BAO + Pantheon+ + Planck r_d)", "",
             f"r_d box lower edge r_lo = {d['r_lo']} Mpc; bound-tightening subset = {d['subset']} SNe.",
             f"Flat LCDM (BAO+SN): Omega_m = {d['lcdm']['bao_sn']['om']:.4f}, h r_d = {d['lcdm']['bao_sn']['h_rd']:.2f} Mpc, "
             f"chi2 = {d['lcdm']['bao_sn']['chi2']:.2f}.", "",
             "| L | Delta | class-min chi2 | LCDM chi2 (if in class) | T | certified H0_max | feasible H0 (lower bound on max) | gap |",
             "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        feas = r["H0_feasible"]
        gap = (r["H0_max_certified"] - feas) if feas is not None else None
        lines.append(f"| {r['L']:g} | {r['Delta']:g} | {r['chi2_class_min']:.2f} | "
                     f"{'-' if r['chi2_lcdm'] is None else f'{r['chi2_lcdm']:.2f}'} | {r['T']:.2f} | "
                     f"**{r['H0_max_certified']:.3f}** | {'-' if feas is None else f'{feas:.3f}'} | "
                     f"{'-' if gap is None else f'{gap:.3f}'} |")
    lines += ["", "Certified: each row's H0_max is backed by a chain of conic dual certificates verified in ball "
              "arithmetic (see FORMULATION.md §4 and results/certificates/<tag>/). The feasible column is a "
              "verified class member with chi2 <= T, so the true maximum lies in [feasible, certified]."]
    (RESULTS / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
