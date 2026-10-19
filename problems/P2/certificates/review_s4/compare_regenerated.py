"""Follow-up (2026-10-19): the certificates regenerated after the review against the versions it reviewed.

The reviewed versions are read from git (commit 0e329217); the current ones from certificates/part1/.
  * S (both sides, all 92 leaves): every field must agree, except printed intervals '[lo, hi]', which may
    change only by the switch from round-to-nearest to outward rounding (review W4): the new interval must
    contain the old printed endpoints and move each by at most one unit in the 8th significant digit.
  * L: leaves_detail must be identical; the summary extremes are compared with the ball endpoints.
  * controls, cover, replay: every changed leaf path is listed.
Usage: compare_regenerated.py OUT.json   (run inside the git checkout)
"""
from __future__ import annotations

import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

P2 = Path(__file__).resolve().parents[2]
CERT = P2 / "certificates" / "part1"
REVIEWED = "0e329217"


def old(name):
    rel = (CERT / name).relative_to(subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True,
                                                   cwd=P2, check=True).stdout.strip())
    return json.loads(subprocess.run(["git", "show", f"{REVIEWED}:{rel}"], capture_output=True, text=True, cwd=P2,
                                     check=True).stdout)


def iv(s):
    s = s.strip()
    if s.startswith("[") and "," in s and "+/-" not in s:
        a, b = s.strip("[]").split(",")
        return Fraction(a.strip()), Fraction(b.strip())
    return None


def ulp8(x: Fraction) -> Fraction:
    a = abs(x)
    if a == 0:
        return Fraction(0)
    k = 0
    while Fraction(10) ** k > a:
        k -= 1
    while Fraction(10) ** (k + 1) <= a:
        k += 1
    return Fraction(10) ** (k - 7)


def walk(o, n, notes, path=""):
    if isinstance(o, dict) and isinstance(n, dict):
        for k in sorted(set(o) | set(n)):
            if k not in o or k not in n:
                notes.append({"path": path + "." + k, "kind": "missing in " + ("old" if k not in o else "new")})
            else:
                walk(o[k], n[k], notes, path + "." + k)
        return
    if isinstance(o, list) and isinstance(n, list):
        if len(o) != len(n):
            notes.append({"path": path, "kind": "length", "old": len(o), "new": len(n)})
            return
        for i, (a, b) in enumerate(zip(o, n)):
            walk(a, b, notes, f"{path}[{i}]")
        return
    if o == n:
        return
    if isinstance(o, str) and isinstance(n, str) and iv(o) and iv(n):
        (lo_o, hi_o), (lo_n, hi_n) = iv(o), iv(n)
        ok = lo_n <= lo_o and hi_n >= hi_o and lo_o - lo_n <= ulp8(lo_o) and hi_n - hi_o <= ulp8(hi_o)
        d = max(lo_o - lo_n, hi_n - hi_o)
        rel = max((lo_o - lo_n) / abs(lo_o) if lo_o else 0, (hi_n - hi_o) / abs(hi_o) if hi_o else 0)
        notes.append({"path": path, "kind": "outward re-rounding" if ok else "VALUE CHANGE", "old": o, "new": n,
                      "abs": float(d), "rel": float(rel)})
        return
    notes.append({"path": path, "kind": "VALUE CHANGE", "old": o, "new": n})


def main():
    out = {}
    for tag in ("plus", "c"):
        o, n = old(f"S_PB_{tag}.json"), json.loads((CERT / f"S_PB_{tag}.json").read_text())
        notes = []
        for a, b in zip(o["leaves"], n["leaves"]):
            walk(a, b, notes, f"leaf{a['zeta_over_OmegaP']}")
        re_ = [x for x in notes if x["kind"] == "outward re-rounding"]
        other = [x for x in notes if x["kind"] != "outward re-rounding"]
        out[f"S_{tag}"] = {"leaves": [len(o["leaves"]), len(n["leaves"])], "outward_rerounded_strings": len(re_),
                           "max_abs_move": max((x["abs"] for x in re_), default=0),
                           "max_rel_move": max((x["rel"] for x in re_), default=0),
                           "other_changes": other,
                           "top_level_changes": {k: [o.get(k), n.get(k)] for k in set(o) | set(n)
                                                 if k != "leaves" and o.get(k) != n.get(k)}}
    for tag in ("plus", "c"):
        o, n = old(f"L_PB_{tag}.json"), json.loads((CERT / f"L_PB_{tag}.json").read_text())
        out[f"L_{tag}"] = {"leaves_detail_identical": o["leaves_detail"] == n["leaves_detail"],
                           "top_level_changes": {k: [o.get(k), n.get(k)] for k in set(o) | set(n)
                                                 if k != "leaves_detail" and o.get(k) != n.get(k)}}
    for name in ("controls_PB.json", "cover_PB.json", "replay_PB.json", "T_PB.json", "A_PB.json", "bulk_PB.json"):
        notes = []
        walk(old(name), json.loads((CERT / name).read_text()), notes)
        out[name] = notes
    print(json.dumps({k: (v if not isinstance(v, list) else f"{len(v)} changed paths") for k, v in out.items()}, indent=1)[:4000])
    Path(sys.argv[1] if len(sys.argv) > 1 else "compare_regenerated.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
