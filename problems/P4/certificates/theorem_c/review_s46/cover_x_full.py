#!/usr/bin/env python3
"""S4-6 review: independent exact x-cover of all 576 strip pieces (no
s43_parts code; extends review_s45/cover_independent.py to the centre).

  * tube pieces: the unique cache step whose float label and step equal the
    plan's (x, h), radius 0, interval [x_k - h_k, x_k] (Fractions); sonic
    pieces [-hi, -lo];
  * the 546 non-centre intervals: union and gaps in exact arithmetic;
  * centre pieces in t = e^x: c00.lo = 0 (x = -inf), consecutive pieces
    touch exactly (hi_i = lo_{i+1} as Fractions), so their union is
    [-inf, log hi(c29)];
  * the junction: log hi(c29) minus the lowest end of the non-centre union,
    in ball arithmetic at 384 bits (must be > 0: no gap at the centre end);
  * every centre_t14 part's cover y-range, tau range, eps cells and status
    against its plan piece.

usage (repository root): python problems/P4/results/theorem_c/review_s46/cover_x_full.py [out.json]
"""
import json
import os
import sys
from fractions import Fraction

from flint import arb, ctx, fmpz

HERE = os.path.dirname(os.path.abspath(__file__))
P4 = os.path.abspath(os.path.join(HERE, "..", "..", ".."))


def farb(q):
    return arb(fmpz(q.numerator)) / arb(fmpz(q.denominator))


def main(out):
    ctx.prec = 384
    plan = json.load(open(os.path.join(P4, "results/theorem_c/parts/plan.json")))
    caches = {}
    for name in ("full", "gap", "deep", "band"):
        d = json.load(open(os.path.join(P4, "data/cache/tube_%s.json" % name)))
        st = []
        for s in d["steps"]:
            m, e, rm, _ = (int(v) for v in s["x"])
            st.append(dict(x=Fraction(m) * Fraction(2) ** e, rad0=(rm == 0),
                           h=Fraction(s["h"]), hf=float(s["h"]),
                           xf=float(arb(fmpz(m)) * arb(2) ** e)))
        caches[name] = st
    iv, errs = {}, []
    for p in plan["pieces"]:
        if p["kind"] == "centre":
            continue
        if p["kind"] == "sonic":
            iv[p["id"]] = (-Fraction(p["hi"]), -Fraction(p["lo"]))
            continue
        st = caches[p["cache"]]
        hits = [k for k, s in enumerate(st) if s["xf"] == p["x"] and s["hf"] == p["h"]]
        if len(hits) != 1 or not st[hits[0]]["rad0"]:
            errs.append((p["id"], hits))
            continue
        s = st[hits[0]]
        iv[p["id"]] = (s["x"] - s["h"], s["x"])
    ints = sorted(iv.values())
    lo_all, cur, gaps = ints[0][0], ints[0][0], []
    for a, b in ints:
        if a > cur:
            gaps.append((float(cur), float(a)))
        cur = max(cur, b)
    cen = sorted((p for p in plan["pieces"] if p["kind"] == "centre"), key=lambda p: p["index"])
    touch = all(Fraction(a["hi"]) == Fraction(b["lo"]) for a, b in zip(cen, cen[1:]))
    junction = farb(Fraction(cen[-1]["hi"])).log() - farb(lo_all)
    # the centre parts against the plan
    bad = []
    for p in cen:
        d = json.load(open(os.path.join(P4, "results/theorem_c/parts/centre_t14", p["id"] + ".json")))
        par = d["result"]["params"]
        cells = par["eps_cells"]
        good = (d["status"] == "pass" and par["la"] == p["lo"] and par["lb"] == p["hi"]
                and par["tau_min"] == 14.0 and d["scope"]["tau"] == [14.0, "inf"]
                and cells[0][0] == 0.0 and cells[-1][1] == 6.0
                and all(a[1] == b[0] for a, b in zip(cells, cells[1:]))
                and d["certified_bound"] < 0
                and (par["s_lo"] <= p["lo"] * 14.0) and (p["lo"] > 0 or par["s_lo"] == 0.0))
        if not good:
            bad.append(p["id"])
    rec = dict(schema="review-s46-cover-x-full/1",
               noncentre_pieces=len(iv), noncentre_errors=errs,
               noncentre_union=[float(lo_all), float(cur)], noncentre_union_lower_exact=str(lo_all),
               noncentre_gaps=gaps,
               centre_pieces=len(cen), c00_lo=cen[0]["lo"], centre_touch_exactly=touch,
               centre_union_x=["-inf", float(farb(Fraction(cen[-1]["hi"])).log())],
               junction_overlap=junction.str(20, radius=True),
               junction_positive=bool(junction > 0),
               centre_parts_mismatch=bad,
               total_pieces=len(iv) + len(cen))
    rec["ok"] = (not errs and not gaps and float(cur) == 0.0 and touch and cen[0]["lo"] == 0.0
                 and rec["junction_positive"] and not bad and rec["total_pieces"] == 576)
    json.dump(rec, open(out, "w"), indent=1)
    print(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "cover_x_full.json"))
