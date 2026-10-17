"""S4-5 review: independent exact x-cover of the 546 non-centre pieces, from
plan.json and the four tube cache files only (no s43_parts code).

Tube-kind piece: the cache step whose float label float(x_k) and step h_k equal
the plan's (x, h) -- searched over ALL steps, required unique -- with radius 0;
its interval is [x_k - h_k, x_k] exactly.  Sonic piece: [-hi, -lo].  Also
checks every layer_t14 part's cover parameters (la, lb, laurent) against the
piece.  Prints: exactness facts, the union's gaps on [-4.56, 0], the smallest
cross-region overlaps, and the overlap the centre must meet.
"""
import json
import math
import sys
from fractions import Fraction

from flint import arb, fmpz, ctx

ctx.prec = 384
P4 = sys.argv[1]
plan = json.load(open(P4 + "/results/theorem_c/parts/plan.json"))
caches = {}
for name in ("full", "gap", "deep", "band"):
    d = json.load(open(P4 + "/data/cache/tube_%s.json" % name))
    st = []
    for s in d["steps"]:
        m, e, rm, re_ = (int(v) for v in s["x"])
        x = Fraction(m) * Fraction(2) ** e
        st.append(dict(x=x, rad0=(rm == 0), h=Fraction(s["h"]), hf=float(s["h"]),
                       xf=float(arb(fmpz(m)) * arb(2) ** e)))
    caches[name] = st
facts = {}
for name, st in caches.items():
    facts[name] = dict(steps=len(st), all_radius0=all(s["rad0"] for s in st),
                       joins_exact=all(st[k + 1]["x"] == st[k]["x"] - st[k]["h"] for k in range(len(st) - 1)),
                       top=float(st[0]["x"]), bottom=float(st[-1]["x"] - st[-1]["h"]))
print("caches:", json.dumps(facts))
iv, errs = {}, []
for p in plan["pieces"]:
    if p["kind"] == "centre":
        continue
    if p["kind"] == "sonic":
        iv[p["id"]] = (-Fraction(p["hi"]), -Fraction(p["lo"]), p["region"])
        continue
    st = caches[p["cache"]]
    hits = [k for k, s in enumerate(st) if s["xf"] == p["x"] and s["hf"] == p["h"]]
    if len(hits) != 1:
        errs.append((p["id"], "hits", hits))
        continue
    s = st[hits[0]]
    if not s["rad0"]:
        errs.append((p["id"], "radius"))
    iv[p["id"]] = (s["x"] - s["h"], s["x"], p["region"])
print("pieces:", len(iv), "errors:", errs)
# union and gaps on [-4.56..., 0]
ints = sorted((a, b, r, i) for i, (a, b, r) in iv.items())
lo_all = ints[0][0]
cur, gaps = lo_all, []
for a, b, r, i in ints:
    if a > cur:
        gaps.append((float(cur), float(a)))
    cur = max(cur, b)
print("union lower end:", float(lo_all), " exact:", lo_all, " upper end:", float(cur), " gaps:", gaps)
# smallest overlap between pieces of different regions that touch / overlap
best = []
for j, (a, b, r, i) in enumerate(ints):
    for a2, b2, r2, i2 in ints[j + 1:j + 400]:
        if r2 == r:
            continue
        ov = min(b, b2) - max(a, a2)
        if ov >= 0:
            best.append((float(ov), i, i2))
best.sort()
print("smallest cross-region overlaps:", best[:6])
# the centre must reach above the deep end
c29 = [p for p in plan["pieces"] if p["id"] == "centre/c29"][0]
print("deep end x =", float(lo_all), "; centre c29 hi in t:", c29["hi"], "log ->",
      arb(Fraction(c29["hi"]).numerator).__truediv__(arb(Fraction(c29["hi"]).denominator)).log().str(20),
      "; overlap with deep end:", float(arb(c29["hi"]).log() - arb(lo_all.numerator) / arb(lo_all.denominator)))
# layer part cover parameters and the step the region_call names
import re
bad = []
stepidx = {}
for p in plan["pieces"]:
    if p["kind"] == "tube":
        st = caches[p["cache"]]
        k_abs = [k for k, s in enumerate(st) if s["xf"] == p["x"] and s["hf"] == p["h"]][0]
        sel = [k for k, s in enumerate(st) if s["xf"] - s["hf"] >= p["x_min"] - 1e-12 and s["xf"] <= p["x_max"] + 1e-12]
        stepidx[p["id"]] = sel.index(k_abs)       # region_call prints the drivers' selection index
for p in plan["pieces"]:
    if p["kind"] == "centre":
        continue
    d = json.load(open(P4 + "/results/theorem_c/parts/layer_t14/%s.json" % p["id"]))
    par = d["result"]["params"]
    la, lb = par.get("la"), par.get("lb")
    rc = d.get("region_call", "")
    if p["kind"] == "tube":
        m = re.search(r"tube_(\w+)\.steps\[(\d+)\]", rc)
        okp = (la == 0 and lb == p["h"] and m is not None and m.group(1) == p["cache"]
               and int(m.group(2)) == stepidx[p["id"]])
    else:
        okp = (la == p["lo"] and lb == p["hi"])
    okp = okp and d["status"] == "pass" and d["scope"]["tau"][0] == 14.0 and d["scope"]["tau"][1] == "inf" \
        and d["scope"]["eps"] == [0.0, 6.0] and d["certified_bound"] < 0 and par["tau_min"] == 14.0 \
        and [c[0] for c in par["eps_cells"]][0] == 0.0 and par["eps_cells"][-1][1] == 6.0 \
        and all(a[1] == b[0] for a, b in zip(par["eps_cells"], par["eps_cells"][1:]))
    if not okp:
        bad.append((p["id"], la, lb, rc, d["status"]))
print("layer_t14 parts with a cover/scope/step mismatch:", len(bad), bad[:3])
