#!/usr/bin/env python3
"""S4-6 review: adversarial points against the logged covers.

Points at the extremes the centre certificate must handle: y -> 0 on the
closing piece (y down to 1e-40), tau -> infinity (up to 1e40), eps at 0, 6
and a cell edge, y at the piece ends, s = y tau at the s-grid edges
(s_lo, s_a = 1e-3, 0.5 = S_SPLIT, 8 = Z_LARGE, s_big = 1e3) and tau at the
pre-split edges (14, t_pre = 2e4).  For each point every ACCEPTED box whose
P(raw) contains it (exact test, as cover_verify.py) is found; there must be
at least one, and the independent direct value (direct_chain.py, midpoint
background, precision raised until tight) must lie in each such box's
enclosure.

usage: PYTHONPATH=problems/P4/src python adversarial_points.py [out.json]
"""
import gzip
import json
import math
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
INF = float("inf")


def q(x):
    return INF if x == "inf" else Fraction(x)


def inside(raw, y, eps, tau):
    la, lb, e0, e1, s0, s1, t0, t1 = (q(v) for v in raw)
    y, eps, tau = Fraction(y), Fraction(eps), Fraction(tau)
    s = y * tau
    return (la <= y <= lb and y > 0 and e0 <= eps <= e1 and t0 <= tau <= t1
            and s0 <= s <= s1)


def main(out):
    import direct_chain as dc
    import p4.validated.tc_energy as te
    from p4.validated.arbseries import precision
    at = te.Atoms(full=True)
    exps, pcache = {}, {}

    def value(y, eps, tau):
        p = dc.needed_prec(y)
        for _ in range(4):
            if p not in exps:
                exps[p] = dc.midpoint_expansion(max(640, p))
            with precision(p):
                if (y, p) not in pcache:
                    pcache[(y, p)] = dc.pairs_at(at, exps[p], y)
                v, _ = dc.sigma_c(pcache[(y, p)], y, eps, tau)
            if v.is_finite() and float(v.rad()) <= 1e-6 * max(1.0, abs(float(v.mid()))):
                return v, p
            p *= 2
        return None, p

    trees = {}
    for name in ("c00_t14", "c00_t224", "c29_t14", "c29_t224"):
        with gzip.open(os.path.join(HERE, "trees", name + ".json.gz"), "rt") as fh:
            trees[name] = json.load(fh)
    pts = []
    c00 = trees["c00_t14"]
    hi0 = c00["hi"]
    for y in (1e-40, 1e-30, 1e-20, 1e-13, hi0 * (1 - 1e-12), hi0):
        for eps in (0.0, 0.03, 6.0):
            for tau in (14.0, 14.0 * (1 + 1e-12), 224.0, 2e4, 1e8, 1e15, 1e40):
                pts.append(("c00", y, eps, tau))
    c29 = trees["c29_t14"]
    lo9, hi9 = c29["lo"], c29["hi"]
    for y in (lo9, lo9 * (1 + 1e-9), math.sqrt(lo9 * hi9), hi9 * (1 - 1e-12), hi9):
        for eps in (0.0, 3.4, 6.0):
            taus = {14.0, 224.0, 1e12, 1e30}
            for s in (lo9 * 14.0, 1e-3, 0.5, 8.0, 1e3):
                if s / y >= 14.0:
                    taus.add(s / y)
            for tau in sorted(taus):
                pts.append(("c29", y, eps, tau))
    rows, nfail, nnone, nund = [], 0, 0, 0
    for piece, y, eps, tau in pts:
        v, p = value(y, eps, tau)
        for name, tr in trees.items():
            if not name.startswith(piece):
                continue
            if tau < tr["tau_min"]:
                continue
            boxes = [n for n in tr["nodes"] if n[4] == "accepted" and inside(n[0], y, eps, tau)]
            row = dict(tree=name, y=y, eps=eps, tau=tau, s=y * tau, boxes=len(boxes),
                       prec=p, value=None if v is None else [float(v.lower()), float(v.upper())])
            if not boxes:
                nnone += 1
                row["fail"] = "no accepted box contains the point"
            elif v is None:
                nund += 1
                row["fail"] = "value not tight"
            else:
                bad = [n[0] for n in boxes
                       if not (float(v.upper()) <= q(n[3]) and float(v.lower()) >= n[2])]
                row["enclosures"] = [[n[2], n[3]] for n in boxes]
                if bad:
                    nfail += 1
                    row["fail"] = "value outside %d enclosures" % len(bad)
            rows.append(row)
    rec = dict(schema="review-s46-adversarial-points/1", points=len(pts), checks=len(rows),
               failures=nfail, uncovered=nnone, undecided=nund,
               max_boxes_per_point=max(r["boxes"] for r in rows),
               sup_value=max(r["value"][1] for r in rows if r["value"]),
               rows=rows, ok=(nfail == 0 and nnone == 0 and nund == 0))
    json.dump(rec, open(out, "w"), indent=1)
    print("points %d, checks %d, failures %d, uncovered %d, undecided %d, max boxes/point %d, "
          "sup value %.6g, ok=%s" % (len(pts), len(rows), nfail, nnone, nund,
                                     rec["max_boxes_per_point"], rec["sup_value"], rec["ok"]))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "adversarial_points.json"))
