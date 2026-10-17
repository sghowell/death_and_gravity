#!/usr/bin/env python3
"""S4-6 review: point values inside the ACCEPTED boxes of logged covers.

For every accepted node of a cover tree (cover_replay.py) it draws points
(y, eps, tau) inside P(raw) (cover_verify.py's semantics) and evaluates
Sigma_c there with direct_chain.py (an independent chain, the midpoint
background of the certified balls, high precision).  Each value must lie in
the node's logged enclosure [lower, upper]; in particular upper < 0 must be
above it.  y is taken from a per-tree grid (log-spaced, plus the piece ends;
one build of the point pairs per grid value) restricted to the box's
y-projection, or the projection's midpoint when no grid value falls inside.
tau is drawn log-uniformly in the box's tau-range at that y (cut at 1e18 for
an unbounded range, and 1e18 itself is used once per such box); eps
uniformly, with the cell ends.  A value whose ball is not tight (radius
> 1e-6 relative) is recomputed at doubled precision, else counted undecided.

usage: python point_check.py trees/c29_t14.json.gz [...] [--per-box N] [--grid N]
       [--seed S] [--out FILE]
"""
import gzip
import json
import math
import os
import random
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
INF = float("inf")


def f(x):
    return INF if x == "inf" else (-INF if x == "-inf" else float(x))


def yproj(b):
    la, lb, e0, e1, s0, s1, t0, t1 = b
    ymin = max(la, (s0 / t1) if t1 < INF else 0.0)
    ymax = min(lb, (s1 / t0) if s1 < INF else INF)
    return ymin, ymax


def main(argv):
    import direct_chain as dc
    import p4.validated.tc_energy as te
    from p4.validated.arbseries import precision
    per_box, ngrid, seed, out_path = 2, 300, 20261017, None
    paths = []
    k = 0
    while k < len(argv):
        a = argv[k]
        if a == "--per-box":
            per_box = int(argv[k + 1]); k += 2
        elif a == "--grid":
            ngrid = int(argv[k + 1]); k += 2
        elif a == "--seed":
            seed = int(argv[k + 1]); k += 2
        elif a == "--out":
            out_path = argv[k + 1]; k += 2
        else:
            paths.append(a); k += 1
    rng = random.Random(seed)
    at = te.Atoms(full=True)
    exps = {}

    def expansion(prec):
        if prec not in exps:
            exps[prec] = dc.midpoint_expansion(prec)
        return exps[prec]

    pair_cache = {}

    def value(y, eps, tau):
        p = dc.needed_prec(y)
        for _ in range(3):
            key = (y, p)
            with precision(p):
                if key not in pair_cache:
                    pair_cache[key] = dc.pairs_at(at, expansion(max(640, p)), y)
                v, _ = dc.sigma_c(pair_cache[key], y, eps, tau)
            if v.is_finite() and float(v.rad()) <= 1e-6 * max(1.0, abs(float(v.mid()))):
                return v
            p *= 2
        return None

    report = dict(schema="review-s46-point-check/1", seed=seed, per_box=per_box, grid=ngrid,
                  trees=[])
    for path in paths:
        tstart = time.time()
        with gzip.open(path, "rt") as fh:
            d = json.load(fh)
        lo, hi = d["lo"], d["hi"]
        glo = lo if lo > 0 else 1e-14
        grid = sorted({math.exp(math.log(glo) + (math.log(hi) - math.log(glo)) * i / (ngrid - 1))
                       for i in range(ngrid)} | {hi} | ({lo} if lo > 0 else set()))
        acc = [n for n in d["nodes"] if n[4] == "accepted"]
        npts = nfail = nund = 0
        fails = []
        worst_ratio = None          # min of value.upper / box upper (both < 0; >= 1 is inside)
        closest = None
        sup_val = None
        for n in acc:
            raw = [f(x) for x in n[0]]
            lower, upper = n[2], f(n[3])
            ymin, ymax = yproj(raw)
            ymin = max(ymin, 1e-300)
            if ymax <= 0 or ymin > ymax:
                continue
            ys = [y for y in grid if ymin <= y <= ymax]
            ys = rng.sample(ys, min(per_box, len(ys))) if ys else \
                [math.sqrt(ymin * ymax) if ymin > 0 else ymax / 2]
            for j, y in enumerate(ys):
                la, lb, e0, e1, s0, s1, t0, t1 = raw
                tlo = max(t0, s0 / y)
                thi = min(t1, s1 / y)
                if tlo > thi:
                    continue
                cap = thi if thi < INF else 1e18
                if cap < tlo:
                    cap = tlo
                taus = [math.exp(rng.uniform(math.log(tlo), math.log(cap)))]
                if thi == INF and j == 0:
                    taus.append(1e18 if 1e18 > tlo else tlo * 1e6)
                for tau in taus:
                    eps = rng.choice([e0, e1, rng.uniform(e0, e1)])
                    v = value(y, eps, tau)
                    npts += 1
                    if v is None:
                        nund += 1
                        continue
                    vu, vl = float(v.upper()), float(v.lower())
                    sup_val = vu if sup_val is None else max(sup_val, vu)
                    bad = not (vu <= upper and (lower is None or vl >= lower))
                    r = vu / upper if upper < 0 else None
                    if r is not None and (worst_ratio is None or r < worst_ratio):
                        worst_ratio = r
                        closest = dict(y=y, eps=eps, tau=tau, value=[vl, vu],
                                       enclosure=[lower, upper], box=n[0])
                    if bad:
                        nfail += 1
                        if len(fails) < 20:
                            fails.append(dict(y=y, eps=eps, tau=tau, value=[vl, vu],
                                              enclosure=[lower, upper], box=n[0]))
        rec = dict(path=os.path.relpath(path, HERE), piece=d["piece"], tau_min=d["tau_min"],
                   accepted_boxes=len(acc), points=npts, failures=nfail, undecided=nund,
                   fail_rows=fails, sup_direct_value=sup_val,
                   min_value_over_box_bound=worst_ratio, closest_point=closest,
                   seconds=time.time() - tstart)
        report["trees"].append(rec)
        print("%s tau>=%g: %d boxes, %d points, %d failures, %d undecided, sup value %.6g, "
              "min value/bound %.9g, %.0fs" % (d["piece"], d["tau_min"], len(acc), npts, nfail,
                                               nund, sup_val, worst_ratio or float("nan"),
                                               time.time() - tstart), flush=True)
    report["ok"] = all(t["failures"] == 0 and t["undecided"] == 0 for t in report["trees"])
    if out_path:
        with open(out_path, "w") as fh:
            json.dump(report, fh, indent=1)


if __name__ == "__main__":
    main(sys.argv[1:])
