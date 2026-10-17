#!/usr/bin/env python3
"""S4-6 review: exact check of a logged centre cover (cover_replay.py).

No p4 code is imported.  A box b = (la, lb, e0, e1, s0, s1, t0, t1) stands
for P(b) = {(y, eps, tau): y in [la, lb], y > 0, eps in [e0, e1],
tau in [t0, t1], tau finite, y tau in [s0, s1]} (tl.box_geometry's
semantics).  Every float is read as the exact dyadic rational it is
(Fraction); inf stays inf.  The check proves, for the piece [lo, hi] and
tau_min T,

    D = {y in [lo, hi], y > 0} x [0, 6] x [T, inf)  is contained in the union
        of P(raw) over the accepted nodes,

by induction over the logged tree:
  1. roots: all have y-range [lo, hi]; their eps cells tile [0, 6]; per eps
     cell the s-cells tile [s_first, inf) with s_first <= lo T; per s-cell
     the tau ranges tile [T, inf).  So D is in the union of the roots' P.
  2. every evaluated node: P(raw) is contained in P(nb) (exact projections
     of P(raw) onto y, tau and s inside nb's ranges, same eps), or, when the
     code called it empty, P(raw) is empty (exactly).
  3. every split: the two children are nb cut at one point m of one
     coordinate (shared endpoint, all else equal), so P(nb) = P(c1) u P(c2).
  4. every leaf is accepted (float upper bound < 0) or empty; nothing
     failed; the tree is consistent (every non-root node is the child of
     exactly one split; the roots are exactly the logged root boxes).
The enclosure of each accepted box is valid on tl.box_geometry's hulls of
its raw box, which contain P(nb) by construction (nb is built from them);
whether the enclosure itself is right is point_check.py's and the code
review's subject.

usage: python cover_verify.py trees/c00_t14.json.gz [...]  [--out FILE]
"""
import gzip
import json
import os
import sys
from fractions import Fraction

INF = float("inf")
HERE = os.path.dirname(os.path.abspath(__file__))


def q(x):
    if x == "inf" or x == INF:
        return INF
    if x == "-inf":
        return -INF
    return Fraction(x)


def div(a, b):
    """a / b for a >= 0, b > 0 (either may be inf); inf/inf does not occur."""
    if b == INF:
        if a == INF:
            raise ValueError("inf/inf")
        return Fraction(0)
    if a == INF:
        return INF
    return a / b


def mul(a, b):
    if a == 0 or b == 0:
        return Fraction(0)
    if a == INF or b == INF:
        return INF
    return a * b


def proj(b):
    """Exact (nonempty, (ymin, ymax), (tmin, tmax), (smin, smax)) of P(b)."""
    la, lb, e0, e1, s0, s1, t0, t1 = b
    if not (e0 <= e1 and t0 <= t1 and s0 <= s1 and t0 < INF and la <= lb):
        return False, None, None, None
    ymin = max(la, div(s0, t1))
    ymax = min(lb, div(s1, t0))
    if ymax <= 0 or ymin > ymax:
        return False, None, None, None
    # y in [ymin, ymax], y > 0, tau in [max(t0, s0/y), min(t1, s1/y)]
    tlo = max(t0, div(s0, ymax))
    thi = min(t1, div(s1, ymin) if ymin > 0 else INF)
    slo = max(s0, mul(ymin, t0))
    shi = min(s1, mul(ymax, t1))
    return True, (ymin, ymax), (tlo, thi), (slo, shi)


def contains(nb, raw):
    ok, y, t, s = proj(raw)
    if not ok:
        return True, "raw empty"
    L0, L1, e0, e1, S0, S1, T0, T1 = nb
    if not (e0 == raw[2] and e1 == raw[3]):
        return False, "eps range changed"
    good = L0 <= y[0] and L1 >= y[1] and T0 <= t[0] and T1 >= t[1] \
        and S0 <= s[0] and S1 >= s[1]
    return good, None if good else "projection outside nb: y %s t %s s %s nb %s" % (
        [float(v) for v in y], [float(v) for v in t], [float(v) for v in s],
        [float(v) for v in nb])


VAR = dict(l=(0, 1), e=(2, 3), s=(4, 5), t=(6, 7))


def verify(path, d=None):
    if d is None:
        with gzip.open(path, "rt") as fh:
            d = json.load(fh)
    lo, hi, T = Fraction(d["lo"]), Fraction(d["hi"]), Fraction(d["tau_min"])
    nodes = d["nodes"]
    raw = [tuple(q(x) for x in n[0]) for n in nodes]
    nbs = [tuple(q(x) for x in n[1]) if n[1] is not None else None for n in nodes]
    errs = []
    # ---- 1. roots
    roots = [tuple(q(x) for x in r) for r in d["roots"]]
    if any((r[0], r[1]) != (lo, hi) for r in roots):
        errs.append("a root's y-range is not [lo, hi]")
    cells = sorted({(r[2], r[3]) for r in roots})
    if not (cells[0][0] == 0 and cells[-1][1] == 6
            and all(a[1] == b[0] for a, b in zip(cells, cells[1:]))):
        errs.append("eps cells do not tile [0, 6]")
    for c in cells:
        rs = [r for r in roots if (r[2], r[3]) == c]
        scells = sorted({(r[4], r[5]) for r in rs})
        if not (scells[0][0] <= lo * T and scells[-1][1] == INF
                and all(a[1] == b[0] for a, b in zip(scells, scells[1:]))):
            errs.append("s-cells of eps cell %s do not tile [<= lo T, inf)" % (c,))
        for sc in scells:
            ts = sorted((r[6], r[7]) for r in rs if (r[4], r[5]) == sc)
            if not (ts[0][0] == T and ts[-1][1] == INF
                    and all(a[1] == b[0] for a, b in zip(ts, ts[1:]))):
                errs.append("tau ranges of cell %s, s %s do not tile [T, inf)" % (c, sc))
    # ---- 4. tree consistency
    parent = {}
    for i, n in enumerate(nodes):
        if n[4] == "split":
            for c in n[6]:
                if c in parent:
                    errs.append("node %d has two parents" % c)
                parent[c] = i
    root_nodes = [i for i in range(len(nodes)) if i not in parent]
    if sorted(raw[i] for i in root_nodes) != sorted(roots):
        errs.append("root nodes differ from the logged root boxes")
    outc = {}
    for n in nodes:
        outc[n[4]] = outc.get(n[4], 0) + 1
    if set(outc) - {"accepted", "split", "empty"}:
        errs.append("outcomes other than accepted/split/empty: %s" % outc)
    # ---- 2. normalization / emptiness, 3. splits, leaves
    n_cont = n_emp = n_split = n_acc = 0
    worst = None
    for i, n in enumerate(nodes):
        if n[4] == "empty":
            ok, *_ = proj(raw[i])
            if ok:
                errs.append("node %d called empty but P(raw) is not" % i)
            n_emp += 1
            continue
        if nbs[i] is None:
            errs.append("node %d has no normalized box" % i)
            continue
        good, why = contains(nbs[i], raw[i])
        if not good:
            errs.append("node %d: %s" % (i, why))
        n_cont += 1
        if n[4] == "accepted":
            up = q(n[3]) if n[3] not in ("inf", "-inf") else q(n[3])
            if not up < 0:
                errs.append("node %d accepted with upper %s" % (i, n[3]))
            worst = up if worst is None else max(worst, up)
            n_acc += 1
        elif n[4] == "split":
            var, kids = n[5], n[6]
            if len(kids) != 2:
                errs.append("node %d: %d children" % (i, len(kids)))
                continue
            a, b = raw[kids[0]], raw[kids[1]]
            i0, i1 = VAR[var]
            nb = nbs[i]
            same = all(a[k] == nb[k] and b[k] == nb[k] for k in range(8) if k not in (i0, i1))
            cut = a[i0] == nb[i0] and b[i1] == nb[i1] and a[i1] == b[i0] \
                and nb[i0] <= a[i1] <= nb[i1]
            if not (same and cut):
                errs.append("node %d: children are not a cut of nb along %s" % (i, var))
            n_split += 1
    res = d["result"]
    if n_acc != res["boxes"]:
        errs.append("accepted %d != boxes %d" % (n_acc, res["boxes"]))
    if worst is not None and float(worst) != res["worst"]:
        errs.append("worst %r != result worst %r" % (float(worst), res["worst"]))
    out = dict(path=os.path.relpath(path, HERE), piece=d["piece"], tau_min=d["tau_min"],
               lo=d["lo"], hi=d["hi"], roots=len(roots), eps_cells=len(cells),
               nodes=len(nodes), accepted=n_acc, split=n_split, empty=n_emp,
               containment_checked=n_cont, worst=(float(worst) if worst is not None else None),
               same_as_part=d["same_as_part"], replay_cpu_s=d["certify_cpu_s"],
               errors=errs[:20], n_errors=len(errs), ok=not errs and all(d["same_as_part"].values()))
    return out


def controls(path):
    """Corrupt a copy of the tree five ways; each must be reported."""
    import copy
    with gzip.open(path, "rt") as fh:
        d0 = json.load(fh)
    out = {}
    acc = [i for i, n in enumerate(d0["nodes"]) if n[4] == "accepted"]
    spl = [i for i, n in enumerate(d0["nodes"]) if n[4] == "split"]
    # 1. an accepted node's normalized box shrunk in y below the raw box's reach
    d = copy.deepcopy(d0)
    i = acc[len(acc) // 2]
    L0, L1 = q(d["nodes"][i][1][0]), q(d["nodes"][i][1][1])
    d["nodes"][i][1][1] = float(L0 + (L1 - L0) / 1000)
    out["nb_shrunk_below_raw"] = verify(path, d)["n_errors"] > 0
    # 2. a split whose second child starts above the cut (a sliver left out)
    d = copy.deepcopy(d0)
    j = spl[0]
    var, kids = d["nodes"][j][5], d["nodes"][j][6]
    i0, i1 = VAR[var]
    b = d["nodes"][kids[1]][0]
    lo_, hi_ = q(b[i0]), q(b[i1])
    b[i0] = float(lo_ + (hi_ - lo_) / 1000) if hi_ != INF else float(lo_ + max(abs(lo_) / 1000, Fraction(1, 10 ** 6)))
    out["split_gap"] = verify(path, d)["n_errors"] > 0
    # 3. an accepted node with upper bound >= 0
    d = copy.deepcopy(d0)
    d["nodes"][acc[0]][3] = 1e-9
    out["accepted_nonnegative"] = verify(path, d)["n_errors"] > 0
    # 4. a nonempty node relabelled empty
    d = copy.deepcopy(d0)
    d["nodes"][acc[1]][4] = "empty"
    out["false_empty"] = verify(path, d)["n_errors"] > 0
    # 5. a root eps cell narrowed (the roots no longer tile [0, 6])
    d = copy.deepcopy(d0)
    d["roots"] = [r for r in d["roots"] if not (r[2] == 0.0)]
    out["root_cell_missing"] = verify(path, d)["n_errors"] > 0
    return out


def main(argv):
    out_path = None
    if "--out" in argv:
        k = argv.index("--out")
        out_path = argv[k + 1]
        argv = argv[:k] + argv[k + 2:]
    ctl = None
    if "--controls" in argv:
        k = argv.index("--controls")
        ctl = controls(argv[k + 1])
        argv = argv[:k] + argv[k + 2:]
        print("controls (True = the corruption was reported):", ctl)
    res = [verify(p) for p in argv]
    for r in res:
        print("%s tau>=%g: ok=%s nodes=%d accepted=%d split=%d empty=%d errors=%d same=%s"
              % (r["piece"], r["tau_min"], r["ok"], r["nodes"], r["accepted"], r["split"],
                 r["empty"], r["n_errors"], r["same_as_part"]))
        for e in r["errors"][:5]:
            print("   ", e)
    if out_path:
        with open(out_path, "w") as fh:
            json.dump(dict(schema="review-s46-cover-verify/1", trees=res, controls=ctl,
                           ok=all(r["ok"] for r in res)
                           and (ctl is None or all(ctl.values()))), fh, indent=1)


if __name__ == "__main__":
    main(sys.argv[1:])
