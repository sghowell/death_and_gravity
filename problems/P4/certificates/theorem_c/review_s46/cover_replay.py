#!/usr/bin/env python3
"""S4-6 review: replay a centre part's cover with the whole box tree logged.

certify_logged() is tc_layer.certify_piece (tc_layer.py sha256 29a864dd...,
lines 611-729) copied statement for statement, with logging added and nothing
else changed: every box evaluated by _upper becomes a node (raw box, the
normalized box or None if empty, the enclosure of Sigma_c, its float upper
bound), and every popped node gets its outcome (accepted, split with the
chosen variable and the children, failed).  The run uses the parts' own
inputs: tc_layer_centre.centre_layer_data(at, ce, hi, raised=True), sigma_fn =
sigma_box_centre, target 0, max_evals 150000, precision 320, tau_min 14 or
224 (scripts/s46_centre.py work).  The counts (worst, boxes, evals, empty)
are compared with the committed part: equality shows the copy is faithful.

cover_verify.py then checks the logged tree in exact rational arithmetic,
and point_check.py samples points inside the accepted boxes.

usage (repository root):
    PYTHONPATH=problems/P4/src python problems/P4/results/theorem_c/review_s46/cover_replay.py c00 14 [c00:224 ...]
Output: review_s46/trees/<piece>_t<tau>.json.gz
"""
import gzip
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
PARTS = os.path.join(ROOT, "problems", "P4", "results", "theorem_c", "parts")
INF = float("inf")


def certify_logged(tl, ld, la, lb, eps_cells=None, s_lo=None, s_a=1e-3, s_big=None,
                   s_ratio=2.0, target=0.0, max_evals=400000, max_depth=60, log=None,
                   t_pre=2e4, t_ratio=1.5, tau_min=None, sigma_fn=None):
    nodes = []                       # [raw, nb, lower, upper(float), outcome, var, children]
    enc = {}

    def upper_logged(box):
        # tl._upper with the enclosure kept (same calls, same order)
        nb, g = tl._normalize(ld, box)
        if nb is None:
            return None, -INF, None
        try:
            S = (tl.sigma_box if sigma_fn is None else sigma_fn)(ld, nb, g)
        except (ArithmeticError, ZeroDivisionError, ValueError):
            return nb, INF, None
        if S is None:
            return None, -INF, None
        u = S.upper()
        return nb, (float(u) if u.is_finite() else INF), S

    def new_node(raw, nb, up, S):
        nodes.append([list(raw), (list(nb) if nb is not None else None),
                      (float(S.lower()) if S is not None else None), up, None, None, None])
        return len(nodes) - 1

    t_start = time.time()
    tmin = tl.TSTAR if tau_min is None else float(tau_min)
    if not tmin > 0:
        raise ValueError("tau_min must be positive")
    if not target <= 0:
        raise ValueError("target must be <= 0")
    eps_cells = tl.eps_cells_default() if eps_cells is None else eps_cells
    if not (eps_cells[0][0] == 0.0 and eps_cells[-1][1] == tl.EPS_MAX and all(
            a[1] == b[0] for a, b in zip(eps_cells, eps_cells[1:]))):
        raise ValueError("eps cells must tile [0, EPS_MAX]")
    xa = ld.xi_of_ell(la)
    if s_lo is None:
        s_lo = xa * tmin * (1 - 1e-12)
    if s_big is None:
        s_big = max(1e3, 1e4 * ld.xi_of_ell(lb))
    grid = []
    if s_lo <= 0:
        grid.append((0.0, s_a))
        s = s_a
    else:
        s = s_lo
    while s < s_big:
        s2 = min(s * s_ratio, s_big)
        grid.append((s, s2))
        s = s2
    grid.append((s_big, INF))
    edges = [tmin]
    while edges[-1] < t_pre:
        edges.append(min(edges[-1] * t_ratio, t_pre))
    slabs = list(zip(edges, edges[1:])) + [(edges[-1], INF)]
    xb = ld.xi_of_ell(lb)
    stack = []
    roots = []
    for e0, e1 in eps_cells:
        for s0, s1 in grid:
            if xa <= 0 and s0 < t_pre * xb:
                for ta, tb in slabs:
                    stack.append(((la, lb, e0, e1, s0, s1, ta, tb), 0, None, None))
            else:
                stack.append(((la, lb, e0, e1, s0, s1, tmin, INF), 0, None, None))
    roots = [list(b[0]) for b in stack]
    nev = ncert = nempty = 0
    worst = -INF
    fails = []
    hist = {}
    while stack:
        box, depth, up, nid = stack.pop()
        if up is None:
            nb, up, S = upper_logged(box)
            nev += 1
            nid = new_node(box, nb, up, S)
            if nb is None:
                nempty += 1
                nodes[nid][4] = "empty"
                continue
            box = nb
        if up < target or (target == 0 and up < 0):
            ncert += 1
            worst = max(worst, up)
            hist[depth] = hist.get(depth, 0) + 1
            nodes[nid][4] = "accepted"
            continue
        if depth >= max_depth or nev >= max_evals:
            fails.append(dict(box=list(box), upper=up, depth=depth))
            nodes[nid][4] = "failed"
            if len(fails) >= 8:
                break
            continue
        best = None
        for var in ("s", "e", "t", "l"):
            ch = tl._children(ld, box, var)
            if ch is None:
                continue
            res = []
            for cb in ch:
                nb, u, S = upper_logged(cb)
                nev += 1
                res.append((nb, u, cb, S))
            score = max(r[1] for r in res)
            if best is None or score < best[0]:
                best = (score, var, res)
        if best is None:
            fails.append(dict(box=list(box), upper=up, depth=depth, why="unsplittable"))
            nodes[nid][4] = "failed"
            if len(fails) >= 8:
                break
            continue
        kids = []
        for nb, u, cb, S in best[2]:
            cid = new_node(cb, nb, u, S)
            kids.append(cid)
            if nb is None:
                nempty += 1
                nodes[cid][4] = "empty"
                continue
            stack.append((nb, depth + 1, u, cid))
        nodes[nid][4], nodes[nid][5], nodes[nid][6] = "split", best[1], kids
    ok = not fails and not stack
    out = dict(ok=ok, worst=(worst if worst > -INF else None), boxes=ncert,
               evals=nev, empty=nempty, fails=fails[:8],
               depth_hist={str(k): v for k, v in sorted(hist.items())},
               time=time.time() - t_start,
               params=dict(s_a=s_a, s_big=s_big, s_ratio=s_ratio, zsplit=tl.ZSPLIT,
                           tstar=tl.TSTAR, tau_min=tmin, eps_cells=eps_cells, la=la, lb=lb,
                           target=target, s_lo=s_lo, t_pre=t_pre, t_ratio=t_ratio,
                           max_evals=max_evals,
                           evaluator=("sigma_box" if sigma_fn is None
                                      else getattr(sigma_fn, "__name__", str(sigma_fn))),
                           combine=bool(getattr(ld, "combine", False))))
    return out, roots, nodes


def jnum(x):
    if x is None:
        return None
    if x == INF:
        return "inf"
    if x == -INF:
        return "-inf"
    return x


def main(specs):
    import p4.validated.tc_energy as te
    import p4.validated.tc_layer as tl
    import p4.validated.tc_layer_centre as tcc
    from p4.validated.arbseries import precision
    from p4.validated.ec_centre import certified_centre_expansion
    at = te.Atoms(full=True)
    ce = certified_centre_expansion(256, K=30)
    assert tcc.check_nb_identity(at)["ok"]
    plan = json.load(open(os.path.join(PARTS, "plan.json")))
    byid = {p["id"].split("/")[1]: p for p in plan["pieces"] if p["kind"] == "centre"}
    os.makedirs(os.path.join(HERE, "trees"), exist_ok=True)
    built = {}
    for spec in specs:
        pid, t = spec.split(":")
        t = float(t)
        p = byid[pid]
        if pid not in built:
            tb = time.process_time()
            built[pid] = (tcc.centre_layer_data(at, ce, p["hi"], label=p["id"], raised=True),
                          time.process_time() - tb)
        ld, build_cpu = built[pid]
        c0 = time.process_time()
        with precision(320):
            res, roots, nodes = certify_logged(tl, ld, p["lo"], p["hi"], target=0.0,
                                               tau_min=t, sigma_fn=tcc.sigma_box_centre,
                                               max_evals=150000)
        cpu = time.process_time() - c0
        part = json.load(open(os.path.join(PARTS, "centre_t%d" % int(t), "centre", pid + ".json")))
        pr = part["result"]
        same = dict(worst=res["worst"] == pr["worst"], boxes=res["boxes"] == pr["boxes"],
                    evals=res["evals"] == pr["evals"], empty=res["empty"] == pr["empty"],
                    depth_hist=res["depth_hist"] == pr["depth_hist"],
                    params=json.loads(json.dumps(res["params"], default=jnum)) ==
                    {k: (v if v != INF else "inf") for k, v in pr["params"].items()})
        rec = dict(schema="review-s46-cover-tree/1", piece=p["id"], lo=p["lo"], hi=p["hi"],
                   tau_min=t, result={k: v for k, v in res.items() if k != "fails"},
                   fails=res["fails"], part_result=dict(worst=pr["worst"], boxes=pr["boxes"],
                                                        evals=pr["evals"], empty=pr["empty"]),
                   same_as_part=same, build_cpu_s=build_cpu, certify_cpu_s=cpu,
                   node_fields=["raw", "nb", "lower", "upper", "outcome", "var", "children"],
                   roots=[[jnum(x) for x in r] for r in roots],
                   nodes=[[[jnum(x) for x in n[0]], ([jnum(x) for x in n[1]] if n[1] else None),
                           n[2], jnum(n[3]), n[4], n[5], n[6]] for n in nodes])
        path = os.path.join(HERE, "trees", "%s_t%d.json.gz" % (pid, int(t)))
        with gzip.open(path, "wt") as fh:
            json.dump(rec, fh, separators=(",", ":"))
        print("%s tau>=%g: ok=%s worst=%r boxes=%d evals=%d empty=%d nodes=%d cpu=%.1fs same=%s"
              % (pid, t, res["ok"], res["worst"], res["boxes"], res["evals"], res["empty"],
                 len(nodes), cpu, same), flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
