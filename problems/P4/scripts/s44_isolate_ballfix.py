#!/usr/bin/env python3
"""Isolate the evaluation-ball fix on centred-form (mv3) layer parts.

The recompute of 2026-10-17 changed two things at once for mv3 parts: the
evaluation balls (ballutil, eba80cb4; bit-neutral wherever no miss occurs)
and the mv3 iota-centre and gradient w-hull (13119db1; changes results by
design).  This script recomputes given pieces on the current code but with
the pre-13119db1 iota-centre and w-hull patched back into tc_layer_mv, with
the engine options of the pre-fix parts, and compares bound, worst, boxes,
evals and ladder attempts with the pre-fix values in
results/theorem_c/layer_recompute/snapshot.json.  Identical results show that
the ball fix alone leaves those parts bit for bit unchanged.  Writes nothing
into parts/; the record goes to --out.

usage: s44_isolate_ballfix.py --ids tube/full-170,gap/gap-000 --out FILE
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import s44_layer as s44  # noqa: E402

INF = float("inf")


def old_centre_patch(mv):
    from flint import acb, arb

    def _centre(g, box):
        tlo, thi = g["tau"]
        ec = (box[2] + box[3]) / 2
        ic = (1.0 / tlo + (0.0 if thi == INF else 1.0 / thi)) / 2
        return ec, arb(ic), arb(1) / arb(ic)

    def _w_dv(g, n, ie, ii):
        Eb, itb, wv = g["E"], g["itau"], g["w"]
        q = acb(1) / acb(Eb * itb, arb(1))
        d = [acb(0)] * n
        d[ie] = -(wv * wv)
        d[ii] = acb(0, 1) * q * q
        return mv.DV(wv, d)
    mv._centre, mv._w_dv = _centre, _w_dv


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True)
    ap.add_argument("--snapshot", default=os.path.join(s44.RECOMPUTE_DIR, "snapshot.json"))
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    snap = json.load(open(a.snapshot))["parts"]
    plan = s44.s43.load_plan()
    byid = {p["id"]: p for p in plan["pieces"]}
    import p4.validated.tc_layer_mv as mv
    old_centre_patch(mv)
    out = dict(schema="s44-isolate-ballfix/1", created=s44.s43.now(),
               git_head=s44.s43.git_head()[0], pieces={})
    ctx = None
    for pid in a.ids.split(","):
        p = byid[pid]
        olds = {dn: snap.get("%s|%s" % (dn, pid)) for dn in ("layer", "layer_t14")}
        eng = (olds["layer"] or olds["layer_t14"])["engine"]
        opts = dict(s44.ENGINE_DEFAULTS, **eng)
        if ctx is None:
            ctx = s44.Ctx(None, opts)
        t0 = time.time()
        ld, la, lb, inputs, call = s44.build_piece(p, ctx, opts["combine"])
        rec = {}
        for dn, tmin in (("layer", 224.0), ("layer_t14", 14.0)):
            o = olds[dn]
            if o is None:
                continue
            r = ctx.tl.run_layer_piece(ld, la, lb, tau_min=tmin, **s44.engine_kwargs(opts))
            new = dict(status="pass" if r["ok"] else "fail", certified_bound=r["certified_bound"],
                       worst=r["worst"], boxes=r["boxes"], evals=r["evals"],
                       attempts=s44._attempts(r["attempts"]))
            oldv = dict((k, o.get(k)) for k in ("status", "certified_bound", "worst", "boxes",
                                                 "evals"))
            oldv["attempts"] = s44._attempts(o.get("attempts"))
            same = {k: oldv[k] == new[k] for k in new}
            rec[dn] = dict(identical=all(same.values()), fields=same, old=oldv, new=new)
            s44.s43.log("%s %s identical=%s bound old %r new %r boxes old %r new %r"
                        % (pid, dn, rec[dn]["identical"], oldv["certified_bound"],
                           new["certified_bound"], oldv["boxes"], new["boxes"]))
        rec["engine"] = opts
        rec["wall_s"] = time.time() - t0
        out["pieces"][pid] = rec
    out["all_identical"] = all(v["identical"] for r in out["pieces"].values()
                               for k, v in r.items() if isinstance(v, dict) and "identical" in v)
    out["code"] = s44.s43.loaded_code()
    s44.s43.write_json_atomic(a.out, out)
    print("all identical:", out["all_identical"], "->", a.out)


if __name__ == "__main__":
    main()
