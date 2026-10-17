#!/usr/bin/env python3
"""S4-4 layer-engine driver: tau >= 224 on the near-sonic band (p4.validated.tc_layer).

The S4-3 commutator certificate (scripts/s43_parts.py) leaves tau > 224 open on
the near-sonic band: the 51 sonic pieces and the gap-tube steps 0-88
(xi <= 0.0357), 140 pieces in all.  For each such piece of
results/theorem_c/parts/plan.json this driver builds the Taylor-model pairs
exactly as the S4-3 region code does (tc_commutator.run_sonic / run_tube),
runs tc_layer.run_layer_piece, and writes

  problems/P4/results/theorem_c/parts/layer/<id>.json      (schema s44-layer/1)

A passing layer part certifies, for every xi in the piece, every
eps in [0, 6] and every tau >= 224 (tau = infinity included through the
pseudo-boxes):  the denominators dT, dB and nB of the commuted problem are
nonzero, and  Sigma = S xi^2 / tau^2 <= certified_bound < 0,  hence S < 0.
Together with the S4-3 part of the same piece (status pass-direct: tau in
[14, 224]) it closes the piece for all tau >= 14.  This driver does not modify
s43_parts.py or its parts; it only reads the plan, the cache records and the
helper functions.

Commands:
  list   [--ids PAT,..] [--regions R,..] [--all-gap]   print the selected piece ids
  work   [--ids PAT,..] [--regions R,..] [--all-gap] [--tube-cache FILE] [--force]
         [--max-pieces N] [--tau-min T]                compute and write layer parts
  status [--tau-min T]                                 layer parts vs the band, and
                                                       the S4-3 status of each piece
--tau-min 14 runs the engine on tau >= 14 (both S4-3 phases) and writes to
parts/layer_t14/ instead.
Selection defaults to the band (sonic/* and gap steps with xi_hi <= XI_BAND).
"""
from __future__ import annotations

import argparse
import fnmatch
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import s43_parts as s43  # noqa: E402  (plan, cache records, hashing helpers)

LAYER = os.path.join(s43.PARTS, "layer")
SCHEMA = "s44-layer/1"
TAU_MIN = 224.0                  # the S4-4 item: tau > 224 on the band
XI_BAND = 0.0357 + 1e-9          # gap steps 0-88: the tau-quadratic phase fails there


def layer_dir(tau_min=TAU_MIN):
    """parts/layer for the default tau >= 224; parts/layer_t<tau_min> for a
    lower tau_min (e.g. layer_t14: the engine covering both S4-3 phases)."""
    if tau_min == TAU_MIN:
        return LAYER
    return os.path.join(s43.PARTS, "layer_t%g" % tau_min)


def part_path(pid, tau_min=TAU_MIN):
    return os.path.join(layer_dir(tau_min), pid + ".json")


def xi_range(p):
    a, b = s43.x_interval(p)
    return (-b, -a)


def band_pieces(plan, all_gap=False):
    out = []
    for p in plan["pieces"]:
        if p["region"] == "sonic":
            out.append(p)
        elif p["region"] == "gap" and (all_gap or xi_range(p)[1] <= XI_BAND):
            out.append(p)
    return out


def select(plan, args):
    ps = band_pieces(plan, getattr(args, "all_gap", False))
    if getattr(args, "regions", None):
        want = set(args.regions.split(","))
        ps = [p for p in ps if p["region"] in want]
    if getattr(args, "ids", None):
        pats = args.ids.split(",")
        ps = [p for p in ps if any(fnmatch.fnmatch(p["id"], q) for q in pats)]
    return ps


def _clean(x):
    """JSON-safe copy (inf -> "inf", tuples -> lists)."""
    if isinstance(x, float):
        return x if math.isfinite(x) else ("inf" if x > 0 else "-inf")
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    return x


class Ctx:
    def __init__(self, tube_cache=None):
        t0 = time.time()
        s43.import_p4()
        import p4.validated.tc_energy as te
        import p4.validated.tc_layer as tl
        self.te, self.tl = te, tl
        self.at = te.Atoms(full=True)
        self.ex = te.sonic_expansion_cached(320)
        self.tube_cache = tube_cache
        self._tube = None
        s43.log("context ready (%.1f s)" % (time.time() - t0))

    def gap_tube(self):
        if self._tube is None:
            from p4.validated.lintube import Tube
            path = self.tube_cache or s43.cache_path("gap")
            rec = json.load(open(s43.cache_record_path("gap")))
            h = s43.tube_content_hash(path)
            if h != rec["content_sha256"]:
                raise RuntimeError("gap tube %s: content hash %s differs from the record %s"
                                   % (path, h, rec["content_sha256"]))
            self._tube = (Tube.load(path), h, path)
        return self._tube


def run_piece(p, ctx, tau_min=TAU_MIN):
    tl = ctx.tl
    t0, c0 = time.time(), time.process_time()
    if p["kind"] == "sonic":
        ld = tl.sonic_layer_data(ctx.at, ctx.ex, p["hi"], label=p["id"])
        la, lb = p["lo"], p["hi"]
        inputs = dict(kind="sonic", series="tc_energy.sonic_expansion_cached(320)",
                      tm_domain=[0.0, p["hi"]])
        call = ("tc_layer.sonic_layer_data(at, ex, %r) ; run_layer_piece(ld, %r, %r)"
                % (p["hi"], la, lb))
    elif p["kind"] == "tube" and p["cache"] == "gap":
        tube, h, path = ctx.gap_tube()
        steps = [sd for sd in tube.steps
                 if float(sd.x) - sd.h >= p["x_min"] - 1e-12 and float(sd.x) <= p["x_max"] + 1e-12]
        sd = steps[p["index"]]
        if float(sd.x) != p["x"] or sd.h != p["h"]:
            raise RuntimeError("step mismatch: plan (%r, %r) vs cache (%r, %r)"
                               % (p["x"], p["h"], float(sd.x), sd.h))
        ld = tl.tube_layer_data(ctx.at, sd, label=p["id"])
        la, lb = 0.0, float(sd.h)
        inputs = dict(kind="tube", cache="gap", cache_content_sha256=h,
                      cache_file=s43.rel(path))
        call = "tc_layer.tube_layer_data(at, tube_gap.steps[%d]) ; run_layer_piece(ld, 0, %r)" \
            % (p["index"], lb)
    else:
        raise RuntimeError("unsupported piece kind for the layer engine: %s" % p["id"])
    tb, cb = time.time() - t0, time.process_time() - c0
    t1, c1 = time.time(), time.process_time()
    res = tl.run_layer_piece(ld, la, lb, tau_min=tau_min)
    tc_, cc = time.time() - t1, time.process_time() - c1
    timing = dict(build_wall_s=tb, build_cpu_s=cb, certify_wall_s=tc_, certify_cpu_s=cc,
                  wall_s=tb + tc_, cpu_s=cb + cc)
    return res, inputs, call, timing


def cmd_list(args):
    plan = s43.load_plan()
    ps = select(plan, args)
    for p in ps:
        print(p["id"])
    print("# %d pieces" % len(ps), file=sys.stderr)


def cmd_work(args):
    plan = s43.load_plan()
    ps = select(plan, args)
    ctx = None
    done = 0
    for p in ps:
        if args.max_pieces and done >= args.max_pieces:
            break
        path = part_path(p["id"], args.tau_min)
        if os.path.exists(path) and not args.force:
            d = json.load(open(path))
            if d.get("schema") == SCHEMA and d.get("spec") == s43.spec_of(p) and \
                    s43.code_matches(d["code"]["files"]):
                continue
        if ctx is None:
            ctx = Ctx(args.tube_cache)
        s43.log("start %s" % p["id"])
        start = s43.now()
        try:
            res, inputs, call, timing = run_piece(p, ctx, args.tau_min)
            status = "pass" if res["ok"] else "fail"
            err = None
        except Exception as e:                       # recorded, never a pass
            import traceback
            res, inputs, call, timing = None, None, None, None
            status, err = "error", traceback.format_exc()
        files, (head, dirty) = s43.loaded_code(), s43.git_head()
        xi = xi_range(p)
        rec = dict(schema=SCHEMA, id=p["id"], region=p["region"], spec=s43.spec_of(p),
                   status=status,
                   scope=dict(tau=[args.tau_min, "inf"], eps=[0.0, 6.0], xi=list(xi),
                              statement=("for all xi in the piece, eps in [0, 6], "
                                         "tau >= %g: dT, dB, nB != 0 and "
                                         "Sigma = S xi^2/tau^2 <= certified_bound < 0"
                                         % args.tau_min)),
                   certified_bound=(res or {}).get("certified_bound"),
                   result=_clean({k: v for k, v in (res or {}).items()}),
                   error=err, inputs=inputs, region_call=call,
                   code=dict(files=files, fingerprint=s43.fingerprint(files),
                             git_head=head, git_dirty_src=dirty),
                   driver_sha256=s43.sha256_file(os.path.abspath(__file__)),
                   env=s43.env_info(),
                   timing=dict(start=start, end=s43.now(), **(timing or {})),
                   provenance="computed by problems/P4/scripts/s44_layer.py work")
        s43.write_json_atomic(path, _clean(rec))
        done += 1
        s43.log("done  %s status=%s bound=%s boxes=%s wall=%.1fs"
                % (p["id"], status, rec["certified_bound"],
                   (res or {}).get("boxes"), (timing or {}).get("wall_s", 0.0)))


def cmd_status(args):
    plan = s43.load_plan()
    band = band_pieces(plan)
    rows = dict(pass_=0, fail=0, error=0, missing=0)
    worst = None
    cpu = 0.0
    for p in band:
        path = part_path(p["id"], args.tau_min)
        d = json.load(open(path)) if os.path.exists(path) else None
        st43 = s43.part_state(p, s43.read_part(p["id"]))
        if d is None:
            rows["missing"] += 1
            st = "missing"
        else:
            st = d["status"]
            rows["pass_" if st == "pass" else st] += 1
            if d.get("certified_bound") is not None:
                worst = d["certified_bound"] if worst is None else max(worst, d["certified_bound"])
            cpu += (d.get("timing") or {}).get("cpu_s", 0.0)
        if args.verbose or st != "missing":
            print("%-14s layer=%-8s bound=%-10s s43=%s"
                  % (p["id"], st, None if d is None else d.get("certified_bound"), st43))
    print("band: %d pieces; layer pass %d, fail %d, error %d, missing %d; "
          "worst certified bound %s; cpu %.1f s"
          % (len(band), rows["pass_"], rows["fail"], rows["error"], rows["missing"],
             worst, cpu))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    def sel(q):
        q.add_argument("--ids", help="comma-separated fnmatch patterns on piece ids")
        q.add_argument("--regions", help="sonic,gap")
        q.add_argument("--all-gap", action="store_true",
                       help="all 107 gap steps instead of the band steps 0-88")

    q = sub.add_parser("list")
    sel(q)
    q.set_defaults(fn=cmd_list)
    q = sub.add_parser("work")
    sel(q)
    q.add_argument("--tube-cache", help="gap tube file (content hash must match the record)")
    q.add_argument("--force", action="store_true")
    q.add_argument("--max-pieces", type=int, default=0)
    q.add_argument("--tau-min", type=float, default=TAU_MIN,
                   help="lower end of the tau range (default 224; 14 covers both "
                        "S4-3 phases; parts then go to parts/layer_t14)")
    q.set_defaults(fn=cmd_work)
    q = sub.add_parser("status")
    q.add_argument("-v", "--verbose", action="store_true")
    q.add_argument("--tau-min", type=float, default=TAU_MIN)
    q.set_defaults(fn=cmd_status)
    args = ap.parse_args(argv)
    args.fn(args)


if __name__ == "__main__":
    main()
