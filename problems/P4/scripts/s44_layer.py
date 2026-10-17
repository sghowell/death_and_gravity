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
[14, 224]) it closes the piece for all tau >= 14; `s43_parts.py merge` joins
the two by piece id.  This driver never writes S4-3 parts; it reads the plan,
the cache records and the helper functions of s43_parts.py.

Commands:
  list   [selection]                       print the selected piece ids
  work   [selection] [engine options] [--tube-cache FILE] [--force] [--max-pieces N]
                                           compute and write layer parts; one pair
                                           build per piece, one cover per tau_min
  launch --workers N [--nice 10] [--tag T] [selection] [engine options]
                                           detached `work` processes sharing one
                                           lock-guarded queue (logs/s44/)
  status [--tau-min T] [-v]                layer parts vs the band, and the S4-3
                                           state of each piece
selection:      --ids PAT,..  --regions R,..  --all-gap  --beyond-band
                (default: the band = sonic/* and the gap steps with xi_hi <=
                XI_BAND; --beyond-band selects from every sonic and tube-kind
                piece of the plan, e.g. --beyond-band --regions tube,deep)
engine options: --tau-mins T1,T2 (default 224; 14 covers both S4-3 phases;
                each tau_min writes to its own directory: parts/layer for 224,
                parts/layer_t<T> otherwise; --tau-min T is the one-mode alias)
                --norm layer|piece  (margin ladder: absolute -0.75/-0.5/-0.25,
                or fractions 0.5/0.25 of the piece's own float scale)
                --combine           (tube steps: derivative blocks combined at
                the polynomial level, tc_layer.LayerData(combine=True))
                --mv mv|mv3         (tube steps: centred enclosure of Sigma from
                p4.validated.tc_layer_mv, in (eps, 1/tau) or (sigma, eps, 1/tau))
                --piece-ladder-evals N --piece-final-evals N  (norm=piece budgets)
A part is skipped on re-run when it is valid: schema, spec equal to the plan
entry, every recorded source hash current, and the same engine options.
Failures are recorded (status fail/error) and also skipped unless --force.
"""
from __future__ import annotations

import argparse
import fnmatch
import json
import math
import os
import signal
import socket
import subprocess
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
LOGS = os.path.join(s43.P4, "logs", "s44")
LOCKS = os.path.join(LOGS, "locks")


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


def in_band(p):
    return p["region"] == "sonic" or (p["region"] == "gap" and xi_range(p)[1] <= XI_BAND)


def select(plan, args):
    if getattr(args, "beyond_band", False):
        ps = [p for p in plan["pieces"] if p["kind"] in ("sonic", "tube")]
    else:
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


def tau_mins_of(args):
    if getattr(args, "tau_mins", None):
        out = [float(v) for v in args.tau_mins.split(",") if v.strip()]
    else:
        out = [float(getattr(args, "tau_min", TAU_MIN))]
    if not out or len(set(out)) != len(out) or not all(t > 0 for t in out):
        raise SystemExit("--tau-mins: distinct positive values expected")
    return out


ENGINE_DEFAULTS = dict(norm="layer", combine=False, mv=None, piece_ladder_evals=None,
                       piece_final_evals=None)


def engine_opts(args):
    out = dict(ENGINE_DEFAULTS)
    out.update(norm=getattr(args, "norm", "layer") or "layer",
               combine=bool(getattr(args, "combine", False)),
               mv=getattr(args, "mv", None) or None,
               piece_ladder_evals=getattr(args, "piece_ladder_evals", None),
               piece_final_evals=getattr(args, "piece_final_evals", None))
    return out


def engine_kwargs(opts):
    """run_layer_piece keyword arguments of the engine options."""
    kw = dict(norm=opts["norm"])
    if opts.get("mv"):
        import p4.validated.tc_layer_mv as tmv
        kw["sigma_fn"] = getattr(tmv, "sigma_box_" + opts["mv"])
    for k in ("piece_ladder_evals", "piece_final_evals"):
        if opts.get(k):
            kw[k] = int(opts[k])
    return kw


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
        self._tubes = {}
        s43.log("context ready (%.1f s)" % (time.time() - t0))

    def tube(self, name):
        """(Tube, content hash, path) of a tube cache, pinned by its record."""
        if name not in self._tubes:
            from p4.validated.lintube import Tube
            path = (self.tube_cache if (self.tube_cache and name == "gap")
                    else s43.cache_path(name))
            rec = json.load(open(s43.cache_record_path(name)))
            h = s43.tube_content_hash(path)
            if h != rec["content_sha256"]:
                raise RuntimeError("%s tube %s: content hash %s differs from the record %s"
                                   % (name, path, h, rec["content_sha256"]))
            self._tubes[name] = (Tube.load(path), h, path)
        return self._tubes[name]

    def gap_tube(self):
        return self.tube("gap")


def build_piece(p, ctx, combine=False):
    """LayerData of a plan piece (built exactly as the S4-3 region code builds
    its pairs) plus (la, lb), the inputs record and the region call."""
    tl = ctx.tl
    if p["kind"] == "sonic":
        if combine:
            raise RuntimeError("--combine is for tube steps (no Laurent offsets)")
        ld = tl.sonic_layer_data(ctx.at, ctx.ex, p["hi"], label=p["id"])
        la, lb = p["lo"], p["hi"]
        inputs = dict(kind="sonic", series="tc_energy.sonic_expansion_cached(320)",
                      tm_domain=[0.0, p["hi"]])
        call = ("tc_layer.sonic_layer_data(at, ex, %r) ; run_layer_piece(ld, %r, %r)"
                % (p["hi"], la, lb))
    elif p["kind"] == "tube":
        tube, h, path = ctx.tube(p["cache"])
        if h != p["cache_content_sha256"]:
            raise RuntimeError("%s: cache hash %s differs from the plan entry" % (p["id"], h))
        steps = [sd for sd in tube.steps
                 if float(sd.x) - sd.h >= p["x_min"] - 1e-12 and float(sd.x) <= p["x_max"] + 1e-12]
        sd = steps[p["index"]]
        if float(sd.x) != p["x"] or sd.h != p["h"]:
            raise RuntimeError("step mismatch: plan (%r, %r) vs cache (%r, %r)"
                               % (p["x"], p["h"], float(sd.x), sd.h))
        ld = tl.tube_layer_data(ctx.at, sd, label=p["id"], combine=combine)
        la, lb = 0.0, float(sd.h)
        inputs = dict(kind="tube", cache=p["cache"], cache_content_sha256=h,
                      cache_file=s43.rel(path))
        call = ("tc_layer.tube_layer_data(at, tube_%s.steps[%d]%s) ; run_layer_piece(ld, 0, %r)"
                % (p["cache"], p["index"], ", combine=True" if combine else "", lb))
    else:
        raise RuntimeError("unsupported piece kind for the layer engine: %s" % p["id"])
    return ld, la, lb, inputs, call


def part_valid(d, p, opts):
    """A stored part may be skipped: right schema and spec, current code, the
    same engine options (a part without the record used the defaults)."""
    if d.get("schema") != SCHEMA or d.get("spec") != s43.spec_of(p):
        return False
    if not s43.code_matches(d["code"]["files"]):
        return False
    have = dict(ENGINE_DEFAULTS)
    have.update(d.get("engine") or {})
    return have == opts


def claim(pid):
    """O_EXCL lock in logs/s44/locks (stale locks of dead local workers are
    reclaimed), so several `work` processes can share one queue."""
    os.makedirs(LOCKS, exist_ok=True)
    lk = os.path.join(LOCKS, pid.replace("/", "__") + ".lock")
    payload = json.dumps(dict(host=socket.gethostname(), pid=os.getpid(), t=s43.now()))
    for _ in range(2):
        try:
            fd = os.open(lk, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
            os.write(fd, payload.encode())
            os.close(fd)
            return lk
        except FileExistsError:
            try:
                info = json.load(open(lk))
            except Exception:
                return None
            if info.get("host") != socket.gethostname():
                return None
            try:
                os.kill(int(info["pid"]), 0)
                return None
            except ProcessLookupError:
                os.remove(lk)
            except PermissionError:
                return None
    return None


def write_part(p, mode_tau, res, err, inputs, call, timing, opts, build_shared):
    files, (head, dirty) = s43.loaded_code(), s43.git_head()
    xi = xi_range(p)
    status = "error" if err else ("pass" if res["ok"] else "fail")
    rec = dict(schema=SCHEMA, id=p["id"], region=p["region"], spec=s43.spec_of(p),
               status=status,
               scope=dict(tau=[mode_tau, "inf"], eps=[0.0, 6.0], xi=list(xi),
                          statement=("for all xi in the piece, eps in [0, 6], "
                                     "tau >= %g: dT, dB, nB != 0 and "
                                     "Sigma = S xi^2/tau^2 <= certified_bound < 0"
                                     % mode_tau)),
               certified_bound=(res or {}).get("certified_bound"),
               engine=opts,
               result=_clean({k: v for k, v in (res or {}).items()}),
               error=err, inputs=inputs, region_call=call,
               code=dict(files=files, fingerprint=s43.fingerprint(files),
                         git_head=head, git_dirty_src=dirty),
               driver_sha256=s43.sha256_file(os.path.abspath(__file__)),
               env=s43.env_info(),
               timing=dict(timing, build_shared_with=build_shared),
               provenance="computed by problems/P4/scripts/s44_layer.py work")
    s43.write_json_atomic(part_path(p["id"], mode_tau), _clean(rec))
    return rec


def cmd_work(args):
    plan = s43.load_plan()
    ps = select(plan, args)
    taus = tau_mins_of(args)
    opts = engine_opts(args)
    ctx = None
    done = 0
    s43.log("layer worker start: %d candidate pieces, tau_min %s, engine %s"
            % (len(ps), taus, opts))

    def todo(p):
        out = []
        for t in taus:
            path = part_path(p["id"], t)
            if os.path.exists(path) and not args.force:
                try:
                    d = json.load(open(path))
                except Exception:
                    d = None
                if d is not None and part_valid(d, p, opts):
                    continue
            out.append(t)
        return out

    for p in ps:
        if args.max_pieces and done >= args.max_pieces:
            break
        if not todo(p):
            continue
        lk = claim("layer/" + p["id"])
        if lk is None:
            continue
        try:
            modes = todo(p)                       # re-check under the lock
            if not modes:
                continue
            if ctx is None:
                ctx = Ctx(args.tube_cache)
            s43.log("start %s (tau_min %s)" % (p["id"], modes))
            start = s43.now()
            t0, c0 = time.time(), time.process_time()
            try:
                if p["kind"] == "sonic" and opts.get("mv"):
                    raise RuntimeError("--mv is for tube steps (sonic pieces use the "
                                       "natural layer form)")
                ld, la, lb, inputs, call = build_piece(p, ctx, opts["combine"])
                berr = None
            except Exception:
                import traceback
                ld, la, lb, inputs, call = None, None, None, None, None
                berr = traceback.format_exc()
            tb, cb = time.time() - t0, time.process_time() - c0
            shared = [t for t in modes] if len(modes) > 1 else None
            for t in modes:
                t1, c1 = time.time(), time.process_time()
                res, err = None, berr
                if err is None:
                    try:
                        res = ctx.tl.run_layer_piece(ld, la, lb, tau_min=t,
                                                     **engine_kwargs(opts))
                    except Exception:                 # recorded, never a pass
                        import traceback
                        res, err = None, traceback.format_exc()
                tc_, cc = time.time() - t1, time.process_time() - c1
                timing = dict(start=start, end=s43.now(), build_wall_s=tb, build_cpu_s=cb,
                              certify_wall_s=tc_, certify_cpu_s=cc,
                              wall_s=tb + tc_, cpu_s=cb + cc)
                rec = write_part(p, t, res, err, inputs, call, timing, opts, shared)
                done += 1
                s43.log("done  %s tau>=%g status=%s bound=%s scale=%s boxes=%s build=%.1fs "
                        "cover=%.1fs" % (p["id"], t, rec["status"], rec["certified_bound"],
                                         (res or {}).get("scale"), (res or {}).get("boxes"),
                                         tb, tc_))
        finally:
            try:
                os.remove(lk)
            except OSError:
                pass
    s43.log("layer worker exit: %d parts written" % done)


def spawn_detached(argv, logfile, nice=10):
    os.makedirs(os.path.dirname(logfile), exist_ok=True)
    out = open(logfile, "ab")

    def pre():
        try:
            os.nice(nice)
        except OSError:
            pass
        signal.signal(signal.SIGHUP, signal.SIG_IGN)
    proc = subprocess.Popen(["nohup", sys.executable, os.path.abspath(__file__)] + argv,
                            cwd=s43.ROOT, stdin=subprocess.DEVNULL, stdout=out, stderr=out,
                            start_new_session=True, preexec_fn=pre,
                            env=dict(os.environ, PYTHONUNBUFFERED="1"))
    return proc.pid


def work_argv(args):
    av = ["work"]
    for k in ("ids", "regions", "tube_cache", "tau_mins", "norm", "mv",
              "piece_ladder_evals", "piece_final_evals"):
        v = getattr(args, k, None)
        if v:
            av.append("--%s=%s" % (k.replace("_", "-"), v))
    for k in ("all_gap", "beyond_band", "combine", "force"):
        if getattr(args, k, False):
            av.append("--" + k.replace("_", "-"))
    if getattr(args, "max_pieces", 0):
        av += ["--max-pieces", str(args.max_pieces)]
    if not getattr(args, "tau_mins", None):
        av.append("--tau-min=%g" % args.tau_min)
    return av


def cmd_launch(args):
    tag = args.tag or time.strftime("s44-%Y%m%d-%H%M%S")
    pids = []
    av = work_argv(args)
    for k in range(args.workers):
        lf = os.path.join(LOGS, "%s-w%d.log" % (tag, k))
        pids.append((spawn_detached(av, lf, args.nice), s43.rel(lf)))
    rec = dict(tag=tag, host=socket.gethostname(), started=s43.now(), nice=args.nice,
               workers=[dict(pid=p, log=l) for p, l in pids], argv=av)
    s43.write_json_atomic(os.path.join(LOGS, "launch-%s.json" % tag), rec)
    for p, l in pids:
        print("pid %d  log %s" % (p, l))


def cmd_list(args):
    plan = s43.load_plan()
    ps = select(plan, args)
    for p in ps:
        print(p["id"])
    print("# %d pieces" % len(ps), file=sys.stderr)


def cmd_status(args):
    plan = s43.load_plan()
    band = band_pieces(plan)
    rows = dict(pass_=0, fail=0, error=0, stale=0, missing=0)
    worst = None
    cpu = 0.0
    for p in band:
        path = part_path(p["id"], args.tau_min)
        d = json.load(open(path)) if os.path.exists(path) else None
        st43 = s43.part_state(p, s43.read_part(p["id"]))
        if d is None:
            rows["missing"] += 1
            st = "missing"
        elif d.get("spec") != s43.spec_of(p) or not s43.code_matches(d["code"]["files"]):
            rows["stale"] += 1
            st = "stale"
        else:
            st = d["status"]
            rows["pass_" if st == "pass" else st] += 1
            if d.get("certified_bound") is not None:
                worst = d["certified_bound"] if worst is None else max(worst, d["certified_bound"])
            cpu += (d.get("timing") or {}).get("cpu_s", 0.0)
        if args.verbose or st not in ("missing", "pass"):
            print("%-14s layer=%-8s bound=%-10s s43=%s"
                  % (p["id"], st, None if d is None else d.get("certified_bound"), st43))
    print("band (tau >= %g): %d pieces; layer pass %d, fail %d, error %d, stale %d, "
          "missing %d; worst certified bound %s; cpu %.1f s"
          % (args.tau_min, len(band), rows["pass_"], rows["fail"], rows["error"],
             rows["stale"], rows["missing"], worst, cpu))
    # layer parts outside the band (e.g. the deep probe)
    extra = []
    byid = {p["id"]: p for p in plan["pieces"]}
    root = layer_dir(args.tau_min)
    for dp, _, fs in os.walk(root):
        for f in sorted(fs):
            if not f.endswith(".json"):
                continue
            pid = os.path.relpath(os.path.join(dp, f), root)[:-5]
            if pid in byid and not in_band(byid[pid]):
                extra.append((pid, json.load(open(os.path.join(dp, f)))))
    for pid, d in sorted(extra):
        r = d.get("result") or {}
        print("beyond band: %-16s status=%s bound=%s scale=%s engine=%s cpu=%.0fs"
              % (pid, d["status"], d.get("certified_bound"), r.get("scale"),
                 d.get("engine"), (d.get("timing") or {}).get("cpu_s", 0.0)))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    def sel(q):
        q.add_argument("--ids", help="comma-separated fnmatch patterns on piece ids")
        q.add_argument("--regions", help="sonic,gap (with --beyond-band also tube,deep,band)")
        q.add_argument("--all-gap", action="store_true",
                       help="all 107 gap steps instead of the band steps 0-88")
        q.add_argument("--beyond-band", action="store_true",
                       help="select from every sonic and tube-kind piece of the plan")

    def eng(q):
        q.add_argument("--tube-cache", help="gap tube file (content hash must match the record)")
        q.add_argument("--force", action="store_true")
        q.add_argument("--max-pieces", type=int, default=0)
        q.add_argument("--tau-min", type=float, default=TAU_MIN,
                       help="lower end of the tau range (default 224; 14 covers both "
                            "S4-3 phases; parts then go to parts/layer_t14)")
        q.add_argument("--tau-mins", help="several tau_min values, one pair build: e.g. 224,14")
        q.add_argument("--norm", default="layer", choices=("layer", "piece"))
        q.add_argument("--combine", action="store_true",
                       help="tube steps: polynomial-level derivative blocks")
        q.add_argument("--mv", choices=("mv", "mv3"),
                       help="tube steps: centred enclosure (p4.validated.tc_layer_mv)")
        q.add_argument("--piece-ladder-evals", type=int)
        q.add_argument("--piece-final-evals", type=int)

    q = sub.add_parser("list")
    sel(q)
    q.set_defaults(fn=cmd_list)
    q = sub.add_parser("work")
    sel(q)
    eng(q)
    q.set_defaults(fn=cmd_work)
    q = sub.add_parser("launch")
    sel(q)
    eng(q)
    q.add_argument("--workers", type=int, required=True)
    q.add_argument("--nice", type=int, default=10)
    q.add_argument("--tag")
    q.set_defaults(fn=cmd_launch)
    q = sub.add_parser("status")
    q.add_argument("-v", "--verbose", action="store_true")
    q.add_argument("--tau-min", type=float, default=TAU_MIN)
    q.set_defaults(fn=cmd_status)
    args = ap.parse_args(argv)
    args.fn(args)


if __name__ == "__main__":
    main()
