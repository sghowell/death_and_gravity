#!/usr/bin/env python3
"""S4-6: the centre variant of the layer engine on the 30 centre pieces.

The centre sigma-function is p4.validated.tc_layer_centre.sigma_box_centre:
graded pivot forms in y (z = y kappa, w = 1/kappa), Sigma_c = S/(1 + s^2),
s = y tau, on layer data whose nB, nB1, nB2 (and so nY) kappa^1..kappa^7
coefficients are raised from Laurent offset 8 to 10 by the exact centre layer
identity (tc_layer_centre.check_nb_identity).  It is plugged into the
unchanged tc_layer.run_layer_piece.  See notes/s4-6-centre.md.

usage:
  s46_centre.py identity [--out FILE]        the exact nB layer identity
  s46_centre.py work [--pieces c29,...] [--tau-mins 224,14] [--max-evals N]
                     [--force]               compute parts (O_EXCL locks per piece)
  s46_centre.py launch --workers N [work options] [--nice 10] [--tag T]
  s46_centre.py status                       parts summary
  s46_centre.py grid [--pieces c29,c00] [--out FILE]   float design grid

Parts: problems/P4/results/theorem_c/parts/centre_t<tau_min>/centre/<cXX>.json
with the import-time source hashes, the write-time check and the driver pin
(s43_parts.code_record).  Launch records, logs and locks: problems/P4/logs/s46/
(machine-local).
"""
from __future__ import annotations

import argparse
import fnmatch
import json
import os
import signal
import socket
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import s43_parts as s43  # noqa: E402

DRIVER = os.path.abspath(__file__)
DRIVER_SHA256 = s43.sha256_file(DRIVER)
LOGS = os.path.join(s43.P4, "logs", "s46")
LOCKS = os.path.join(LOGS, "locks")
IDENTITY = os.path.join(s43.TC_DIR, "centre_nb_identity.json")
GRID = os.path.join(s43.TC_DIR, "centre_grid.json")
SCHEMA = "s46-centre/1"
EVALUATOR = "sigma_box_centre"
INF = float("inf")


def part_dir(tau_min):
    return os.path.join(s43.PARTS, "centre_t%d" % int(tau_min))


def part_path(pid, tau_min):
    return os.path.join(part_dir(tau_min), pid + ".json")


def _clean(x):
    if isinstance(x, float):
        if x != x:
            return "nan"
        if x in (INF, -INF):
            return "inf" if x > 0 else "-inf"
        return x
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    return x


class Ctx:
    def __init__(self):
        t0 = time.time()
        s43.import_p4()
        import p4.validated.tc_energy as te
        import p4.validated.tc_layer as tl
        import p4.validated.tc_layer_centre as tcc
        from p4.validated.ec_centre import certified_centre_expansion
        s43.snapshot_code()                      # the sources as loaded
        self.te, self.tl, self.tcc = te, tl, tcc
        self.at = te.Atoms(full=True)
        # the certified (a*, mu*) balls (ec_centre), as tc_commutator._centre_expansion
        self.ce = certified_centre_expansion(256, K=30)
        self.identity = tcc.check_nb_identity(self.at)
        if not self.identity["ok"]:
            raise SystemExit("the centre nB layer identity failed")
        self.setup_s = time.time() - t0

    def build(self, p):
        t0, c0 = time.time(), time.process_time()
        ld = self.tcc.centre_layer_data(self.at, self.ce, p["hi"], label=p["id"], raised=True)
        return ld, time.time() - t0, time.process_time() - c0


def centre_pieces(plan, pattern=None):
    ps = [p for p in plan["pieces"] if p["kind"] == "centre"]
    if pattern:
        pats = [x.strip() for x in pattern.split(",") if x.strip()]
        ps = [p for p in ps if any(fnmatch.fnmatch(p["id"].split("/")[1], q) or
                                   fnmatch.fnmatch(p["id"], q) for q in pats)]
    # the largest y first: those pieces are the most expensive
    return sorted(ps, key=lambda p: -p["index"])


# -------------------------------------------------------------- identity
def cmd_identity(args):
    ctx = Ctx()
    idn = ctx.identity
    rec = dict(schema="s46-centre-nb-identity/1", created=s43.now(), ok=idn["ok"],
               identity=_clean(idn),
               use=("tc_layer_centre.raise_nb: nB, nB1, nB2 kappa^1..kappa^7 coefficients "
                    "re-expressed at Laurent offset 10 (TM.drop_lead(2): the dropped "
                    "coefficients enclose 0 and lie in the exact-coefficient range); nY "
                    "follows (nY = nB - u2 kappa^2 dB with u2 dB_{j-2} of offset >= 10)"),
               method=("exact fmpq arithmetic on the compiled atoms of tc_energy.Atoms, "
                       "centre grading p = -e_N + 2 e_W + e_V with A -> 1, as E2 "
                       "(tc_commutator.check_layer_identity)"),
               consistent_with=("certificates/theorem_c/review_s45/centre_layers.json: "
                                "val(U6) = 10, val(denU) = 8, no odd layers of U6"),
               setup_s=ctx.setup_s, env=s43.env_info(),
               code=s43.code_record({DRIVER: DRIVER_SHA256}), driver_sha256=DRIVER_SHA256)
    out = args.out or IDENTITY
    s43.write_json_atomic(out, rec)
    print("identity ok=%s -> %s" % (idn["ok"], s43.rel(out)))


# ------------------------------------------------------------------ work
def claim(name):
    os.makedirs(LOCKS, exist_ok=True)
    lk = os.path.join(LOCKS, name.replace("/", "__") + ".lock")
    payload = json.dumps(dict(host=socket.gethostname(), pid=os.getpid(), t=s43.now()))
    for _ in range(2):
        try:
            fd = os.open(lk, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            with os.fdopen(fd, "w") as fh:
                fh.write(payload)
            return lk
        except FileExistsError:
            try:
                d = json.load(open(lk))
                if d.get("host") == socket.gethostname():
                    os.kill(int(d["pid"]), 0)
                    return None             # held by a live local worker
            except ProcessLookupError:
                os.remove(lk)               # stale lock of a dead local worker
                continue
            except Exception:
                return None
            return None
    return None


def part_valid(d, p, tau_min):
    return (d.get("schema") == SCHEMA and d.get("id") == p["id"]
            and d.get("spec") == s43.spec_of(p) and float(d.get("tau_min", -1)) == tau_min
            and d.get("driver_sha256") == DRIVER_SHA256 and s43.write_check_ok(d)
            and s43.code_matches(d["code"]["files"]))


def cmd_work(args):
    plan = s43.load_plan()
    ps = centre_pieces(plan, args.pieces)
    taus = [float(t) for t in args.tau_mins.split(",")]
    ctx = None
    done = 0

    def todo(p):
        out = []
        for t in taus:
            path = part_path(p["id"], t)
            if os.path.exists(path) and not args.force:
                try:
                    if part_valid(json.load(open(path)), p, t):
                        continue
                except Exception:
                    pass
            out.append(t)
        return out

    s43.log("centre worker start: %d candidate pieces, tau_min %s" % (len(ps), taus))
    for p in ps:
        if not todo(p):
            continue
        lk = claim("centre/" + p["id"])
        if lk is None:
            continue
        try:
            modes = todo(p)
            if not modes:
                continue
            if ctx is None:
                ctx = Ctx()
            start = s43.now()
            ld, tb, cb = ctx.build(p)
            s43.log("built %s in %.1f s" % (p["id"], tb))
            for t in modes:
                t1, c1 = time.time(), time.process_time()
                res, err = None, None
                try:
                    from p4.validated.arbseries import precision
                    with precision(320):
                        res = ctx.tl.run_layer_piece(
                            ld, p["lo"], p["hi"], tau_min=t,
                            sigma_fn=ctx.tcc.sigma_box_centre, ladder=(),
                            max_evals=args.max_evals)
                except Exception:
                    import traceback
                    err = traceback.format_exc()
                timing = dict(start=start, end=s43.now(), build_wall_s=tb, build_cpu_s=cb,
                              certify_wall_s=time.time() - t1,
                              certify_cpu_s=time.process_time() - c1)
                rec = write_part(ctx, ld, p, t, res, err, timing, args)
                done += 1
                s43.log("done  %s tau>=%g status=%s bound=%s boxes=%s evals=%s %.1fs"
                        % (p["id"], t, rec["status"], rec["certified_bound"],
                           (res or {}).get("boxes"), (res or {}).get("evals"),
                           timing["certify_wall_s"]))
        finally:
            try:
                os.remove(lk)
            except OSError:
                pass
    s43.log("centre worker exit: %d parts written" % done)


def write_part(ctx, ld, p, t, res, err, timing, args):
    status = "error" if err else ("pass" if res["ok"] else "fail")
    tcc = ctx.tcc
    rec = dict(
        schema=SCHEMA, id=p["id"], region="centre", spec=s43.spec_of(p), status=status,
        tau_min=t,
        scope=dict(y=[p["lo"], p["hi"]], eps=[0.0, 6.0], tau=[t, "inf"],
                   statement=("for all y in (0, %r] with y >= %r, eps in [0, 6], tau >= %g: "
                              "dT, dB, nB != 0 and Sigma_c = S/(1 + (y tau)^2) <= "
                              "certified_bound < 0" % (p["hi"], p["lo"], t))),
        certified_bound=(res or {}).get("certified_bound"),
        result=_clean({k: v for k, v in (res or {}).items() if k != "forms"}),
        forms=_clean((res or {}).get("forms")), error=err,
        engine=dict(evaluator=EVALUATOR, cover="tc_layer.run_layer_piece", ladder=[],
                    max_evals=args.max_evals, raised=True, kmax=tcc.KMAX,
                    s_split=tcc.S_SPLIT,
                    centring={"%s/%s" % k: v for k, v in ld.centring.items()}),
        inputs=dict(centre_series="ec_centre.certified_centre_expansion(256, K=30) "
                                  "(the certified (a*, mu*) balls)",
                    components="tc_tm.centre_components(ce, hi), y in [0, hi]",
                    pairs="tc_commutator.build_pairs(AtomKP(at, comps, sd_shift=0))",
                    raise_nb="tc_layer_centre.raise_nb (offset 8 -> 10, kappa^1..kappa^7)",
                    tm_domain=[0.0, p["hi"]]),
        identity=dict(ok=ctx.identity["ok"], raise_by=ctx.identity["raise_by"],
                      kappa=ctx.identity["kappa"]),
        code=s43.code_record({DRIVER: DRIVER_SHA256}), driver_sha256=DRIVER_SHA256,
        env=s43.env_info(), timing=timing,
        provenance="computed by problems/P4/scripts/s46_centre.py work")
    s43.write_json_atomic(part_path(p["id"], t), _clean(rec))
    return rec


def cmd_launch(args):
    tag = args.tag or time.strftime("s46-centre-%Y%m%d-%H%M%S")
    os.makedirs(LOGS, exist_ok=True)
    workers = []
    for k in range(args.workers):
        log = os.path.join(LOGS, "%s-w%d.log" % (tag, k))
        argv = [sys.executable, DRIVER, "work", "--tau-mins", args.tau_mins,
                "--max-evals", str(args.max_evals)]
        if args.pieces:
            argv += ["--pieces", args.pieces]
        if args.force:
            argv.append("--force")

        def pre():
            os.nice(args.nice)
            signal.signal(signal.SIGHUP, signal.SIG_IGN)
        with open(log, "ab") as fh:
            proc = subprocess.Popen(["nohup"] + argv, cwd=s43.ROOT, stdin=subprocess.DEVNULL,
                                    stdout=fh, stderr=fh, start_new_session=True,
                                    preexec_fn=pre, env=dict(os.environ, PYTHONUNBUFFERED="1"))
        workers.append(dict(pid=proc.pid, log=s43.rel(log), argv=argv))
    rec = dict(tag=tag, host=socket.gethostname(), started=s43.now(), nice=args.nice,
               workers=workers, git_head=s43.git_head()[0], driver_sha256=DRIVER_SHA256,
               note="resumable: work skips parts valid for the current code and driver; "
                    "O_EXCL locks per piece in logs/s46/locks")
    s43.write_json_atomic(os.path.join(LOGS, "launch-%s.json" % tag), rec)
    for w in workers:
        print("pid %d  log %s" % (w["pid"], w["log"]))


def cmd_status(args):
    plan = s43.load_plan()
    rows = []
    for t in (224.0, 14.0):
        tot = dict(pass_=0, fail=0, error=0, missing=0, invalid=0)
        cpu = 0.0
        worst = None
        for p in centre_pieces(plan):
            path = part_path(p["id"], t)
            if not os.path.exists(path):
                tot["missing"] += 1
                continue
            d = json.load(open(path))
            if not part_valid(d, p, t):
                tot["invalid"] += 1
                continue
            tot["pass_" if d["status"] == "pass" else d["status"]] += 1
            cpu += (d.get("timing") or {}).get("certify_cpu_s") or 0.0
            b = d.get("certified_bound")
            if b is not None:
                worst = b if worst is None else max(worst, b)
            rows.append((t, p["id"], d["status"], b, (d.get("result") or {}).get("boxes"),
                         (d.get("result") or {}).get("evals"),
                         (d.get("timing") or {}).get("certify_cpu_s")))
        print("tau >= %g: pass %d fail %d error %d missing %d invalid %d; worst bound %s; "
              "certify cpu %.0f s" % (t, tot["pass_"], tot["fail"], tot["error"],
                                      tot["missing"], tot["invalid"], worst, cpu))
    if args.verbose:
        for r in rows:
            print("  t%-3d %-11s %-5s bound=%s boxes=%s evals=%s cpu=%s" % r)


# ------------------------------------------------------------------ grid
GRID_EPS = (0.0, 0.03, 0.2, 1.0, 2.0, 3.4, 4.6, 6.0)


def grid_points(p, tau_min):
    lo, hi = p["lo"], p["hi"]
    if lo > 0:
        ys = [lo, (lo * hi) ** 0.5, hi]
    else:
        ys = [hi, hi / 10, hi / 1e3, hi / 1e6, hi / 1e9]
    taus = [tau_min * 2.0 ** k for k in range(0, 40, 3)]
    return ys, taus


def cmd_grid(args):
    ctx = Ctx()
    tcc = ctx.tcc
    from p4.validated.arbseries import precision
    plan = s43.load_plan()
    out = dict(schema="s46-centre-grid/2", created=s43.now(),
               quantity="Sigma_c = S/(1 + s^2), s = y tau (sigma_box_centre at points, "
                        "raised layer data)",
               note=("'ball': the centre sigma at the point (certified-ball Taylor models); "
                     "'direct': S/(1 + s^2) by the direct chain (tc_layer.point_S_direct) "
                     "from point pairs, where its ball is tight"),
               pieces={})
    for p in centre_pieces(plan, args.pieces):
        ld, tb, cb = ctx.build(p)
        rows = []
        with precision(320):
            for tau_min in (224.0, 14.0):
                ys, taus = grid_points(p, tau_min)
                for y in ys:
                    pr = tcc.point_pairs_centre(ctx.at, ctx.ce, y) if y >= 1e-6 else None
                    for e in GRID_EPS:
                        for tau in taus:
                            r = dict(tau_mode=tau_min, y=y, eps=e, tau=tau, s=y * tau)
                            try:
                                v = tcc.sigma_box_centre(ld, tcc.point_box(y, e, tau))
                                r["ball"], r["ball_rad"] = float(v.mid()), float(v.rad())
                            except (ArithmeticError, ZeroDivisionError, ValueError) as ex:
                                r["ball"], r["ball_err"] = None, str(ex)[:80]
                            if pr is not None:
                                d, S = tcc.point_sigma_direct(pr, y, e, tau)
                                if float(d.rad()) < 1e-3 * abs(float(d.mid())):
                                    r["direct"] = float(d.mid())
                            rows.append(r)
        summ = {}
        for tau_min in (224.0, 14.0):
            rr = [r for r in rows if r["tau_mode"] == tau_min and r.get("ball") is not None]
            w = max(rr, key=lambda r: r["ball"]) if rr else None
            summ["tau>=%g" % tau_min] = dict(
                points=sum(1 for r in rows if r["tau_mode"] == tau_min),
                max_sigma_c=w and w["ball"],
                at=w and dict(y=w["y"], eps=w["eps"], tau=w["tau"], s=w["s"]),
                undefined=sum(1 for r in rows if r["tau_mode"] == tau_min
                              and r.get("ball") is None),
                max_abs_ball_minus_direct=max((abs(r["ball"] - r["direct"]) for r in rr
                                              if "direct" in r), default=None))
        out["pieces"][p["id"]] = dict(lo=p["lo"], hi=p["hi"], summary=summ, rows=rows)
        s43.log("grid %s: %s" % (p["id"], json.dumps(summ)))
    out["code"] = s43.code_record({DRIVER: DRIVER_SHA256})
    path = args.out or GRID
    s43.write_json_atomic(path, _clean(out))
    print("grid -> %s" % s43.rel(path))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    q = sub.add_parser("identity")
    q.add_argument("--out")
    q.set_defaults(fn=cmd_identity)
    for name, fn in (("work", cmd_work), ("launch", cmd_launch)):
        q = sub.add_parser(name)
        q.add_argument("--pieces", default=None)
        q.add_argument("--tau-mins", default="224,14")
        q.add_argument("--max-evals", type=int, default=150000)
        q.add_argument("--force", action="store_true")
        if name == "launch":
            q.add_argument("--workers", type=int, required=True)
            q.add_argument("--nice", type=int, default=10)
            q.add_argument("--tag")
        q.set_defaults(fn=fn)
    q = sub.add_parser("status")
    q.add_argument("--verbose", "-v", action="store_true")
    q.set_defaults(fn=cmd_status)
    q = sub.add_parser("grid")
    q.add_argument("--pieces", default="c29,c00")
    q.add_argument("--out")
    q.set_defaults(fn=cmd_grid)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
