#!/usr/bin/env python3
"""S4-6 probe: the centre variant of the layer engine on centre pieces.

The centre sigma-function is p4.validated.tc_layer_centre.sigma_box_centre
(graded pivot forms in y with z = y kappa, w = 1/kappa; Sigma_c = S/(1 + s^2),
s = y tau), plugged into the unchanged tc_layer.run_layer_piece.  Probe
stage: the natural enclosure only (no centred form), pieces c29 and c00,
tau >= 224 and tau >= 14.  See notes/s4-6-centre.md.

usage:
  s46_centre.py grid [--pieces c29,c00] [--out FILE]     float design grid
  s46_centre.py run [--pieces c29,c00] [--tau-mins 224,14] [--max-evals N]
                    [--ladder-evals N] [--force]          enclosure parts
  s46_centre.py launch [run options] [--nice 10] [--tag T] detached run
  s46_centre.py status                                      parts summary

Parts: problems/P4/results/theorem_c/centre_probe/t<tau_min>/<piece>.json,
with the import-time source hashes and the write-time check of
scripts/s43_parts.py (code_record).  Launch records and logs:
problems/P4/logs/s46/ (machine-local).
"""
from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import s43_parts as s43  # noqa: E402

DRIVER = os.path.abspath(__file__)
DRIVER_SHA256 = s43.sha256_file(DRIVER)
OUT_DIR = os.path.join(s43.TC_DIR, "centre_probe")
LOGS = os.path.join(s43.P4, "logs", "s46")
SCHEMA = "s46-centre-probe/1"
INF = float("inf")


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
        self.setup_s = time.time() - t0

    def piece(self, pid):
        plan = s43.load_plan()
        p = [q for q in plan["pieces"] if q["id"] == "centre/" + pid]
        if not p:
            raise SystemExit("no plan piece centre/%s" % pid)
        return p[0]

    def build(self, p):
        t0, c0 = time.time(), time.process_time()
        ld = self.tcc.centre_layer_data(self.at, self.ce, p["hi"], label=p["id"])
        return ld, time.time() - t0, time.process_time() - c0


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


def midpoint_layer_data(ld):
    """A copy of the layer data with every Taylor-model coefficient replaced
    by its midpoint (float design check only: the near-zero leading
    coefficients of nB, nB1, nB2 and nY become exact zeros)."""
    import copy
    m = copy.copy(ld)
    m.forms = {}
    for k, F in ld.forms.items():
        G = copy.copy(F)
        G.tm = {j: _mid_tm(t) for j, t in F.tm.items()}
        G._cache = {}
        m.forms[k] = G
    for nm in ("tt", "ttp", "u2"):
        lc = getattr(ld, nm)
        setattr(m, nm, type(lc)(_mid_tm(lc.tm), lc.off))
    m._rcache = {}
    return m


def _mid_tm(t):
    from flint import arb
    return type(t)([arb(c.mid()) for c in t.c], t.R, t.ve)


def cmd_grid(args):
    ctx = Ctx()
    tcc = ctx.tcc
    from p4.validated.arbseries import precision
    out = dict(schema="s46-centre-grid/1", created=s43.now(),
               quantity="Sigma_c = S/(1 + s^2), s = y tau (tc_layer_centre.sigma_graded at points)",
               note=("'ball': the graded centre sigma at the point with the Taylor models of the "
                     "certified (a*, mu*) balls (radius included; None where a pivot form's "
                     "leading coefficient is a ball containing 0); 'mid': the same with every "
                     "Taylor-model coefficient at its midpoint (a float design check).  "
                     "'direct': S/(1 + s^2) by the direct chain (tc_layer.point_S_direct) "
                     "from point pairs, where its ball is tight."),
               pieces={})
    for pid in args.pieces.split(","):
        p = ctx.piece(pid)
        ld, tb, cb = ctx.build(p)
        md = midpoint_layer_data(ld)
        rows = []
        with precision(320):
            for tau_min in (224.0, 14.0):
                ys, taus = grid_points(p, tau_min)
                for y in ys:
                    pr = tcc.point_pairs_centre(ctx.at, ctx.ce, y) if y >= 1e-6 else None
                    for e in GRID_EPS:
                        for tau in taus:
                            box = tcc.point_box(y, e, tau)
                            r = dict(tau_mode=tau_min, y=y, eps=e, tau=tau, s=y * tau)
                            for nm, L in (("ball", ld), ("mid", md)):
                                try:
                                    v = tcc.sigma_box_centre(L, box)
                                    r[nm] = float(v.mid())
                                    r[nm + "_rad"] = float(v.rad())
                                except (ArithmeticError, ZeroDivisionError, ValueError) as ex:
                                    r[nm] = None
                                    r[nm + "_err"] = str(ex)[:80]
                            if pr is not None:
                                d, S = tcc.point_sigma_direct(pr, y, e, tau)
                                if float(d.rad()) < 1e-3 * abs(float(d.mid())):
                                    r["direct"] = float(d.mid())
                            rows.append(r)
        ok = [r for r in rows if r.get("mid") is not None]
        worst = max(ok, key=lambda r: r["mid"]) if ok else None
        summ = {}
        for tau_min in (224.0, 14.0):
            rr = [r for r in ok if r["tau_mode"] == tau_min]
            w = max(rr, key=lambda r: r["mid"]) if rr else None
            near0 = [r for r in rr if r["eps"] <= 0.03]
            summ["tau>=%g" % tau_min] = dict(
                points=len(rr), max_sigma_c_mid=w and w["mid"],
                at=w and dict(y=w["y"], eps=w["eps"], tau=w["tau"], s=w["s"]),
                max_sigma_c_mid_eps_le_0_03=max((r["mid"] for r in near0), default=None),
                ball_undefined=sum(1 for r in rows if r["tau_mode"] == tau_min
                                   and r.get("ball") is None),
                max_abs_ball_minus_direct=max((abs(r["ball"] - r["direct"]) for r in rr
                                              if r.get("ball") is not None and "direct" in r),
                                             default=None))
        out["pieces"][pid] = dict(id=p["id"], lo=p["lo"], hi=p["hi"], build_wall_s=tb,
                                  summary=summ, worst=worst, rows=rows)
        s43.log("grid %s: %s" % (pid, json.dumps(summ)))
    out["code"] = s43.code_record({DRIVER: DRIVER_SHA256})
    path = args.out or os.path.join(OUT_DIR, "grid.json")
    s43.write_json_atomic(path, _clean(out))
    print("grid -> %s" % s43.rel(path))


# ------------------------------------------------------------------- run
def part_path(pid, tau_min):
    return os.path.join(OUT_DIR, "t%d" % int(tau_min), pid + ".json")


def cmd_run(args):
    ctx = None
    taus = [float(t) for t in args.tau_mins.split(",")]
    for pid in args.pieces.split(","):
        todo = [t for t in taus if args.force or not _done(pid, t)]
        if not todo:
            s43.log("skip %s (parts present for the current code)" % pid)
            continue
        if ctx is None:
            ctx = Ctx()
        p = ctx.piece(pid)
        start = s43.now()
        ld, tb, cb = ctx.build(p)
        s43.log("built %s in %.1f s" % (pid, tb))
        for t in todo:
            t1, c1 = time.time(), time.process_time()
            res, err = None, None
            try:
                from p4.validated.arbseries import precision
                with precision(320):
                    res = ctx.tl.run_layer_piece(
                        ld, p["lo"], p["hi"], tau_min=t, sigma_fn=ctx.tcc.sigma_box_centre,
                        ladder=(-1.0,), ladder_evals=args.ladder_evals,
                        max_evals=args.max_evals)
            except Exception:
                import traceback
                err = traceback.format_exc()
            timing = dict(start=start, end=s43.now(), build_wall_s=tb, build_cpu_s=cb,
                          certify_wall_s=time.time() - t1, certify_cpu_s=time.process_time() - c1)
            rec = write_part(p, t, res, err, timing, args)
            s43.log("done %s tau>=%g status=%s bound=%s boxes=%s evals=%s %.1fs"
                    % (pid, t, rec["status"], rec["certified_bound"],
                       (res or {}).get("boxes"), (res or {}).get("evals"),
                       timing["certify_wall_s"]))


def _done(pid, t):
    path = part_path(pid, t)
    if not os.path.exists(path):
        return False
    try:
        d = json.load(open(path))
    except Exception:
        return False
    return (d.get("driver_sha256") == DRIVER_SHA256 and s43.write_check_ok(d)
            and s43.code_matches(d["code"]["files"]))


def write_part(p, t, res, err, timing, args):
    status = "error" if err else ("pass" if res["ok"] else "fail")
    rec = dict(
        schema=SCHEMA, id=p["id"], spec=s43.spec_of(p), status=status,
        statement=("for all y in [%r, %r] (y > 0), eps in [0, 6], tau >= %g: dT, dB, nB != 0 "
                   "and Sigma_c = S/(1 + (y tau)^2) <= certified_bound < 0 "
                   "(tc_layer_centre.sigma_box_centre, natural enclosure)" % (p["lo"], p["hi"], t)),
        tau_min=t, certified_bound=(res or {}).get("certified_bound"),
        result=_clean({k: v for k, v in (res or {}).items() if k != "forms"}),
        forms=_clean((res or {}).get("forms")), error=err,
        inputs=dict(centre_series="ec_centre.certified_centre_expansion(256, K=30) "
                                  "(certified (a*, mu*) balls)",
                    components="tc_tm.centre_components(ce, hi), y in [0, hi]",
                    pairs="tc_commutator.build_pairs(AtomKP(at, comps, sd_shift=0))",
                    tm_domain=[0.0, p["hi"]]),
        engine=dict(sigma_fn="tc_layer_centre.sigma_box_centre", cover="tc_layer.run_layer_piece",
                    ladder=[-1.0], ladder_evals=args.ladder_evals, max_evals=args.max_evals,
                    kmax=None, s_split=None),
        code=s43.code_record({DRIVER: DRIVER_SHA256}), driver_sha256=DRIVER_SHA256,
        env=s43.env_info(), timing=timing,
        provenance="computed by problems/P4/scripts/s46_centre.py run")
    try:
        import p4.validated.tc_layer_centre as tcc
        rec["engine"].update(kmax=tcc.KMAX, s_split=tcc.S_SPLIT)
    except Exception:
        pass
    s43.write_json_atomic(part_path(p["id"].split("/")[1], t), _clean(rec))
    return rec


def cmd_launch(args):
    tag = args.tag or time.strftime("s46-centre-%Y%m%d-%H%M%S")
    os.makedirs(LOGS, exist_ok=True)
    log = os.path.join(LOGS, tag + ".log")
    argv = [sys.executable, DRIVER, "run", "--pieces", args.pieces, "--tau-mins", args.tau_mins,
            "--max-evals", str(args.max_evals), "--ladder-evals", str(args.ladder_evals)]
    if args.force:
        argv.append("--force")

    def pre():
        os.nice(args.nice)
        signal.signal(signal.SIGHUP, signal.SIG_IGN)
    with open(log, "ab") as fh:
        proc = subprocess.Popen(["nohup"] + argv, cwd=s43.ROOT, stdin=subprocess.DEVNULL,
                                stdout=fh, stderr=fh, start_new_session=True, preexec_fn=pre,
                                env=dict(os.environ, PYTHONUNBUFFERED="1"))
    rec = dict(tag=tag, pid=proc.pid, started=s43.now(), nice=args.nice, argv=argv,
               log=s43.rel(log), git_head=s43.git_head()[0], driver_sha256=DRIVER_SHA256,
               note="resumable: run skips parts present for the current code and driver")
    s43.write_json_atomic(os.path.join(LOGS, "launch-%s.json" % tag), rec)
    print("pid %d  log %s" % (proc.pid, s43.rel(log)))


def cmd_status(args):
    for t in (224, 14):
        d = os.path.join(OUT_DIR, "t%d" % t)
        if not os.path.isdir(d):
            continue
        for f in sorted(os.listdir(d)):
            r = json.load(open(os.path.join(d, f)))
            res = r.get("result") or {}
            print("t%-3d %-4s %-5s bound=%s boxes=%s evals=%s fails=%d valid=%s"
                  % (t, f[:-5], r["status"], r.get("certified_bound"), res.get("boxes"),
                     res.get("evals"), len(res.get("fails") or []),
                     s43.write_check_ok(r) and s43.code_matches(r["code"]["files"])))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    q = sub.add_parser("grid")
    q.add_argument("--pieces", default="c29,c00")
    q.add_argument("--out")
    q.set_defaults(fn=cmd_grid)
    for name, fn in (("run", cmd_run), ("launch", cmd_launch)):
        q = sub.add_parser(name)
        q.add_argument("--pieces", default="c29,c00")
        q.add_argument("--tau-mins", default="224,14")
        q.add_argument("--max-evals", type=int, default=20000)
        q.add_argument("--ladder-evals", type=int, default=5000)
        q.add_argument("--force", action="store_true")
        if name == "launch":
            q.add_argument("--nice", type=int, default=10)
            q.add_argument("--tag")
        q.set_defaults(fn=fn)
    q = sub.add_parser("status")
    q.set_defaults(fn=cmd_status)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
