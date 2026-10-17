#!/usr/bin/env python3
"""S4-3 region-part driver for the commutator strip certificate (tc_commutator).

Persistent replacement for the September scratchpad drivers (worker.py,
gap_worker.py, gapd_worker.py, sonic_worker2.py, run_centre*.py, tchunk.py,
dchunk.py, tband2/3.py, tubeband.py, build_*_tube.py, curate_merge.py and the
cloud's cworker.py), all of which were lost from /tmp in the 2026-09-07 reboot.
Nothing in p4.validated is modified: every piece is evaluated by the committed
region functions (run_tube / run_sonic with a one-piece slice, and run_centre's
per-piece body centre_piece), so a part computed here is a replay of the code
whose source hashes it records.

Layout (all paths relative to the repository root, so the script runs
unchanged on the Mac and on the cloud node):

  problems/P4/data/cache/tube_<name>.json            tube caches (gitignored; pinned
                                                     by content hash in parts/caches/)
  problems/P4/results/theorem_c/parts/caches/<name>.json  cache build records (committed)
  problems/P4/results/theorem_c/parts/plan.json      every piece of the certificate (committed)
  problems/P4/results/theorem_c/parts/<region>/<id>.json  one result per piece (committed)
  problems/P4/results/theorem_c/parts/checks.json    exact identities E1/E2/E3/E5 (committed)
  problems/P4/results/theorem_c/parts/lemma_d*.json  Lemma D re-runs (committed)
  problems/P4/results/theorem_c/tc_commutator.json   merged certificate (only when complete)
  problems/P4/logs/s43/                              worker logs, launch records, locks
                                                     (gitignored)

Commands (python problems/P4/scripts/s43_parts.py <cmd> -h for options):

  build-caches [--only full,gap,deep,band] [--detach] [--nice 10]
  cache-hash FILE...                 content hash of a tube cache (steps + K, not info)
  plan                               enumerate all pieces -> parts/plan.json
                                     (per piece: tq = False for the direct-only
                                     pieces, tstar_fac for the centre variant)
  checks                             E1/E2/E3/E5 -> parts/checks.json
  lemma-d [--centre-balls certified|midpoint] [--replay]
                                     Lemma D -> parts/lemma_d_<balls>.json
                                     (--replay: parts/replays/lemma_d_<balls>.json,
                                     compared with the stored record)
  p43-centre                         P4-3 centre tail (tc_monotone.run_centre) with
                                     the certified balls, and the midpoint-ball
                                     replay of the old record -> parts/p43_centre.json
  bench --ids ID,.. [--modes tq,direct] [--detach --nice 10 --tag T]
                                     time analyze_piece with tq=True / tq=False on
                                     single pieces (phase split recorded) ->
                                     results/theorem_c/bench/<id>__<mode>.json
                                     (resumable: finished modes are skipped)
  list [selection options]           print the selected piece ids (machine splits)
  work [selection options] [--order cost|plan|deep-first]
                                     foreground worker loop (skips valid results)
  launch --workers N [--nice 10] [selection options]
                                     detached workers (new session, nohup-safe),
                                     sharing one lock-guarded queue
  selection options: --regions R,.. --ids PAT,.. --ids-file F --exclude-ids-file F
                     --x-range LO,HI --shard K/M
  status [-v]                        progress per region
  import-cloud DIR --code-manifest F [--tube-hash name=sha ...]
                                     Branch A: convert the node's cw_*.json chunk
                                     reports into per-piece parts
  cloud-manifest --out F             Branch A: {source: sha256} the cworkers executed
  replay ID                          recompute one piece, compare with its part
                                     (record in parts/replays/, part untouched)
  merge [--accept-code FILE=SHA ...] [--accept-reason TEXT]
                                     coverage-checked assembly; writes
                                     tc_commutator.json only when complete,
                                     tc_commutator_partial.json otherwise

A piece's result is "valid" (skipped on re-run) iff its file parses, its spec
matches the plan entry, its tube cache content hash matches the local cache,
and every source file it recorded still has the recorded sha256.  Pieces that
ended with a certified failure are also skipped (the computation is
deterministic) unless --retry-failed is given.  Results are written
atomically; a worker claims a piece with an O_EXCL lock file, so any number of
workers on one machine can share the queue.  Different machines must work on
disjoint selections (--regions / --ids / --shard) and exchange part files.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import fnmatch
import glob
import hashlib
import json
import math
import os
import platform
import signal
import socket
import subprocess
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
P4 = os.path.join(ROOT, "problems", "P4")
SRC = os.path.join(P4, "src")
CACHE_DIR = os.path.join(P4, "data", "cache")
TC_DIR = os.path.join(P4, "results", "theorem_c")
PARTS = os.path.join(TC_DIR, "parts")
LOGS = os.path.join(P4, "logs", "s43")
LOCKS = os.path.join(LOGS, "locks")
PLAN = os.path.join(PARTS, "plan.json")
CERT = os.path.join(TC_DIR, "tc_commutator.json")
CERT_PARTIAL = os.path.join(TC_DIR, "tc_commutator_partial.json")
SCHEMA = "s43-part/1"

if SRC not in sys.path:
    sys.path.insert(0, SRC)

B1 = 14.0

# ----------------------------------------------------------------- utilities


def rel(path):
    return os.path.relpath(os.path.abspath(path), ROOT)


def now():
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def tube_content_hash(path):
    """sha256 over the canonical JSON of (K, steps) of a Tube.save file.  The
    info block is excluded (it carries the wall-clock build time), so two
    builds of the same tube agree.  Needs only json + hashlib, so the same
    function can be pasted into a shell on the cloud node."""
    d = json.load(open(path))
    payload = json.dumps(dict(K=d["K"], steps=d["steps"]), sort_keys=True,
                         separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def write_json_atomic(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = "%s.tmp.%d" % (path, os.getpid())
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, default=str)
        fh.write("\n")
    os.replace(tmp, path)


def log(msg):
    print("[%s %s:%d] %s" % (now(), socket.gethostname().split(".")[0],
                             os.getpid(), msg), flush=True)


def git_head():
    try:
        h = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"],
                           capture_output=True, text=True, timeout=20).stdout.strip()
        dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain",
                                "--", "problems/P4/src"], capture_output=True,
                               text=True, timeout=20).stdout.strip()
        return h or None, bool(dirty)
    except Exception:          # no git on the node checkout: provenance only
        return None, None


def env_info():
    try:
        import flint
        fv = flint.__version__
    except Exception:
        fv = None
    return dict(host=socket.gethostname(), platform=platform.platform(),
                machine=platform.machine(), cpu_count=os.cpu_count(),
                python=sys.version.split()[0], flint=fv)


# ------------------------------------------------------------ code identity
EAGER = ("tc_commutator", "tc_tm", "tc_energy", "tc_reduce", "lintube", "centre",
         "sonic", "shootsys", "tailbound", "matching", "recursion", "tmint",
         "linstep", "polysys", "systems", "linsys", "arbseries")
_HASH_CACHE = {}


def import_p4():
    import importlib
    for m in EAGER:
        importlib.import_module("p4.validated." + m)


def loaded_code():
    """{repo-relative path: sha256} of every p4.* module loaded in this process."""
    out = {}
    for name, mod in sorted(sys.modules.items()):
        if not (name == "p4" or name.startswith("p4.")):
            continue
        f = getattr(mod, "__file__", None)
        if not f or not f.endswith(".py"):
            continue
        r = rel(f)
        if r not in _HASH_CACHE:
            _HASH_CACHE[r] = sha256_file(f)
        out[r] = _HASH_CACHE[r]
    return out


def fingerprint(files):
    h = hashlib.sha256()
    for k in sorted(files):
        h.update(("%s=%s\n" % (k, files[k])).encode())
    return h.hexdigest()


def current_hash(relpath):
    p = os.path.join(ROOT, relpath)
    if not os.path.exists(p):
        return None
    if relpath not in _HASH_CACHE:
        _HASH_CACHE[relpath] = sha256_file(p)
    return _HASH_CACHE[relpath]


def code_matches(files, accept=None):
    """True iff every recorded source file still has its recorded hash (or a
    hash explicitly accepted via accept = {relpath: {sha, ...}})."""
    accept = accept or {}
    for r, sha in files.items():
        cur = current_hash(r)
        if sha == cur:
            continue
        if sha in accept.get(r, ()) or sha in accept.get(os.path.basename(r), ()):
            continue
        return False
    return True


# ----------------------------------------------------------------- caches
V0C, V0W = "0.1124394013880983", 1e-16     # = tc_energy.V0C/V0W = modecount.V0_EC/W_V0


def certified_centre_balls():
    """(a, mu) balls of the certified EC solution: a3_midpoint.json midpoint +
    its stored Krawczyk offset K(X) - m.  Equal to the constants used by
    modecount.certified_centre and by the September tube builds.  Since
    2026-10-16 the construction lives in p4.validated.ec_centre (shared with
    tc_commutator._centre_expansion and tc_monotone.run_centre); at the
    working precisions 256/320/384 it returns the same balls as the driver's
    original construction (float radius; the endpoint guard does not fire)."""
    from p4.validated.ec_centre import certified_centre_balls as ccb
    return ccb()


def midpoint_centre_balls():
    """(a, mu) exactly as tc_commutator._centre_expansion and
    tc_monotone.run_centre built them before 2026-10-16 (midpoint only, radius
    <= 1e-41; does NOT contain the certified (a*, mu*)).  Kept only to replay
    the old records."""
    from flint import arb
    aj = json.load(open(os.path.join(P4, "results", "a3_midpoint.json")))
    ball = lambda s: arb(s[1:-1].split(" +/-")[0]) + arb(0, float(s[1:-1].split(" +/- ")[1]))
    return ball(aj["a"]), ball(aj["mu"])


def centre_expansion(balls="certified", prec=256):
    from p4.validated import centre
    from p4.validated.arbseries import precision
    if balls == "certified":
        from p4.validated.ec_centre import certified_centre_expansion
        return certified_centre_expansion(prec, K=30)
    with precision(prec):
        a, mu = midpoint_centre_balls()
        ce = centre.centre_expansion(mu * (2 * a).exp(), nhat=(-a).exp(), K=30)
        if not ce.certify().ok:
            raise RuntimeError("centre certificate failed (%s balls)" % balls)
    return ce


CACHES = {
    # name: (description, recipe) -- the exact September constructions
    "full": ("main tube: A3 tube [-3, -0.05] (244 steps, hmax 0.02) + centre "
             "extension to -8 at h = 0.05 (100 steps); = scratchpad tube_full.json "
             "(built 2026-08-31, used by S4-2 and S4-3)",
             "Tube.build(V0C, 1e-16, x_c=-3.0); with precision(256): "
             "ce = certified centre; tube.extend_centre(ce, -8.0, hmax=0.05)"),
    "gap": ("gap tube [-0.055, -0.0085] (107 steps)",
            "Tube.build(V0C, 1e-16, x0=-0.0085, x_c=-0.055, K=28, hmax=0.004, prec=384)"),
    "deep": ("deep re-step [-4.56, -2.95] at h = 0.0125 (129 steps; replaces the "
             "h = 0.05 steps at x <= -3.05 whose chain TM inverses fail)",
             "with precision(384): ce = certified centre; Tube([], 28, "
             "info=dict(x_c=-2.95)).extend_centre(ce, -4.56, hmax=0.0125, prec=384)"),
    "band": ("band re-step: A3 tube [-0.17, -0.05] at hmax 0.002 (replaces the 4 "
             "h = 0.0075 steps in [-0.169, -0.139] whose chain TM inverses fail)",
             "Tube.build(V0C, 1e-16, x0=-0.05, x_c=-0.17, K=28, hmax=0.002, prec=384)"),
}


def cache_path(name):
    return os.path.join(CACHE_DIR, "tube_%s.json" % name)


def cache_record_path(name):
    return os.path.join(PARTS, "caches", "%s.json" % name)


def build_one_cache(name):
    from flint import arb  # noqa: F401
    from p4.validated import lintube
    from p4.validated.arbseries import precision
    import_p4()
    t0, c0 = time.time(), time.process_time()
    log("building cache %s: %s" % (name, CACHES[name][1]))
    if name == "full":
        tube = lintube.Tube.build(V0C, V0W, x_c=-3.0)
        with precision(256):
            ce = centre_expansion("certified", prec=256)
            tube.extend_centre(ce, -8.0, hmax=0.05)
    elif name == "gap":
        tube = lintube.Tube.build(V0C, V0W, x0=-0.0085, x_c=-0.055, K=28,
                                  hmax=0.004, prec=384)
    elif name == "deep":
        with precision(384):
            ce = centre_expansion("certified", prec=384)
            tube = lintube.Tube([], 28, info=dict(x_c=-2.95))
            tube.extend_centre(ce, -4.56, hmax=0.0125, prec=384)
    elif name == "band":
        tube = lintube.Tube.build(V0C, V0W, x0=-0.05, x_c=-0.17, K=28,
                                  hmax=0.002, prec=384)
    else:
        raise SystemExit("unknown cache " + name)
    os.makedirs(CACHE_DIR, exist_ok=True)
    out = cache_path(name)
    tmp = out + ".tmp.%d" % os.getpid()
    tube.save(tmp)
    os.replace(tmp, out)
    xs = [float(s.x) for s in tube.steps]
    rec = dict(schema="s43-cache/1", name=name, description=CACHES[name][0],
               recipe=CACHES[name][1], file=rel(out),
               content_sha256=tube_content_hash(out), file_sha256=sha256_file(out),
               steps=len(tube.steps), x_first=xs[0],
               x_last_end=float(tube.steps[-1].x) - tube.steps[-1].h,
               h_min=min(s.h for s in tube.steps), h_max=max(s.h for s in tube.steps),
               info={k: str(v) for k, v in tube.info.items()},
               build_wall_s=time.time() - t0, build_cpu_s=time.process_time() - c0,
               built=now(), env=env_info(), code=loaded_code(),
               driver_sha256=sha256_file(os.path.abspath(__file__)))
    write_json_atomic(cache_record_path(name), rec)
    log("cache %s: %d steps, %.0f s, content %s" % (name, rec["steps"],
                                                    rec["build_wall_s"],
                                                    rec["content_sha256"][:16]))
    return rec


def cache_ok(name):
    """The local cache exists and matches its committed build record."""
    p, r = cache_path(name), cache_record_path(name)
    if not (os.path.exists(p) and os.path.exists(r)):
        return False
    rec = json.load(open(r))
    return tube_content_hash(p) == rec["content_sha256"]


# ---------------------------------------------------------------- regions
# The September region parameters (transcript of agent a7accaae8fcc3cf60):
#   tube: worker.py / tchunk.py / cworker.py (degrade=True from 09-01 18:12 on)
#   band: tubeband.py (x in [-0.1665, -0.134]) + tband3.py ([-0.1705, -0.1655])
#   deep: cworker.py chunks C0-C5 (x_min=-4.6, x_max=-2.9)
#   gap:  gap_worker.py / gapd_worker.py (x_min=-0.056, x_max=-0.008)
#   sonic: sonic_worker2.py (run_sonic, 51 pieces below XI_SONIC_TOP)
#   centre: run_centre_final.py (run_centre, 30 pieces)
# Per-piece phase settings (recorded in every plan entry and part spec):
#  * tq = False ("direct_only"): the tau-quadratic phase is not attempted; the
#    piece is at best pass-direct (tau in [b1, tstar] certified, tau > tstar
#    open).  Used where the September tau-quadratic phase failed on every
#    piece: gap idx 0-88 (xi in [0.0085, 0.0357]) and all 51 sonic pieces;
#    the S4-4 layer engine is to carry tau > tstar there (recovery note 8(i)).
#  * tstar_fac (centre): the tau-quadratic phase starts at tstar = tstar_fac*b1;
#    tstar_fac = 1 makes it cover every tau >= b1 with no direct boxes
#    (recovery note 8(iii); the 04:18Z September code closed the centre so).
REGIONS = {
    "centre": dict(kind="centre", y0=1e-8, ratio=1.6, prec=320, ncells=6,
                   centre_balls="certified", tstar_fac=1.0),
    "tube": dict(kind="tube", cache="full", x_min=-4.5, x_max=-0.05, degrade=True,
                 exclude=dict(band=(-0.1705, -0.134), deep_below=-3.05)),
    "band": dict(kind="tube", cache="band", x_min=-0.1705, x_max=-0.134, degrade=True),
    "deep": dict(kind="tube", cache="deep", x_min=-4.6, x_max=-2.9, degrade=True),
    "gap": dict(kind="tube", cache="gap", x_min=-0.056, x_max=-0.008, degrade=True,
                direct_only=dict(index_max=88, why=(
                    "xi <= 0.0357: the tau-quadratic phase failed on every one of "
                    "idx 0-88 in September (bounds +8e8 ... +7.8e10); direct "
                    "phase only, tau > tstar left to the S4-4 layer engine"))),
    "sonic": dict(kind="sonic", xi0=1e-6, prec=320,
                  direct_only=dict(index_max=50, why=(
                      "all 51 pieces passed only the direct phase in September; "
                      "direct phase only, tau > tstar left to the S4-4 layer "
                      "engine"))),
}
REGION_ORDER = ("centre", "tube", "band", "deep", "gap", "sonic")


def piece_tq(region, index):
    do = REGIONS[region].get("direct_only")
    return not (do is not None and index <= do["index_max"])


def est_core_hours(region, x=None, idx=None, tq=True, tstar_fac=None):
    """Rough Mac-performance-core hours per piece from the September logs
    (see notes/s4-3-recovery.md, inventory).  Used for ordering and ETAs only."""
    if region == "centre":
        if tstar_fac is not None and tstar_fac <= 1:
            return 0.05     # tau-quadratic from b1, no direct boxes (04:18Z code:
                            # 30 pieces in 1022 s)
        return 2.0          # >= 0.85 (25.3 CPU-h for < 30 pieces, never finished); the
                            # 2-cell unit-test piece did not finish in 15 min (2026-10-16)
    if region == "sonic":
        return 0.08
    if region == "band":
        return 0.03
    if region == "gap":
        return 0.2 if tq else 0.06   # gapd: direct phase ~19179 s / 90 steps
    if region == "deep":
        return 2.0          # cloud C chunks (21-22 steps) unfinished at 65 h
    if region == "tube":
        if x is None:
            return 0.1
        if x >= -0.2:
            return 0.12
        if x >= -0.8:
            return 0.04
        if x >= -1.52:
            return 0.03
        if x >= -1.82:
            return 1.8      # cloud A1: 29.3 vCPU-h / 8 steps, vCPU taken as 0.5 Mac core
        return 4.0          # cloud B chunks (8 steps) unfinished at 65 h: >= 8 vCPU-h/step
    return 0.1


def build_plan():
    import p4.validated.tc_commutator as tc
    from p4.validated.lintube import Tube
    pieces, excluded = [], []
    caches = {}
    for name in CACHES:
        if not cache_ok(name):
            raise SystemExit("cache %s missing or not matching its record; run "
                             "build-caches first" % name)
        caches[name] = json.load(open(cache_record_path(name)))
    tubes = {}
    for region in REGION_ORDER:
        rd = REGIONS[region]
        if rd["kind"] == "centre":
            tf = rd.get("tstar_fac")
            for i, (lo, hi) in enumerate(tc.y_subranges(rd["y0"], rd["ratio"])):
                ent = dict(id="centre/c%02d" % i, region=region, kind="centre",
                           index=i, lo=lo, hi=hi,
                           x_cover=[None if lo <= 0 else math.log(lo), math.log(hi)],
                           est_core_h=est_core_hours(region, tstar_fac=tf))
                if tf is not None:
                    ent["tstar_fac"] = tf
                pieces.append(ent)
        elif rd["kind"] == "sonic":
            subs = [t for t in tc.xi_subranges(rd["xi0"]) if t[0] < tc.XI_SONIC_TOP]
            for i, (lo, hi) in enumerate(subs):
                tq = piece_tq(region, i)
                pieces.append(dict(id="sonic/s%02d" % i, region=region, kind="sonic",
                                   index=i, lo=lo, hi=hi, x_cover=[-hi, -lo], tq=tq,
                                   est_core_h=est_core_hours(region, tq=tq)))
        else:
            name = rd["cache"]
            if name not in tubes:
                tubes[name] = Tube.load(cache_path(name))
            steps = [sd for sd in tubes[name].steps
                     if float(sd.x) - sd.h >= rd["x_min"] - 1e-12
                     and float(sd.x) <= rd["x_max"] + 1e-12]
            for i, sd in enumerate(steps):
                x, h = float(sd.x), float(sd.h)
                tq = piece_tq(region, i)
                ent = dict(id="%s/%s-%03d" % (region, name, i), region=region,
                           kind="tube", cache=name,
                           cache_content_sha256=caches[name]["content_sha256"],
                           index=i, x=x, h=h, x_cover=[x - h, x],
                           x_min=rd["x_min"], x_max=rd["x_max"], degrade=rd["degrade"],
                           tq=tq, est_core_h=est_core_hours(region, x=x, idx=i, tq=tq))
                ex = rd.get("exclude")
                if ex:
                    blo, bhi = ex["band"]
                    if x - h >= blo - 1e-12 and x <= bhi + 1e-12:
                        excluded.append(dict(ent, reason="inside the band window "
                                             "[-0.1705, -0.134]: September TM-inverse "
                                             "failure at h = 0.0075, replaced by band/*"))
                        continue
                    if x <= ex["deep_below"] + 1e-9:
                        excluded.append(dict(ent, reason="x <= -3.05 (h = 0.05 centre "
                                             "extension): September TM-inverse failure, "
                                             "replaced by deep/*"))
                        continue
                pieces.append(ent)
    import_p4()
    plan = dict(schema="s43-plan/1", created=now(), b1=B1,
                tstar=tc.TSTAR_FAC * B1, regions=REGIONS, caches=caches,
                pieces=pieces, excluded=excluded,
                counts={r: sum(1 for p in pieces if p["region"] == r) for r in REGION_ORDER},
                code=loaded_code(), driver_sha256=sha256_file(os.path.abspath(__file__)),
                note=("Coverage requirement checked by merge: the union of the "
                      "certified pieces' x-intervals must contain (-inf, 0]. "
                      "centre: y = e^x in [lo, hi]; tube-like: [x - h, x]; sonic: "
                      "xi = -x in [lo, hi]."))
    return plan


def load_plan():
    if not os.path.exists(PLAN):
        raise SystemExit("no plan: run `plan` first")
    return json.load(open(PLAN))


def part_path(pid):
    return os.path.join(PARTS, pid + ".json")


# ------------------------------------------------------------ piece runner
class Ctx:
    def __init__(self, centre_balls="certified"):
        t0 = time.time()
        import_p4()
        import p4.validated.tc_commutator as tc
        from p4.validated.tc_energy import Atoms
        self.tc = tc
        self.at = Atoms(full=True)
        self.atoms_s = time.time() - t0
        self._tubes, self._ce = {}, None
        self.centre_balls = centre_balls
        log("context ready (atoms %.1f s)" % self.atoms_s)

    def tube(self, name):
        if name not in self._tubes:
            from p4.validated.lintube import Tube
            if not cache_ok(name):
                raise RuntimeError("cache %s missing or not matching its record" % name)
            self._tubes[name] = Tube.load(cache_path(name))
        return self._tubes[name]

    def centre_ce(self):
        if self._ce is None:
            if self.centre_balls == "certified":
                self._ce = self.tc._centre_expansion()   # the committed code path
            else:
                self._ce = centre_expansion(self.centre_balls)
        return self._ce


def classify(res):
    if not res["ok"]:
        return "error" if res.get("err") else "fail"
    return "pass-direct" if res.get("open") else "pass-both"


def run_piece(p, ctx):
    tc = ctx.tc
    tq = p.get("tq", True)
    if p["kind"] == "tube":
        tube = ctx.tube(p["cache"])
        rep = tc.run_tube(ctx.at, tube, B1, x_min=p["x_min"], x_max=p["x_max"],
                          step_slice=(p["index"], p["index"] + 1),
                          degrade=p["degrade"], tq=tq)
        if rep["pieces"] != 1:
            raise RuntimeError("slice returned %d pieces" % rep["pieces"])
        d = rep["detail"][0]
        if d["x"] != p["x"] or d["h"] != p["h"]:
            raise RuntimeError("step mismatch: plan (%r, %r) vs cache (%r, %r)"
                               % (p["x"], p["h"], d["x"], d["h"]))
        return d, dict(region_call="run_tube(at, tube_%s, %g, x_min=%g, x_max=%g, "
                       "step_slice=(%d, %d), degrade=%s, tq=%s)"
                       % (p["cache"], B1, p["x_min"], p["x_max"], p["index"],
                          p["index"] + 1, p["degrade"], tq), tq=tq)
    if p["kind"] == "sonic":
        rep = tc.run_sonic(ctx.at, B1, piece_slice=(p["index"], p["index"] + 1), tq=tq)
        if rep["pieces"] != 1:
            raise RuntimeError("slice returned %d pieces" % rep["pieces"])
        d = rep["detail"][0]
        if d["lo"] != p["lo"] or d["hi"] != p["hi"]:
            raise RuntimeError("sonic piece mismatch")
        return d, dict(region_call="run_sonic(at, %g, piece_slice=(%d, %d), tq=%s)"
                       % (B1, p["index"], p["index"] + 1, tq), tq=tq)
    if p["kind"] == "centre":
        # run_centre's per-piece body, verbatim: the committed
        # _centre_expansion (certified (a*, mu*) balls since 2026-10-16) and
        # centre_piece with the piece's tstar_fac (None = TSTAR_FAC)
        rd = REGIONS["centre"]
        ce = ctx.centre_ce()
        tf = p.get("tstar_fac")
        try:
            rep = tc.centre_piece(ctx.at, ce, p["lo"], p["hi"], B1,
                                  tc.eps_cells_default(rd["ncells"]), rd["prec"],
                                  tstar_fac=tf)
        except (ArithmeticError, RuntimeError) as e:
            rep = dict(ok=False, err=str(e), cells=[])
        mg = max((c.get("margin") for c in rep.get("cells", [])
                  if c.get("margin") is not None), default=None)
        d = dict(lo=p["lo"], hi=p["hi"], ok=rep["ok"], margin=mg, err=rep.get("err"),
                 errs=[c.get("err") for c in rep.get("cells", []) if c.get("err")][:2] or None)
        return d, dict(region_call="centre_piece(at, _centre_expansion() [%s balls], %r, "
                       "%r, %g, eps_cells_default(), %d, tstar_fac=%r)"
                       % (ctx.centre_balls, p["lo"], p["hi"], B1, rd["prec"], tf),
                       centre_balls=ctx.centre_balls, tstar_fac=tf,
                       cells=rep.get("cells"))
    raise RuntimeError("unknown kind " + p["kind"])


def piece_tstar(p, tc):
    """tstar of a piece: tstar_fac * b1 (default TSTAR_FAC)."""
    tf = p.get("tstar_fac")
    return (tc.TSTAR_FAC if tf is None else tf) * B1


SPEC_KEYS = ("id", "region", "kind", "index", "x", "h", "lo", "hi", "cache",
             "cache_content_sha256", "x_min", "x_max", "degrade", "tq", "tstar_fac")


def spec_of(p):
    return {k: p[k] for k in SPEC_KEYS if k in p}


def read_part(pid):
    path = part_path(pid)
    if not os.path.exists(path):
        return None
    try:
        d = json.load(open(path))
    except Exception:
        return None
    return d if d.get("schema") == SCHEMA else None


def part_state(p, d, accept=None):
    """missing | stale | pass-both | pass-direct | fail | error"""
    if d is None:
        return "missing"
    if d.get("spec") != spec_of(p):
        return "stale"
    if p["kind"] == "tube":
        rec = json.load(open(cache_record_path(p["cache"])))
        if rec["content_sha256"] != p["cache_content_sha256"]:
            return "stale"
    if p["kind"] == "centre" and d.get("extra", {}).get("centre_balls") != \
            REGIONS["centre"]["centre_balls"]:
        return "stale"
    if not code_matches(d["code"]["files"], accept):
        return "stale"
    return d["status"]


def claim(pid):
    os.makedirs(LOCKS, exist_ok=True)
    lk = os.path.join(LOCKS, pid.replace("/", "__") + ".lock")
    payload = json.dumps(dict(host=socket.gethostname(), pid=os.getpid(), t=now()))
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
                return None                       # live owner
            except ProcessLookupError:
                os.remove(lk)                     # stale lock from a dead worker
            except PermissionError:
                return None
    return None


def _read_ids(path):
    return {ln.strip() for ln in open(path) if ln.strip() and not ln.startswith("#")}


def select(plan, args):
    """Filters in this order: --regions, --ids, --ids-file, --exclude-ids-file,
    --x-range (pieces whose x-cover lies inside [lo, hi]), --shard k/m (every
    m-th of the remaining pieces in plan order), then --order."""
    ps = plan["pieces"]
    if getattr(args, "regions", None):
        want = set(args.regions.split(","))
        ps = [p for p in ps if p["region"] in want]
    if getattr(args, "ids", None):
        pats = args.ids.split(",")
        ps = [p for p in ps if any(fnmatch.fnmatch(p["id"], q) for q in pats)]
    if getattr(args, "ids_file", None):
        keep = _read_ids(args.ids_file)
        ps = [p for p in ps if p["id"] in keep]
    if getattr(args, "exclude_ids_file", None):
        drop = _read_ids(args.exclude_ids_file)
        ps = [p for p in ps if p["id"] not in drop]
    if getattr(args, "x_range", None):
        lo, hi = (float(v) for v in args.x_range.split(","))
        ps = [p for p in ps if x_interval(p)[0] >= lo - 1e-12 and x_interval(p)[1] <= hi + 1e-12]
    if getattr(args, "shard", None):
        k, m = (int(v) for v in args.shard.split("/"))
        ps = [p for i, p in enumerate(ps) if i % m == k]
    order = getattr(args, "order", "cost") or "cost"
    if order == "cost":
        ps = sorted(ps, key=lambda p: (p["est_core_h"], p["id"]))
    elif order == "deep-first":
        ps = sorted(ps, key=lambda p: (-p["est_core_h"], p["id"]))
    return ps


def cmd_work(args):
    plan = load_plan()
    pieces = select(plan, args)
    accept = parse_accept(args.accept_code)
    ctx = None
    done = 0
    log("worker start: %d candidate pieces (order=%s)" % (len(pieces), args.order))
    for p in pieces:
        if args.max_pieces and done >= args.max_pieces:
            break
        st = part_state(p, read_part(p["id"]), accept)
        if st in ("pass-both", "pass-direct") or (st in ("fail", "error")
                                                  and not args.retry_failed):
            continue
        lk = claim(p["id"])
        if lk is None:
            continue
        try:
            st = part_state(p, read_part(p["id"]), accept)     # re-check under lock
            if st in ("pass-both", "pass-direct") or (st in ("fail", "error")
                                                      and not args.retry_failed):
                continue
            if ctx is None:
                ctx = Ctx(REGIONS["centre"]["centre_balls"])
            log("start %s (est %.2f core-h)" % (p["id"], p["est_core_h"]))
            t0, c0, start = time.time(), time.process_time(), now()
            res, extra = run_piece(p, ctx)
            files, (head, dirty) = loaded_code(), git_head()
            rec = dict(schema=SCHEMA, id=p["id"], region=p["region"], spec=spec_of(p),
                       status=classify(res), result=res, b1=B1,
                       tstar=piece_tstar(p, ctx.tc), extra=extra,
                       code=dict(files=files, fingerprint=fingerprint(files),
                                 git_head=head, git_dirty_src=dirty),
                       driver_sha256=sha256_file(os.path.abspath(__file__)),
                       env=env_info(),
                       timing=dict(start=start, end=now(), wall_s=time.time() - t0,
                                   cpu_s=time.process_time() - c0,
                                   loadavg=list(os.getloadavg())),
                       provenance="computed by problems/P4/scripts/s43_parts.py work")
            write_json_atomic(part_path(p["id"]), rec)
            done += 1
            log("done  %s status=%s margin=%s wall=%.0fs cpu=%.0fs"
                % (p["id"], rec["status"], res.get("margin"), rec["timing"]["wall_s"],
                   rec["timing"]["cpu_s"]))
        except Exception:
            log("EXCEPTION on %s (no result written; will be retried):\n%s"
                % (p["id"], traceback.format_exc()))
        finally:
            try:
                os.remove(lk)
            except OSError:
                pass
    log("worker exit: %d pieces computed" % done)


# ------------------------------------------------------------- detaching
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
                            cwd=ROOT, stdin=subprocess.DEVNULL, stdout=out, stderr=out,
                            start_new_session=True, preexec_fn=pre,
                            env=dict(os.environ, PYTHONUNBUFFERED="1"))
    return proc.pid


def selection_argv(args):
    av = []
    for k in ("regions", "ids", "ids_file", "exclude_ids_file", "x_range", "shard",
              "order"):
        v = getattr(args, k, None)
        if v:
            if k in ("ids_file", "exclude_ids_file"):
                v = os.path.abspath(v)
            av += ["--" + k.replace("_", "-") + "=" + v]
    if getattr(args, "max_pieces", 0):
        av += ["--max-pieces", str(args.max_pieces)]
    if getattr(args, "retry_failed", False):
        av.append("--retry-failed")
    for a in getattr(args, "accept_code", None) or []:
        av += ["--accept-code", a]
    return av


def cmd_launch(args):
    tag = args.tag or time.strftime("%Y%m%d-%H%M%S")
    pids = []
    for k in range(args.workers):
        av = ["work"] + selection_argv(args)          # workers share one lock-guarded queue
        lf = os.path.join(LOGS, "%s-w%d.log" % (tag, k))
        pids.append((spawn_detached(av, lf, args.nice), rel(lf)))
    rec = dict(tag=tag, host=socket.gethostname(), started=now(), nice=args.nice,
               workers=[dict(pid=p, log=l) for p, l in pids],
               argv=["work"] + selection_argv(args))
    write_json_atomic(os.path.join(LOGS, "launch-%s.json" % tag), rec)
    for p, l in pids:
        print("pid %d  log %s" % (p, l))


def cmd_list(args):
    plan = load_plan()
    ps = select(plan, args)
    for p in ps:
        print(p["id"])
    print("# %d pieces, est. %.1f core-h" % (len(ps), sum(p["est_core_h"] for p in ps)),
          file=sys.stderr)


def cmd_build_caches(args):
    names = args.only.split(",") if args.only else list(CACHES)
    todo = [n for n in names if args.force or not cache_ok(n)]
    if not todo:
        print("all caches present and matching their records")
        return
    if args.detach:
        tag = args.tag or time.strftime("caches-%Y%m%d-%H%M%S")
        out = []
        for n in todo:
            lf = os.path.join(LOGS, "%s-%s.log" % (tag, n))
            pid = spawn_detached(["build-caches", "--only", n] + (["--force"] if args.force else []),
                                 lf, args.nice)
            out.append(dict(cache=n, pid=pid, log=rel(lf)))
            print("pid %d  cache %s  log %s" % (pid, n, rel(lf)))
        write_json_atomic(os.path.join(LOGS, "launch-%s.json" % tag),
                          dict(tag=tag, host=socket.gethostname(), started=now(),
                               nice=args.nice, builds=out))
        return
    for n in todo:
        build_one_cache(n)


def cmd_cache_hash(args):
    for f in args.files:
        print("%s  content=%s  file=%s" % (f, tube_content_hash(f), sha256_file(f)))


def cmd_plan(args):
    plan = build_plan()
    write_json_atomic(PLAN, plan)
    tot = sum(p["est_core_h"] for p in plan["pieces"])
    print("plan: %d pieces %s; %d excluded; est. %.0f core-h"
          % (len(plan["pieces"]), plan["counts"], len(plan["excluded"]), tot))


# ----------------------------------------------------- checks and Lemma D
def cmd_checks(args):
    import_p4()
    import p4.validated.tc_commutator as tc
    from p4.validated.tc_energy import Atoms
    t0 = time.time()
    at = Atoms(full=True)
    ta = time.time() - t0
    res, times = {}, {}
    for name, fn in (("E1_U6top", tc.check_u6_top_identity),
                     ("E2_layers", tc.check_layer_identity),
                     ("E3_sonic_zeros", tc.check_sonic_zeros),
                     ("E5_sonic_vanishing", tc.check_sonic_vanishing)):
        t1 = time.time()
        res[name] = fn(at)
        times[name] = time.time() - t1
        log("%s -> %r (%.1f s)" % (name, res[name], times[name]))
    ok = (res["E1_U6top"] is True and res["E2_layers"][1] is True
          and res["E2_layers"][0] >= 14 and all(res["E3_sonic_zeros"])
          and all(res["E5_sonic_vanishing"]))
    rec = dict(schema="s43-checks/1", ok=ok, results=res, seconds=times, atoms_s=ta,
               created=now(), env=env_info(),
               code=dict(files=loaded_code(), fingerprint=fingerprint(loaded_code()),
                         git_head=git_head()[0], git_dirty_src=git_head()[1]),
               note=("E2 returns (pmin, ok): the lowest centre layer present in Z and "
                     "whether the identity holds (no layers below the target, no odd "
                     "layers); the note's 'val >= 16' is the valuation of Z o u after "
                     "A - 1 = O(t^2), see check_layer_identity."))
    write_json_atomic(os.path.join(PARTS, "checks.json"), rec)
    print("checks ok=%s %s" % (ok, res))


def cmd_lemma_d(args):
    """certified: tc_commutator.run_lemma_d as committed (its _centre_expansion
    builds the certified (a*, mu*) balls since 2026-10-16; before that the
    driver substituted them).  midpoint: the pre-2026-10-16 construction
    (bare a3 midpoint), substituted here only to replay the old record."""
    import_p4()
    import p4.validated.tc_commutator as tc
    from p4.validated.tc_energy import Atoms
    at = Atoms(full=True)
    balls = args.centre_balls
    orig = tc._centre_expansion
    if balls == "midpoint":
        tc._centre_expansion = lambda prec=256: centre_expansion("midpoint", prec)
    t0, c0 = time.time(), time.process_time()
    try:
        r = tc.run_lemma_d(at)
    finally:
        tc._centre_expansion = orig
    rec_json = json.load(open(os.path.join(TC_DIR, "tc_energy.json"))).get("lemma_d")
    nbad = sum(1 for s in r["subranges"] if not s["ok"])
    rec = dict(schema="s43-lemma-d/1", ok=r["ok"], centre_balls=balls,
               centre_balls_note=("certified: a3_midpoint midpoint + K_minus_m "
                                  "(encloses the certified (a*, mu*)), the committed "
                                  "tc_commutator._centre_expansion since 2026-10-16; "
                                  "midpoint: the pre-2026-10-16 construction (radius "
                                  "<= 1e-41, does not enclose (a*, mu*))"),
               code_path=("tc_commutator.run_lemma_d as committed" if balls == "certified"
                          else "run_lemma_d with _centre_expansion replaced by the "
                               "midpoint construction"),
               c0=r["c0"], ycut=r["ycut"], layer_identity=r["layer_identity"],
               mins=r["mins"], nsub=len(r["subranges"]), failing_subranges=nbad,
               wall_s=time.time() - t0, cpu_s=time.process_time() - c0,
               recorded_in_tc_energy=rec_json, created=now(), env=env_info(),
               code=dict(files=loaded_code(), fingerprint=fingerprint(loaded_code()),
                         git_head=git_head()[0], git_dirty_src=git_head()[1]))
    out = os.path.join(PARTS, "lemma_d_%s.json" % balls)
    if args.replay:
        prev = json.load(open(out)) if os.path.exists(out) else None
        keys = ("ok", "c0", "ycut", "layer_identity", "mins", "nsub", "failing_subranges")
        rec["stored"] = ({k: prev.get(k) for k in keys} if prev else None)
        rec["stored_file"] = rel(out)
        rec["stored_code_fingerprint"] = (prev or {}).get("code", {}).get("fingerprint")
        rec["identical"] = bool(prev) and all(rec[k] == prev.get(k) for k in keys)
        out = os.path.join(PARTS, "replays", "lemma_d_%s.json" % balls)
    write_json_atomic(out, rec)
    print("lemma_d[%s] ok=%s mins=%s failing=%d wall=%.1fs%s -> %s"
          % (balls, r["ok"], r["mins"], nbad, rec["wall_s"],
             (" identical=%s" % rec["identical"]) if args.replay else "", rel(out)))


def cmd_p43_centre(args):
    """P4-3's centre tail t = e^x in (0, e^{-8}] (tc_monotone.run_centre): the
    committed code (certified (a*, mu*) balls since 2026-10-16), plus a replay
    of the pre-fix construction (bare a3 midpoint) for comparison with the
    record in results/theorem_c/tc_monotone.json."""
    import_p4()
    import p4.validated.ec_centre as ecc
    import p4.validated.tc_monotone as tm
    t0, c0 = time.time(), time.process_time()
    cert = tm.run_centre()
    t1 = time.time()
    orig = ecc.certified_centre_balls
    ecc.certified_centre_balls = midpoint_centre_balls
    try:
        mid = tm.run_centre()
    finally:
        ecc.certified_centre_balls = orig
    stored = json.load(open(os.path.join(TC_DIR, "tc_monotone.json"))).get("centre")
    cmp_keys = sorted(cert["mins"])
    rec = dict(schema="s43-p43-centre/1", ok=cert["ok"], x_tail=cert["x_tail"],
               certified=dict(ok=cert["ok"], mins=cert["mins"], wall_s=t1 - t0,
                              balls="a3_midpoint + K_minus_m (ec_centre."
                                    "certified_centre_balls; encloses (a*, mu*))"),
               midpoint_replay=dict(ok=mid["ok"], mins=mid["mins"],
                                    balls="bare a3 midpoint (pre-2026-10-16; does not "
                                          "enclose (a*, mu*))"),
               stored_record=stored,
               midpoint_replay_identical_to_stored=bool(stored) and all(
                   mid["mins"][k] == stored["mins"].get(k) for k in cmp_keys),
               min_margin=min(cert["mins"].values()),
               change_vs_midpoint={k: cert["mins"][k] - mid["mins"][k] for k in cmp_keys},
               claim=("P4-3 centre tail: v_rel > 0, v_rel' > 0, -rho^' > 0, "
                      "1/3 - v_rel^2 > 0, -t^2 D > 0, D' > 0 on t in (0, e^-8] "
                      "(all t-factored margins must be > 0)"),
               created=now(), env=env_info(),
               wall_s=time.time() - t0, cpu_s=time.process_time() - c0,
               code=dict(files=loaded_code(), fingerprint=fingerprint(loaded_code()),
                         git_head=git_head()[0], git_dirty_src=git_head()[1]))
    write_json_atomic(os.path.join(PARTS, "p43_centre.json"), rec)
    print("p43-centre: ok=%s mins=%s\n  midpoint replay ok=%s identical_to_stored=%s"
          % (cert["ok"], cert["mins"], mid["ok"], rec["midpoint_replay_identical_to_stored"]))


# ------------------------------------------------------------------- bench
BENCH = os.path.join(TC_DIR, "bench")
BENCH_MODES = {"tq": True, "direct": False}


def bench_path(pid, mode):
    return os.path.join(BENCH, "%s__%s.json" % (pid.replace("/", "__"), mode))


class PhaseTimer:
    """Times tc_commutator's phase functions during one analyze_piece call by
    wrapping the module globals that analyze_piece looks up at call time
    (sup_ratio_boxes and tau_quad_check: the tau-quadratic phase;
    direct_S_boxes: the direct phase).  Arguments and results pass through
    unchanged; the wrappers are removed on exit."""
    NAMES = ("sup_ratio_boxes", "tau_quad_check", "direct_S_boxes")

    def __init__(self, tc):
        self.tc, self.orig, self.stats, self.calls = tc, {}, {}, []

    def __enter__(self):
        for n in self.NAMES:
            f = getattr(self.tc, n)
            self.orig[n] = f
            setattr(self.tc, n, self._wrap(n, f))
        return self

    def __exit__(self, *exc):
        for n, f in self.orig.items():
            setattr(self.tc, n, f)
        return False

    def _wrap(self, name, f):
        def g(*a, **kw):
            t0, c0 = time.time(), time.process_time()
            try:
                out = f(*a, **kw)
            except Exception as e:
                self._note(name, a, t0, c0, None, "%s: %s" % (type(e).__name__, e))
                raise
            self._note(name, a, t0, c0, out, None)
            return out
        return g

    def _note(self, name, a, t0, c0, out, err):
        dt, dc = time.time() - t0, time.process_time() - c0
        s = self.stats.setdefault(name, dict(calls=0, wall_s=0.0, cpu_s=0.0))
        s["calls"] += 1
        s["wall_s"] += dt
        s["cpu_s"] += dc
        e = dict(phase=name, wall_s=dt, cpu_s=dc, err=err)
        if name == "direct_S_boxes":
            e["eps_cells"] = a[8] if len(a) > 8 else None
            if out is not None:
                e.update(ok=out[0], worst=out[1])
        elif name == "sup_ratio_boxes":
            e["families"] = sorted(a[0]) if a else None
        elif name == "tau_quad_check" and out is not None:
            e.update(ok=out[0], margin=out[1].get("margin"), err=out[1].get("err"))
        self.calls.append(e)
        log("  %s %s %.0f s%s" % (name, e.get("eps_cells") or e.get("families") or "",
                                  dt, "" if err is None else " ERR " + err))


def bench_one(p, mode, ctx):
    """run_tube's per-step body for one piece with tq = BENCH_MODES[mode],
    without degrade: analyze_piece is timed as called."""
    tc = ctx.tc
    from p4.validated.arbseries import precision
    from p4.validated.tc_tm import AtomKP, step_components
    tube = ctx.tube(p["cache"])
    steps = [sd for sd in tube.steps if float(sd.x) - sd.h >= p["x_min"] - 1e-12
             and float(sd.x) <= p["x_max"] + 1e-12]
    sd = steps[p["index"]]
    if float(sd.x) != p["x"] or float(sd.h) != p["h"]:
        raise RuntimeError("step mismatch for %s" % p["id"])
    tq = BENCH_MODES[mode]
    start, la0 = now(), list(os.getloadavg())
    t0, c0 = time.time(), time.process_time()
    err, rep = None, None
    with precision(320):                                  # run_tube's prec
        comps = step_components(sd)
        ak = AtomKP(ctx.at, comps, sd_shift=0)
        pr = tc.build_pairs(ak)
        cm = tc.commute(pr, pr["dT"][-1][0].trim().inv())
        t_setup, c_setup = time.time() - t0, time.process_time() - c0
        log("bench %s %s: setup %.0f s; analyze_piece(tq=%s) ..." % (p["id"], mode,
                                                                    t_setup, tq))
        with PhaseTimer(tc) as pt:
            t1, c1 = time.time(), time.process_time()
            try:
                rep = tc.analyze_piece(pr, cm, 0.0, float(sd.h), B1,
                                       tc.eps_cells_default(), 0, (0, 0, 0), tq=tq)
            except (ArithmeticError, RuntimeError) as e:
                err = "%s: %s" % (type(e).__name__, e)
            t_an, c_an = time.time() - t1, time.process_time() - c1
    rep = rep or {}
    cells = rep.get("cells", [])
    margin = max((c.get("margin") for c in cells if c.get("margin") is not None),
                 default=None)
    tq_s = sum(pt.stats.get(n, {}).get("wall_s", 0.0)
               for n in ("sup_ratio_boxes", "tau_quad_check"))
    dir_s = pt.stats.get("direct_S_boxes", {}).get("wall_s", 0.0)
    files = loaded_code()
    head, dirty = git_head()
    return dict(
        schema="s43-bench/1", id=p["id"], spec=spec_of(p), mode=mode, tq=tq,
        call=("analyze_piece(pr, cm, 0.0, h, %g, eps_cells_default(), 0, (0, 0, 0), "
              "tq=%s): run_tube's per-step body at prec 320, no degrade" % (B1, tq)),
        result=dict(ok=bool(rep.get("ok")) and err is None, margin=margin, err=err,
                    open=rep.get("tau_quad"), cells=cells,
                    bounds={k: rep.get(k) for k in ("Rbar", "thU", "DW", "DWdec", "DWp",
                                                     "DB", "DBdec", "LL") if k in rep}),
        timing=dict(start=start, end=now(), setup_wall_s=t_setup, setup_cpu_s=c_setup,
                    analyze_wall_s=t_an, analyze_cpu_s=c_an,
                    tau_quadratic_phase_wall_s=tq_s, direct_phase_wall_s=dir_s,
                    loadavg_start=la0, loadavg_end=list(os.getloadavg()),
                    nice=os.nice(0)),
        phases=pt.stats, phase_calls=pt.calls,
        code=dict(files=files, fingerprint=fingerprint(files), git_head=head,
                  git_dirty_src=dirty),
        driver_sha256=sha256_file(os.path.abspath(__file__)), env=env_info(),
        note=("timing benchmark (recovery note 8(ii)); not a part: merge does not "
              "read results/theorem_c/bench/"))


def cmd_bench(args):
    plan = load_plan()
    ids = [s for s in args.ids.split(",") if s]
    modes = [m for m in args.modes.split(",") if m]
    for m in modes:
        if m not in BENCH_MODES:
            raise SystemExit("unknown mode %s (tq | direct)" % m)
    pieces = {p["id"]: p for p in plan["pieces"]}
    for pid in ids:
        if pid not in pieces or pieces[pid]["kind"] != "tube":
            raise SystemExit("bench needs planned tube-kind pieces: %s" % pid)
    if args.detach:
        tag = args.tag or time.strftime("bench-%Y%m%d-%H%M%S")
        out = []
        for pid in ids:                     # one detached process per piece
            lf = os.path.join(LOGS, "%s-%s.log" % (tag, pid.replace("/", "__")))
            av = ["bench", "--ids", pid, "--modes", ",".join(modes)]
            if args.force:
                av.append("--force")
            wp = spawn_detached(av, lf, args.nice)
            out.append(dict(piece=pid, pid=wp, log=rel(lf), argv=av,
                            results=[rel(bench_path(pid, m)) for m in modes]))
            print("pid %d  bench %s (%s)  log %s" % (wp, pid, ",".join(modes), rel(lf)))
        write_json_atomic(os.path.join(LOGS, "launch-%s.json" % tag),
                          dict(tag=tag, host=socket.gethostname(), started=now(),
                               nice=args.nice, bench=out))
        return
    ctx = None
    for pid in ids:
        for mode in modes:
            out = bench_path(pid, mode)
            if os.path.exists(out) and not args.force:
                log("bench %s %s: result present (%s), skipped" % (pid, mode, rel(out)))
                continue
            lk = claim("bench/%s/%s" % (pid, mode))
            if lk is None:
                log("bench %s %s: claimed by a live process, skipped" % (pid, mode))
                continue
            try:
                if ctx is None:
                    ctx = Ctx(REGIONS["centre"]["centre_balls"])
                rec = bench_one(pieces[pid], mode, ctx)
                write_json_atomic(out, rec)
                log("bench %s %s done: ok=%s margin=%s analyze %.0f s (tau-quadratic "
                    "phase %.0f s, direct phase %.0f s) -> %s"
                    % (pid, mode, rec["result"]["ok"], rec["result"]["margin"],
                       rec["timing"]["analyze_wall_s"],
                       rec["timing"]["tau_quadratic_phase_wall_s"],
                       rec["timing"]["direct_phase_wall_s"], rel(out)))
            except Exception:
                log("EXCEPTION in bench %s %s (no result written):\n%s"
                    % (pid, mode, traceback.format_exc()))
            finally:
                try:
                    os.remove(lk)
                except OSError:
                    pass
    log("bench exit")


# ------------------------------------------------------------------ status
def cmd_status(args):
    plan = load_plan()
    accept = parse_accept(getattr(args, "accept_code", None))
    rows, locks = {}, set()
    if os.path.isdir(LOCKS):
        locks = {f[:-5].replace("__", "/") for f in os.listdir(LOCKS) if f.endswith(".lock")}
    for p in plan["pieces"]:
        d = read_part(p["id"])
        st = part_state(p, d, accept)
        if st == "missing" and p["id"] in locks:
            st = "running"
        r = rows.setdefault(p["region"], dict(planned=0, cpu_h=0.0, est_left=0.0))
        r["planned"] += 1
        r[st] = r.get(st, 0) + 1
        if d is not None and st not in ("stale",):
            r["cpu_h"] += d["timing"].get("cpu_s", 0) / 3600.0
        if st in ("missing", "running", "stale"):
            r["est_left"] += p["est_core_h"]
        if args.verbose and st not in ("missing",):
            print("%-22s %-11s margin=%s" % (p["id"], st,
                                             (d or {}).get("result", {}).get("margin")))
    cols = ("planned", "pass-both", "pass-direct", "fail", "error", "stale",
            "running", "missing")
    print("%-7s " % "region" + " ".join("%11s" % c for c in cols)
          + "   cpu_h_done  est_core_h_left")
    for reg in REGION_ORDER:
        r = rows.get(reg, {})
        print("%-7s " % reg + " ".join("%11d" % r.get(c, 0) for c in cols)
              + "   %10.1f  %15.0f" % (r.get("cpu_h", 0), r.get("est_left", 0)))
    for f in sorted(glob.glob(os.path.join(PARTS, "*.json"))):
        if os.path.basename(f) in ("plan.json",):
            continue
        d = json.load(open(f))
        print("%-28s ok=%s" % (os.path.basename(f), d.get("ok")))
    for f in sorted(glob.glob(os.path.join(BENCH, "*.json"))):
        d = json.load(open(f))
        t = d["timing"]
        print("bench %-22s %-6s ok=%s margin=%s analyze=%.0fs (tq-phase %.0fs, direct "
              "%.0fs)" % (d["id"], d["mode"], d["result"]["ok"], d["result"]["margin"],
                          t["analyze_wall_s"], t["tau_quadratic_phase_wall_s"],
                          t["direct_phase_wall_s"]))
    for lk in sorted(locks):
        if lk.startswith("bench/"):
            print("bench running: %s" % lk[len("bench/"):])


def parse_accept(items):
    acc = {}
    for it in items or []:
        f, sha = it.split("=", 1)
        acc.setdefault(f, set()).add(sha)
    return acc


# ------------------------------------------------------------ cloud import
def cmd_import_cloud(args):
    """Branch A: the node's ~/out/cw_<tag>.json are run_tube reports (detail per
    step) from /tmp/cworker.py (tube_full.json chunks A0-A1, B0-B7; tube_deep.json
    chunks C0-C5).  Each detail entry becomes one part, matched to the plan by
    (cache, x, h).  --code-manifest is a JSON {relpath: sha256} of the sources the
    node EXECUTED (see notes/s4-3-recovery.md, Branch A, for how to derive it);
    --tube-hash full=<content sha>, deep=<content sha> are the node caches' content
    hashes (tube_content_hash), which must equal the local caches'."""
    plan = load_plan()
    code = json.load(open(args.code_manifest))
    th = dict(kv.split("=", 1) for kv in args.tube_hash or [])
    by_key = {}
    for p in plan["pieces"]:
        if p["kind"] == "tube":
            by_key[(p["cache"], p["x"], p["h"])] = p
    n_new = n_skip = 0
    unmatched = []
    for f in sorted(glob.glob(os.path.join(args.dir, "cw_*.json"))):
        tag = os.path.basename(f)[3:-5]
        cache = "deep" if tag.startswith("C") else "full"
        if cache in th and th[cache] != plan["caches"][cache]["content_sha256"]:
            raise SystemExit("node cache %s content hash differs from the local cache: "
                             "the imported steps would not be the planned ones" % cache)
        d = json.load(open(f))
        logf = os.path.join(args.dir, "cw_%s.log" % tag)
        for e in d["detail"]:
            p = by_key.get((cache, e["x"], e["h"]))
            if p is None:
                unmatched.append((tag, e["x"], e["h"]))
                continue
            if read_part(p["id"]) is not None and not args.overwrite:
                n_skip += 1
                continue
            rec = dict(schema=SCHEMA, id=p["id"], region=p["region"], spec=spec_of(p),
                       status=classify(e), result=e, b1=B1, tstar=16 * B1,
                       extra=dict(region_call="cworker.py chunk %s: run_tube(at, %s, 14.0, "
                                  "degrade=True)" % (tag, cache),
                                  chunk=dict(tag=tag, steps=d.get("pieces"),
                                             wall_s=d.get("time"), ok=d.get("ok"),
                                             worst=d.get("worst_margin"))),
                       code=dict(files=code, fingerprint=fingerprint(code),
                                 git_head=None, git_dirty_src=None),
                       driver_sha256=None, env=dict(host="cloud " + args.node),
                       timing=dict(start=None, end=None,
                                   wall_s=(d.get("time") or 0) / max(1, d.get("pieces", 1)),
                                   cpu_s=(d.get("time") or 0) / max(1, d.get("pieces", 1)),
                                   note="chunk wall time / steps (one dedicated vCPU)"),
                       provenance="imported from %s (%s); node log %s" % (
                           os.path.basename(f), sha256_file(f),
                           os.path.basename(logf) if os.path.exists(logf) else "absent"))
            write_json_atomic(part_path(p["id"]), rec)
            n_new += 1
    print("imported %d parts (%d already present); %d unmatched steps %s"
          % (n_new, n_skip, len(unmatched), unmatched[:6]))


PROGRESS_BLOCK_MARK = 'if os.environ.get("P4_PROGRESS"):'


def strip_progress_blocks(src):
    """tc_commutator.py as the September cloud workers executed it: the
    committed file minus the two `if os.environ.get("P4_PROGRESS"): print(...)`
    blocks that were added on 2026-09-05T01:13Z, after the cworker processes
    (started 2026-09-04T01:30Z) had imported the module."""
    lines = src.splitlines(keepends=True)
    out, i, removed = [], 0, 0
    while i < len(lines):
        ln = lines[i]
        if ln.strip() == PROGRESS_BLOCK_MARK:
            ind = len(ln) - len(ln.lstrip())
            i += 1
            while i < len(lines) and (len(lines[i]) - len(lines[i].lstrip()) > ind
                                      or not lines[i].strip()):
                i += 1
            removed += 1
            continue
        out.append(ln)
        i += 1
    return "".join(out), removed


def cmd_cloud_manifest(args):
    """Code manifest for import-cloud: hashes of the committed p4 sources the
    cworker processes loaded, with tc_commutator.py replaced by its
    logging-stripped variant (written next to the manifest for diffing)."""
    import_p4()
    files = loaded_code()
    tcp = os.path.join(SRC, "p4", "validated", "tc_commutator.py")
    var, removed = strip_progress_blocks(open(tcp).read())
    if removed != 2:
        raise SystemExit("expected 2 P4_PROGRESS blocks, found %d" % removed)
    os.makedirs(LOGS, exist_ok=True)
    vpath = os.path.join(LOGS, "tc_commutator.executed-2026-09-04.py")   # for diffing
    with open(vpath, "w") as fh:
        fh.write(var)
    files[rel(tcp)] = hashlib.sha256(var.encode()).hexdigest()
    write_json_atomic(args.out, files)
    print("manifest %s: %d files; executed tc_commutator.py variant sha256 %s (%s)"
          % (rel(args.out), len(files), files[rel(tcp)], rel(vpath)))


def cmd_replay(args):
    """Recompute one piece and compare with its stored part (never overwrites
    the part); the replay record goes to parts/replays/<id>.json."""
    plan = load_plan()
    ps = [p for p in plan["pieces"] if p["id"] == args.id]
    if not ps:
        raise SystemExit("unknown piece " + args.id)
    p = ps[0]
    stored = read_part(p["id"])
    ctx = Ctx(REGIONS["centre"]["centre_balls"])
    t0, c0 = time.time(), time.process_time()
    res, extra = run_piece(p, ctx)
    same = stored is not None and all(res.get(k) == stored["result"].get(k)
                                      for k in ("ok", "margin", "open", "err"))
    rec = dict(schema="s43-replay/1", id=p["id"], result=res, extra=extra,
               stored=(stored or {}).get("result"),
               stored_provenance=(stored or {}).get("provenance"), identical=same,
               code=dict(files=loaded_code()), env=env_info(),
               timing=dict(wall_s=time.time() - t0, cpu_s=time.process_time() - c0),
               created=now())
    write_json_atomic(os.path.join(PARTS, "replays", p["id"] + ".json"), rec)
    print("replay %s: identical=%s result=%s stored=%s"
          % (p["id"], same, res, (stored or {}).get("result")))


# ------------------------------------------------------------------- merge
def x_interval(p):
    if p["kind"] == "centre":
        return (-math.inf if p["lo"] <= 0 else math.log(p["lo"]), math.log(p["hi"]))
    if p["kind"] == "sonic":
        return (-p["hi"], -p["lo"])
    return (p["x"] - p["h"], p["x"])


def uncovered(intervals, lo=-math.inf, hi=0.0, tol=1e-12):
    """Sub-intervals of [lo, hi] not covered by the union of closed intervals."""
    gaps, cur = [], lo
    for a, b in sorted(intervals):
        if b < cur:
            continue
        if a > cur + tol:
            gaps.append((cur, a))
        cur = max(cur, b)
        if cur >= hi:
            break
    if cur < hi - tol:
        gaps.append((cur, hi))
    return gaps


def cmd_merge(args):
    import_p4()
    import p4.validated.tc_commutator as tc
    plan = load_plan()
    accept = parse_accept(args.accept_code)
    states, used = {}, []
    code_union = {}
    for p in plan["pieces"]:
        d = read_part(p["id"])
        st = part_state(p, d, accept)
        states[p["id"]] = st
        if st in ("pass-both", "pass-direct", "fail", "error"):
            for r, sha in d["code"]["files"].items():
                code_union.setdefault(r, set()).add(sha)
            if st.startswith("pass"):
                used.append((p, d))
    # several hashes for one source file among the used parts are acceptable only
    # when every non-current hash was explicitly accepted (--accept-code)
    conflicts = {}
    for r, s in code_union.items():
        cur = current_hash(r)
        bad = [x for x in s if x != cur and x not in accept.get(r, ())
               and x not in accept.get(os.path.basename(r), ())]
        if len(s) > 1 or bad:
            conflicts[r] = dict(hashes=sorted(s), current=cur, unaccepted=bad)
    unaccepted = {r: c for r, c in conflicts.items() if c["unaccepted"]}
    direct_iv = [x_interval(p) for p, d in used]
    both_iv = [x_interval(p) for p, d in used if d["status"] == "pass-both"]
    gaps_direct = uncovered(direct_iv)
    gaps_both = uncovered(both_iv)
    missing = sorted(i for i, s in states.items() if s in ("missing", "stale"))
    failed = sorted(i for i, s in states.items() if s in ("fail", "error"))
    checks = json.load(open(os.path.join(PARTS, "checks.json"))) \
        if os.path.exists(os.path.join(PARTS, "checks.json")) else None
    checks_ok = bool(checks and checks["ok"] and code_matches(checks["code"]["files"], accept))
    complete = (not gaps_direct and not missing and not failed and checks_ok
                and not unaccepted)
    opens = [d for p, d in used if d["status"] == "pass-direct"]
    xi_open = max([-a for a, b in gaps_both], default=0.0) if gaps_both else 0.0
    detail = []
    for p, d in used:
        e = dict(d["result"])
        e.update(id=p["id"], region=p["region"])
        detail.append(e)
    worst = max((e["margin"] for e in detail if e.get("margin") is not None), default=None)
    res = dict(
        ok=complete, status="complete" if complete else "partial",
        scope=dict(full_strip_tau=[B1, tc.TSTAR_FAC * B1],
                   tau_le_tstar_uncovered_x=gaps_direct,
                   tau_gt_tstar_uncovered_x=gaps_both,
                   xi_open_max=xi_open,
                   tau_above_tstar=("certified everywhere" if not gaps_both else
                                    "OPEN on x in %s (direct phase only there)" % gaps_both)),
        b1=B1, m=tc.MC, c=[tc.CNUM, tc.CDEN], eps_max=tc.EPS_MAX, tstar_fac=tc.TSTAR_FAC,
        phase_settings=dict(
            tstar_fac_by_region={r: REGIONS[r]["tstar_fac"] for r in REGIONS
                                 if "tstar_fac" in REGIONS[r]},
            direct_only={r: REGIONS[r]["direct_only"] for r in REGIONS
                         if "direct_only" in REGIONS[r]},
            note=("pass-both pieces are certified for every tau >= b1 (centre with "
                  "tstar_fac = 1: the tau-quadratic phase alone); pass-direct pieces "
                  "only for tau in [b1, tstar_fac * b1]")),
        open_pieces=len(opens),
        checks=(dict(E1_U6top=checks["results"]["E1_U6top"],
                     E2_layers=checks["results"]["E2_layers"],
                     E3_sonic_zeros=checks["results"]["E3_sonic_zeros"],
                     E5_sonic_vanishing=checks["results"]["E5_sonic_vanishing"])
                if checks else None),
        checks_valid_for_code=checks_ok,
        pieces=len(detail), worst_margin=worst, detail=detail,
        missing=missing, failed=failed, code_conflicts=conflicts,
        accepted_code_variants={k: sorted(v) for k, v in accept.items()},
        accept_reason=args.accept_reason,
        provenance=dict(plan_sha256=sha256_file(PLAN), caches=plan["caches"],
                        code=loaded_code(), git_head=git_head()[0], merged=now(),
                        driver_sha256=sha256_file(os.path.abspath(__file__)),
                        parts_dir=rel(PARTS)))
    out = CERT if complete else CERT_PARTIAL
    write_json_atomic(out, res)
    print("%s: status=%s used=%d missing=%d failed=%d gaps(tau<=tstar)=%s "
          "gaps(tau>tstar)=%s conflicts=%s -> %s"
          % ("merge", res["status"], len(used), len(missing), len(failed), gaps_direct,
             gaps_both, list(conflicts), rel(out)))


# -------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    def sel(p):
        p.add_argument("--regions")
        p.add_argument("--ids", help="comma-separated fnmatch patterns on piece ids")
        p.add_argument("--ids-file", help="file with one piece id per line")
        p.add_argument("--exclude-ids-file")
        p.add_argument("--x-range", help="lo,hi: pieces whose x-cover lies inside")
        p.add_argument("--shard", help="k/m: every m-th piece starting at k")
        p.add_argument("--order", default="cost", choices=("cost", "plan", "deep-first"))
        p.add_argument("--retry-failed", action="store_true")
        p.add_argument("--accept-code", action="append",
                       help="FILE=SHA: treat this source hash as current")
        p.add_argument("--max-pieces", type=int, default=0)

    p = sub.add_parser("build-caches")
    p.add_argument("--only")
    p.add_argument("--force", action="store_true")
    p.add_argument("--detach", action="store_true")
    p.add_argument("--nice", type=int, default=10)
    p.add_argument("--tag")
    p.set_defaults(fn=cmd_build_caches)
    p = sub.add_parser("cache-hash")
    p.add_argument("files", nargs="+")
    p.set_defaults(fn=cmd_cache_hash)
    p = sub.add_parser("plan")
    p.set_defaults(fn=cmd_plan)
    p = sub.add_parser("checks")
    p.set_defaults(fn=cmd_checks)
    p = sub.add_parser("lemma-d")
    p.add_argument("--centre-balls", default="certified", choices=("certified", "midpoint"))
    p.add_argument("--replay", action="store_true",
                   help="write parts/replays/lemma_d_<balls>.json and compare")
    p.set_defaults(fn=cmd_lemma_d)
    p = sub.add_parser("p43-centre")
    p.set_defaults(fn=cmd_p43_centre)
    p = sub.add_parser("bench")
    p.add_argument("--ids", required=True)
    p.add_argument("--modes", default="tq,direct")
    p.add_argument("--detach", action="store_true")
    p.add_argument("--nice", type=int, default=10)
    p.add_argument("--tag")
    p.add_argument("--force", action="store_true")
    p.set_defaults(fn=cmd_bench)
    p = sub.add_parser("work")
    sel(p)
    p.set_defaults(fn=cmd_work)
    p = sub.add_parser("launch")
    sel(p)
    p.add_argument("--workers", type=int, required=True)
    p.add_argument("--nice", type=int, default=10)
    p.add_argument("--tag")
    p.set_defaults(fn=cmd_launch)
    p = sub.add_parser("list", help="print the selected piece ids (for --ids-file splits)")
    sel(p)
    p.set_defaults(fn=cmd_list)
    p = sub.add_parser("status")
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument("--accept-code", action="append")
    p.set_defaults(fn=cmd_status)
    p = sub.add_parser("import-cloud")
    p.add_argument("dir")
    p.add_argument("--code-manifest", required=True)
    p.add_argument("--tube-hash", action="append")
    p.add_argument("--node", default="ubuntu@216.86.163.117")
    p.add_argument("--overwrite", action="store_true")
    p.set_defaults(fn=cmd_import_cloud)
    p = sub.add_parser("cloud-manifest")
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_cloud_manifest)
    p = sub.add_parser("replay")
    p.add_argument("id")
    p.set_defaults(fn=cmd_replay)
    p = sub.add_parser("merge")
    p.add_argument("--accept-code", action="append")
    p.add_argument("--accept-reason")
    p.set_defaults(fn=cmd_merge)
    args = ap.parse_args(argv)
    args.fn(args)


if __name__ == "__main__":
    main()
