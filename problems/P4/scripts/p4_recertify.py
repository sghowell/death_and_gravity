#!/usr/bin/env python3
"""Re-run the P4 certificates affected by the evaluation-ball fix
(p4.validated.ballutil; review docs/review-2026-10-17-p4-layer-engine-adversarial.md,
item 2f and section 5) on the current code, write the outputs to a new
directory, and compare their certified values, margins and counts with the
stored outputs.  Nothing stored is overwritten.

Certificates and their documented commands (notes/s2-theorem-b.md section 4,
notes/s4-3-recovery.md section 11.2, the module mains):
  checks          scripts/s43_parts.py checks             parts/checks.json
  lemma_d         scripts/s43_parts.py lemma-d --replay   parts/lemma_d_certified.json
  p43_centre      scripts/s43_parts.py p43-centre         parts/p43_centre.json
  tc_monotone     python -m p4.validated.tc_monotone TUBE results/theorem_c/tc_monotone.json
  tc_energy       python -m p4.validated.tc_energy all TUBE 6.0 (+ lemma_d)
                                                          results/theorem_c/tc_energy.json
  krawczyk        the snippet of s2-theorem-b.md section 4 results/theorem_b/krawczyk.json
  analyticity     python -m p4.validated.analyticity TUBE_B --workers 1
                                                          results/theorem_b/analyticity_R.json
  winding         python -m p4.validated.modecount TUBE_B --workers 1
                                                          results/theorem_b/winding_R.json
TUBE is the pinned S4-3 cache data/cache/tube_full.json (the tube tc_monotone.json
records); TUBE_B is the Theorem B tube, rebuilt by its documented command
(Tube.build(V0_EC, W_V0, x_c=-3.0), 244 steps) into data/cache/tube_theorem_b.json
(not committed; its sha256 is in OUTDIR/tube_theorem_b.check.json) and checked
against the first 244 steps of TUBE.

usage: p4_recertify.py run OUTDIR [--only a,b,..]     (one process at a time)
       p4_recertify.py compare OUTDIR --out FILE
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import s43_parts as s43  # noqa: E402

PY = sys.executable
P4 = s43.P4
SRC = s43.SRC
TUBE = s43.cache_path("full")
ORDER = ("checks", "lemma_d", "p43_centre", "tc_monotone", "krawczyk", "analyticity",
         "tc_energy", "winding")
STORED = dict(
    checks=os.path.join(s43.PARTS, "checks.json"),
    lemma_d=os.path.join(s43.PARTS, "lemma_d_certified.json"),
    p43_centre=os.path.join(s43.PARTS, "p43_centre.json"),
    tc_monotone=os.path.join(s43.TC_DIR, "tc_monotone.json"),
    tc_energy=os.path.join(s43.TC_DIR, "tc_energy.json"),
    krawczyk=os.path.join(P4, "results", "theorem_b", "krawczyk.json"),
    analyticity=os.path.join(P4, "results", "theorem_b", "analyticity_R.json"),
    winding=os.path.join(P4, "results", "theorem_b", "winding_R.json"))


def _env():
    return dict(os.environ, PYTHONPATH=SRC, PYTHONUNBUFFERED="1")


def _run(argv, log):
    s43.log("run: %s" % " ".join(argv))
    t0 = time.time()
    with open(log, "ab") as fh:
        r = subprocess.run(argv, cwd=s43.ROOT, env=_env(), stdout=fh, stderr=subprocess.STDOUT)
    s43.log("exit %d after %.0f s" % (r.returncode, time.time() - t0))
    if r.returncode:
        raise SystemExit("failed: %s (log %s)" % (" ".join(argv), log))


TUBE_B = os.path.join(s43.CACHE_DIR, "tube_theorem_b.json")   # a cache (data/cache, not in git)


def tube_b(outdir):
    """The Theorem B tube by its documented command, kept in data/cache (4 MB,
    not committed); checked against the first 244 steps of the pinned S4-3
    cache."""
    path = TUBE_B
    if not os.path.exists(path):
        code = ("from p4.validated import lintube, modecount; "
                "lintube.Tube.build(modecount.V0_EC, modecount.W_V0, x_c=-3.0).save(%r)" % path)
        _run([PY, "-c", code], os.path.join(outdir, "tube_theorem_b.log"))
    a, b = json.load(open(path)), json.load(open(TUBE))
    same = (a["K"] == b["K"] and len(a["steps"]) == 244
            and a["steps"] == b["steps"][:244])
    return path, same


def run(args):
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    only = args.only.split(",") if args.only else ORDER
    if not s43.cache_ok("full"):
        raise SystemExit("tube cache 'full' missing or not matching its record")
    s43s = os.path.join(HERE, "s43_parts.py")
    for name in ORDER:
        if name not in only:
            continue
        out = os.path.join(outdir, name + ".json")
        if os.path.exists(out) and not args.force:
            s43.log("skip %s (exists)" % name)
            continue
        log = os.path.join(outdir, name + ".log")
        if name == "checks":
            _run([PY, s43s, "checks", "--out", out], log)
        elif name == "lemma_d":
            _run([PY, s43s, "lemma-d", "--replay", "--out", out], log)
        elif name == "p43_centre":
            _run([PY, s43s, "p43-centre", "--out", out], log)
        elif name == "tc_monotone":
            code = ("import json, time; from p4.validated import tc_monotone as tm; "
                    "from p4.validated.lintube import Tube; t0 = time.time(); res = {}; "
                    "res['sonic'] = tm.run_sonic(); res['centre'] = tm.run_centre(); "
                    "tube = Tube.load(%r); "
                    "res['tube_info'] = {k: str(v) for k, v in tube.info.items()}; "
                    "res['tube'] = tm.run_tube(tube); "
                    "res['ok'] = all(res[r]['ok'] for r in ('sonic', 'centre', 'tube')); "
                    "res['total_time'] = time.time() - t0; "
                    "json.dump(res, open(%r, 'w'), indent=1)" % (TUBE, out))
            _run([PY, "-c", code], log)
        elif name == "tc_energy":
            code = ("import json, time; from p4.validated import tc_energy as te; "
                    "from p4.validated.tc_commutator import run_lemma_d; "
                    "from p4.validated.lintube import Tube; t0 = time.time(); c0 = 6.0; "
                    "at = te.Atoms(full=True); "
                    "res = dict(c0=c0, rbar=te.RBAR, ttail=te.TTAIL, atoms_time=time.time() - t0); "
                    "res['sonic'] = te.run_sonic(at, c0); "
                    "ld = run_lemma_d(at); ld['mins'] = {k: float(v) for k, v in ld['mins'].items()}; "
                    "ld.pop('subranges', None); res['lemma_d'] = ld; "
                    "tube = Tube.load(%r); "
                    "res['tube_info'] = {k: str(v) for k, v in tube.info.items()}; "
                    "res['tube'] = te.run_tube(at, tube, c0); "
                    "res['all_ok'] = all(res[r]['ok'] for r in ('sonic', 'lemma_d', 'tube')); "
                    "res['total_time'] = time.time() - t0; "
                    "json.dump(res, open(%r, 'w'), indent=1, default=str)" % (TUBE, out))
            _run([PY, "-c", code], log)
        elif name in ("krawczyk", "analyticity", "winding"):
            tb, same = tube_b(outdir)
            json.dump(dict(tube=s43.rel(tb), sha256=s43.sha256_file(tb),
                           same_as_first_244_steps_of_the_s43_cache=same,
                           s43_cache=s43.rel(TUBE), s43_cache_sha256=s43.sha256_file(TUBE)),
                      open(os.path.join(outdir, "tube_theorem_b.check.json"), "w"), indent=1)
            if name == "analyticity":
                _run([PY, "-m", "p4.validated.analyticity", tb, "--workers", "1", "--out", out], log)
            elif name == "winding":
                _run([PY, "-m", "p4.validated.modecount", tb, "--workers", "1", "--out", out], log)
            else:
                code = KRAWCZYK % dict(tube=tb, out=out)
                _run([PY, "-c", code], log)
    s43.log("run done -> %s" % outdir)


KRAWCZYK = r'''
import json, time
from flint import arb, acb
from p4.validated import modecount, krawczyk_kappa, lintube
from p4.validated.arbseries import abs_upper, precision
ctx = modecount.Context(lintube.Tube.load(%(tube)r))
out, zs = {}, {}
for key, kc, w in (("kappa1", modecount.KAPPA1, 2e-7), ("kbar", modecount.KGAUGE, 2e-5)):
    t0 = time.time()
    k = krawczyk_kappa.krawczyk(ctx, arb(kc), w)
    with precision(ctx.prec):
        z = zs[key] = k["zero"]
        out[key] = dict(ok=k["ok"], kappa_c=str(k["kappa_c"]), w_B=k["w_B"], m=k["model"].m, w=w,
                        E=str(k["E"]), E_abs_upper=float(abs(k["E"]).upper()),
                        E_abs_upper_arb=float(abs_upper(k["E"])), E_enc=modecount.enc(k["E"].real)
                        + modecount.enc(k["E"].imag),
                        dE=str(k["dE"]), dE_box=str(k["dE_box"]), rem=k["rem"],
                        zero=str(z), zero_re_interval=[str(z.real.lower()), str(z.real.upper())],
                        zero_im_abs=float(abs_upper(z.imag)), info=k["info"],
                        time=time.time() - t0)
with precision(ctx.prec):
    g = 1 / zs["kappa1"]
    out["gamma"] = dict(ball=str(g), re_interval=[str(g.real.lower()), str(g.real.upper())],
                        im_abs=float(abs_upper(g.imag)))
    kbar = krawczyk_kappa.gauge_value(ctx.bg)
    out["gauge"] = dict(kappa_bar=str(kbar),
                        contained=bool(zs["kbar"].real.contains(kbar)
                                       and zs["kbar"].imag.contains(arb(0))))
    out["disjoint"] = not zs["kappa1"].overlaps(zs["kbar"])
out["gauge_checks"] = str(krawczyk_kappa.gauge_checks(ctx.bg))[:4000]
json.dump(out, open(%(out)r, "w"), indent=1, default=str)
'''


def _load(p):
    return json.load(open(p)) if os.path.exists(p) else None


def compare(args):
    """Certified values, margins and counts of each new output against the
    stored one; identical means equal as stored (floats bit-equal, ball
    strings equal)."""
    outdir = os.path.abspath(args.outdir)
    rows = {}

    def cmp(name, pairs):
        diffs = {k: dict(stored=a, new=b) for k, (a, b) in pairs.items() if a != b}
        rows[name] = dict(identical=not diffs, compared=sorted(pairs), differences=diffs)

    for name in ORDER:
        old, new = _load(STORED[name]), _load(os.path.join(outdir, name + ".json"))
        if new is None:
            rows[name] = dict(identical=None, note="not run")
            continue
        if name == "checks":
            cmp(name, {k: (old[k], new[k]) for k in ("ok", "results")})
        elif name == "lemma_d":
            keys = ("ok", "c0", "ycut", "layer_identity", "mins", "nsub", "failing_subranges")
            cmp(name, {k: (old.get(k), new.get(k)) for k in keys})
        elif name == "p43_centre":
            cmp(name, {"ok": (old["ok"], new["ok"]),
                       "certified.mins": (old["certified"]["mins"], new["certified"]["mins"]),
                       "midpoint_replay.mins": (old["midpoint_replay"]["mins"],
                                                new["midpoint_replay"]["mins"]),
                       "min_margin": (old["min_margin"], new["min_margin"])})
        elif name == "tc_monotone":
            pairs = {"ok": (old["ok"], new["ok"])}
            for r in ("sonic", "centre", "tube"):
                for k in ("ok", "mins", "fails"):
                    if k in old[r] or k in new[r]:
                        pairs["%s.%s" % (r, k)] = (old[r].get(k), new[r].get(k))
            pairs["tube_info"] = (old.get("tube_info"), new.get("tube_info"))
            cmp(name, pairs)
        elif name == "tc_energy":
            pairs = {"all_ok": (old["all_ok"], new["all_ok"])}
            for k in ("ok", "mins", "boxes", "fails"):
                pairs["tube." + k] = (old["tube"].get(k), new["tube"].get(k))
            for k in ("ok", "mins", "boxes", "box_fails", "known_answer_sigma", "ED",
                      "xT_affine", "xT_rem_bound", "dc"):
                pairs["sonic." + k] = (old["sonic"].get(k), new["sonic"].get(k))
            for k in ("ok", "c0", "ycut", "layer_identity", "mins"):
                pairs["lemma_d." + k] = (old["lemma_d"].get(k), new["lemma_d"].get(k))
            pairs["c0"] = (old["c0"], new["c0"])
            cmp(name, pairs)
            # the stored lemma_d block predates the 2026-10-16 centre-ball fix
            # (ed4cab14: the centre series on the certified (a*, mu*) balls);
            # the certified-ball Lemma D is parts/lemma_d_certified.json
            ldc = _load(os.path.join(s43.PARTS, "lemma_d_certified.json")) or {}
            same_ldc = all(ldc.get(k) == new["lemma_d"].get(k)
                           for k in ("ok", "c0", "ycut", "layer_identity", "mins"))
            rows[name]["lemma_d_vs_parts_lemma_d_certified"] = dict(
                identical=same_ldc, file="parts/lemma_d_certified.json",
                note=("the stored tc_energy.json lemma_d block is the pre-2026-10-16 "
                      "midpoint-ball record (s4-3-recovery.md sections 3 and 11); the "
                      "re-run's block is the certified-ball Lemma D and equals "
                      "parts/lemma_d_certified.json"))
            rows[name]["note"] = ("tube_info differs by design: the stored record used a "
                                  "September build of the A3 tube (time %s); this run used the "
                                  "pinned S4-3 cache (time %s), the same construction; lemma_d "
                                  "nsub is not in this run's block (the Lemma D re-run records "
                                  "1100 subranges)"
                                  % (old["tube_info"].get("time"), new["tube_info"].get("time")))
        elif name == "krawczyk":
            pairs = {}
            for key in ("kappa1", "kbar"):
                for k in ("ok", "w_B", "w", "m", "E", "E_abs_upper", "dE", "dE_box", "rem",
                          "zero", "zero_re_interval", "zero_im_abs"):
                    pairs["%s.%s" % (key, k)] = (old[key].get(k), new[key].get(k))
                pairs[key + ".kappa_c"] = (old[key]["kappa_c"], new[key]["kappa_c"])
            pairs["gamma.ball"] = (old["gamma"]["ball"], new["gamma"]["ball"])
            pairs["gamma.re_interval"] = (old["gamma"]["re_interval"], new["gamma"]["re_interval"])
            pairs["gauge.kappa_bar"] = (old["gauge"]["kappa_bar"], new["gauge"]["kappa_bar"])
            pairs["gauge.contained"] = (old["gauge"]["contained"], new["gauge"]["contained"])
            pairs["disjoint"] = (old["disjoint"], new["disjoint"])
            cmp(name, pairs)
        elif name == "analyticity":
            keys = ("rect", "w0", "enlarge", "x_d", "constraint", "tube_regular", "n_steps",
                    "max_Pinv", "n_boxes", "min_w", "min_nu", "max_re_sigma", "max_Z_sonic",
                    "max_Z_centre", "max_g_centre", "tube_info")
            pairs = {k: (old.get(k), new.get(k)) for k in keys}
            # the pool returns the boxes in completion order: compare as multisets
            srt = lambda bs: sorted(json.dumps(b, sort_keys=True) for b in bs)
            pairs["boxes (as a multiset)"] = (srt(old["boxes"]), srt(new["boxes"]))
            cmp(name, pairs)
            if "tube_info" in rows[name]["differences"]:
                rows[name]["note"] = "tube_info.time is the wall time of the tube build"
        elif name == "winding":
            keys = ("rect", "m", "wmax", "winding", "certified", "winding_ball",
                    "n_tiles", "min_abs_E", "tube_info")
            pairs = {k: (old.get(k), new.get(k)) for k in keys}
            # per tile, matched by (kc, w): every field but the wall times
            # (time, t_prop) and 'width', a float diagnostic of the
            # propagation's set width (not a certified quantity)
            tk = lambda r: json.dumps([r["kc"], r.get("w")])
            skip = ("time", "t_prop", "width")
            to = {tk(r): {k: v for k, v in r.items() if k not in skip} for r in old["segments"]}
            tn = {tk(r): {k: v for k, v in r.items() if k not in skip} for r in new["segments"]}
            pairs["segments (all fields but time, t_prop, width)"] = (to, tn)
            wo = {tk(r): r.get("width") for r in old["segments"]}
            wn = {tk(r): r.get("width") for r in new["segments"]}
            relw = max((abs(wo[k] - wn[k]) / abs(wo[k]) for k in wo if k in wn and wo[k]),
                       default=0.0)
            # the pool sums the tile increments in completion order: compare
            # the exact sums in one fixed order (and report the stored total)
            from flint import arb
            from p4.validated.arbseries import precision
            from p4.validated.modecount import dec, enc
            with precision(256):
                fsum = lambda segs: enc(sum((dec(r["inc"]) for r in sorted(
                    segs, key=lambda r: json.dumps(r["inc"]))), arb(0)))
                pairs["total_arg (sum of the tile increments in a fixed order)"] = (
                    fsum(old["segments"]), fsum(new["segments"]))
            c_old = {k: v for k, v in old["constants"].items() if k != "workers"}
            c_new = {k: v for k, v in new["constants"].items() if k != "workers"}
            pairs["constants (without workers)"] = (c_old, c_new)
            cmp(name, pairs)
            rows[name]["width_max_relative_difference"] = relw
            probe = os.path.join(outdir, "winding_width_probe.json")
            if os.path.exists(probe):
                rows[name]["width_probe"] = dict(file=s43.rel(probe),
                                                 result=_load(probe).get("result"))
            rows[name]["total_arg_as_stored"] = dict(stored=old.get("total_arg"),
                                                     new=new.get("total_arg"))
            rows[name]["note_total_arg"] = (
                "the stored total_arg was summed by 8 pool workers in completion order, "
                "this one by 1; the midpoints agree and the radii differ in their last "
                "bits; summed in one fixed order the two sets of tile increments give "
                "bit-identical totals")
    # tube_info.time is a build time, not a certified value
    for name in ("analyticity", "winding"):
        r = rows.get(name) or {}
        d = (r.get("differences") or {}).get("tube_info")
        if d and d["stored"] and d["new"]:
            a = {k: v for k, v in d["stored"].items() if k != "time"}
            b = {k: v for k, v in d["new"].items() if k != "time"}
            if a == b:
                del r["differences"]["tube_info"]
                r["identical"] = not r["differences"]
                r["note"] = "tube_info equal except its build wall time"
    tb = _load(os.path.join(outdir, "tube_theorem_b.check.json"))
    rec = dict(schema="p4-recertify-compare/1", created=s43.now(), outdir=s43.rel(outdir),
               git_head=s43.git_head()[0], code=s43.loaded_code(),
               theorem_b_tube=tb, certificates=rows,
               all_identical=all(r.get("identical") for r in rows.values()
                                 if r.get("identical") is not None))
    s43.write_json_atomic(args.out, rec)
    for name, r in rows.items():
        print("%-12s identical=%s %s" % (name, r.get("identical"),
                                         sorted((r.get("differences") or {}).keys())))
    print("-> %s" % args.out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    q = sub.add_parser("run")
    q.add_argument("outdir")
    q.add_argument("--only")
    q.add_argument("--force", action="store_true")
    q.set_defaults(fn=run)
    q = sub.add_parser("compare")
    q.add_argument("outdir")
    q.add_argument("--out", required=True)
    q.set_defaults(fn=compare)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
