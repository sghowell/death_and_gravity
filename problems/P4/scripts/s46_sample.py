#!/usr/bin/env python3
"""S4-6 check: independent point samples of the centre Sigma_c against the
committed centre parts.  Not part of the certificate; notes/s4-6-centre.md.

At random points (y, eps, tau) of each centre piece (y log-uniform in the
piece, from 1e-12 on the closing piece; eps uniform in [0, 6] plus both
ends; tau log-uniform in [14, 1e9] plus 14 and 224), three enclosures of
Sigma_c = S/(1 + (y tau)^2), s = y tau:

  graded      tc_layer_centre.sigma_box_centre on the point box: the
              certificate's evaluator (pivot forms, grading, raised nB) on
              the layer data the parts use (certified (a*, mu*) balls);
  direct      tc_layer_centre.point_sigma_direct: the direct chain of
              tc_layer.point_S_direct (A0 = nT/dT, B0 = nB/dB, ...; no pivot
              forms, no grading, no raise) from point pairs at the certified
              balls.  Rigorous, wide where s is small;
  direct_mid  the same chain at the background of the midpoints of the
              certified (a*, mu*) balls (a point of the certified box, so
              every enclosure over the box covers it).  Tight.

Per point it checks that graded meets direct and meets direct_mid (they
enclose the same value), and that direct_mid lies below the piece's
certified bound: the centre_t14 bound, and for tau >= 224 also the
centre_t224 bound.  A failure of any check is a bug in the certificate or
in the check.

usage: s46_sample.py [--pieces c29,c00,...] [--seed 20261017] [--ny 6]
                     [--neps 6] [--ntau 8] [--out FILE]
Output: problems/P4/results/theorem_c/centre_sampling.json (default).
"""
from __future__ import annotations

import argparse
import copy
import json
import math
import os
import random
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import s43_parts as s43  # noqa: E402

DRIVER = os.path.abspath(__file__)
DRIVER_SHA256 = s43.sha256_file(DRIVER)
OUT = os.path.join(s43.TC_DIR, "centre_sampling.json")
SCHEMA = "s46-centre-sampling/1"
Y_FLOOR = 1e-12          # the closing piece is sampled on [Y_FLOOR, hi]
TAU_MAX = 1e9


def midpoint_expansion(prec=256, K=30):
    """The centre expansion at the midpoints of the certified (a*, mu*)
    balls (ec_centre.certified_centre_balls): one background inside the
    certified box, with point coefficients."""
    from flint import arb
    from p4.validated import centre
    from p4.validated.arbseries import precision
    from p4.validated.ec_centre import certified_centre_balls
    with precision(prec):
        a, mu = certified_centre_balls()
        a, mu = arb(a.mid()), arb(mu.mid())
        ce = centre.centre_expansion(mu * (2 * a).exp(), nhat=(-a).exp(), K=K)
        if not ce.certify().ok:
            raise RuntimeError("midpoint centre certificate failed")
    return ce


def sample_points(p, rng, ny, neps, ntau):
    lo, hi = p["lo"], p["hi"]
    a = math.log(lo) if lo > 0 else math.log(Y_FLOOR)
    ys = [math.exp(rng.uniform(a, math.log(hi))) for _ in range(ny)]
    eps = [0.0, 6.0] + [rng.uniform(0.0, 6.0) for _ in range(neps)]
    taus = [14.0, 224.0] + [math.exp(rng.uniform(math.log(14.0), math.log(TAU_MAX)))
                            for _ in range(ntau)]
    return ys, eps, taus


def _f(x):
    return float(x)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--pieces", default=None)
    ap.add_argument("--seed", type=int, default=20261017)
    ap.add_argument("--ny", type=int, default=6)
    ap.add_argument("--neps", type=int, default=6)
    ap.add_argument("--ntau", type=int, default=8)
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    t0 = time.time()
    s43.import_p4()
    import p4.validated.tc_energy as te
    import p4.validated.tc_layer_centre as tcc
    from p4.validated.arbseries import precision
    from p4.validated.ec_centre import certified_centre_expansion
    s43.snapshot_code()
    at = te.Atoms(full=True)
    ce = certified_centre_expansion(256, K=30)
    cm = midpoint_expansion(256, K=30)
    if not tcc.check_nb_identity(at)["ok"]:
        raise SystemExit("the centre nB layer identity failed")
    plan = s43.load_plan()
    ps = [p for p in plan["pieces"] if p["kind"] == "centre"]
    if args.pieces:
        want = {x.strip() for x in args.pieces.split(",")}
        ps = [p for p in ps if p["id"].split("/")[1] in want or p["id"] in want]
    rng = random.Random(args.seed)
    out = dict(schema=SCHEMA, created=s43.now(), seed=args.seed,
               design=dict(ny=args.ny, neps=args.neps + 2, ntau=args.ntau + 2,
                           y_floor=Y_FLOOR, tau_max=TAU_MAX),
               checks=("graded meets direct; graded meets direct_mid; direct_mid upper "
                       "<= the centre_t14 bound, and for tau >= 224 <= the centre_t224 "
                       "bound"),
               pieces={})
    tot = dict(points=0, failures=0)
    for p in sorted(ps, key=lambda q: q["index"]):
        pid = p["id"]
        b14 = (s43.read_layer_part(s43.CENTRE_DIR, pid) or {}).get("certified_bound")
        b224 = (s43.read_layer_part(s43.CENTRE_DIR_224, pid) or {}).get("certified_bound")
        tb = time.time()
        ld = tcc.centre_layer_data(at, ce, p["hi"], label=pid, raised=True)
        tb = time.time() - tb
        ys, eps, taus = sample_points(p, rng, args.ny, args.neps, args.ntau)
        rows, fails = [], []
        worst = dict(direct_mid=None, graded_upper=None)
        n_tight = 0
        with precision(320):
            for y in ys:
                pr = tcc.point_pairs_centre(at, ce, y)
                pm = tcc.point_pairs_centre(at, cm, y)
                for e in eps:
                    for tau in taus:
                        r = dict(y=y, eps=e, tau=tau, s=y * tau)
                        try:
                            g = tcc.sigma_box_centre(ld, tcc.point_box(y, e, tau))
                        except (ArithmeticError, ZeroDivisionError, ValueError) as ex:
                            g, r["graded_error"] = None, str(ex)[:80]
                        d, _ = tcc.point_sigma_direct(pr, y, e, tau)
                        m, _ = tcc.point_sigma_direct(pm, y, e, tau)
                        why = []
                        if g is not None:
                            r.update(graded=_f(g.mid()), graded_rad=_f(g.rad()))
                            if not g.overlaps(d):
                                why.append("graded and direct disjoint")
                            if not g.overlaps(m):
                                why.append("graded and direct_mid disjoint")
                            gu = _f(g.upper())
                            if worst["graded_upper"] is None or gu > worst["graded_upper"]:
                                worst["graded_upper"] = gu
                        else:
                            why.append("graded undefined at a point")
                        if d.is_finite():
                            r.update(direct=_f(d.mid()), direct_rad=_f(d.rad()))
                        r.update(direct_mid=_f(m.mid()), direct_mid_rad=_f(m.rad()))
                        if m.is_finite() and _f(m.rad()) < 1e-6:
                            n_tight += 1
                        mu_ = _f(m.upper()) if m.is_finite() else math.inf
                        if b14 is None or not mu_ <= b14:
                            why.append("direct_mid above the centre_t14 bound")
                        if tau >= 224.0 and (b224 is None or not mu_ <= b224):
                            why.append("direct_mid above the centre_t224 bound")
                        if worst["direct_mid"] is None or mu_ > worst["direct_mid"]:
                            worst["direct_mid"] = mu_
                        if why:
                            r["fail"] = why
                            fails.append(r)
                        rows.append(r)
        tot["points"] += len(rows)
        tot["failures"] += len(fails)
        out["pieces"][pid] = dict(
            lo=p["lo"], hi=p["hi"], bound_t14=b14, bound_t224=b224, points=len(rows),
            direct_mid_tight=n_tight, failures=len(fails), fail_rows=fails[:20],
            max_direct_mid=worst["direct_mid"], max_graded_upper=worst["graded_upper"],
            build_s=tb, rows=rows)
        s43.log("sample %s: %d points, %d failures, max direct_mid %.6g, bound %s, "
                "build %.1fs" % (pid, len(rows), len(fails), worst["direct_mid"], b14, tb))
    out.update(points=tot["points"], failures=tot["failures"], ok=tot["failures"] == 0,
               wall_s=time.time() - t0, env=s43.env_info(),
               code=s43.code_record({DRIVER: DRIVER_SHA256}), driver_sha256=DRIVER_SHA256)
    path = args.out or OUT
    s43.write_json_atomic(path, out)
    print("sampling: %d points, %d failures -> %s" % (tot["points"], tot["failures"],
                                                       s43.rel(path)))


if __name__ == "__main__":
    main()
