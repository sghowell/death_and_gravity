#!/usr/bin/env python3
"""Diagnose the tau-quadratic check on one centre piece (S4-3 recovery note,
section 8 (iii): the centre with the tau-quadratic phase from tau = b1).

Runs tc_commutator.centre_piece on the committed code with the given
tstar_fac and records, per eps cell, what tau_quad_check sees: the offsets of
K2 (the tau^2 coefficient) and R1 (the tau-linear remainder slack), the
y-scaled coefficients Ah, Bh, the vertex Bh / (2|Ah|), K0 + R0 and the
margin.  tau_quad_check is wrapped only to read its inputs; the computation
is unchanged.

  python problems/P4/scripts/s43_centre_tq_diag.py centre/c29 --tstar-fac 1
      -> results/theorem_c/parts/experiments/centre-tstar_fac-1/c29_tau_quad_diag.json
"""
import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s43_parts as sp  # noqa: E402  (sets sys.path for p4)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("id")
    ap.add_argument("--tstar-fac", type=float, default=1.0)
    args = ap.parse_args(argv)
    plan = sp.load_plan()
    p = [q for q in plan["pieces"] if q["id"] == args.id]
    if not p or p[0]["kind"] != "centre":
        raise SystemExit("needs a planned centre piece")
    p = p[0]
    from flint import arb  # noqa: F401
    import p4.validated.tc_commutator as tc
    from p4.validated.ballutil import strip_ball
    from p4.validated.tc_energy import Atoms
    at = Atoms(full=True)
    ce = tc._centre_expansion()
    rows, orig = [], tc.tau_quad_check

    def wrap(K2, R1, K0, R0, b1, lo, hi, so):
        K2t = K2.trim()
        el = strip_ball(lo, hi)
        hK2, hR1 = K2t.off + so, R1.min_off() + so
        h = min(hK2, hR1)
        A = K2t.tm.eval(lo, hi)
        Ah = (A * abs(el) ** (hK2 - h)).upper()
        Bh = arb(0)
        for v, o in R1.t:
            Bh += arb(v) * abs(el) ** (o + so - h)
        Bh = Bh.upper()
        ok, m = orig(K2, R1, K0, R0, b1, lo, hi, so)
        rows.append(dict(K2_offset=hK2, R1_min_offset=hR1, h=h,
                         K2hat_upper=float(A.upper()), Ah=float(Ah), Bh=float(Bh),
                         vertex_tau=float(Bh) / (2 * abs(float(Ah))) if float(Ah) < 0 else None,
                         K0R0=m.get("K0R0"), margin=m.get("margin"), ok=ok,
                         R1_terms=[[float(arb(v).upper()), o] for v, o in R1.t]))
        return ok, m
    tc.tau_quad_check = wrap
    rd = sp.REGIONS["centre"]
    t0 = time.time()
    try:
        rep = tc.centre_piece(at, ce, p["lo"], p["hi"], sp.B1,
                              tc.eps_cells_default(rd["ncells"]), rd["prec"],
                              tstar_fac=args.tstar_fac)
    finally:
        tc.tau_quad_check = orig
    wall = time.time() - t0
    for r, c in zip(rows, rep["cells"]):
        r["eps"] = c["eps"]
    files = sp.loaded_code()
    head, dirty = sp.git_head()
    rec = dict(schema="s43-centre-tq-diag/1", id=p["id"], lo=p["lo"], hi=p["hi"],
               b1=sp.B1, tstar_fac=args.tstar_fac, tstar=args.tstar_fac * sp.B1,
               ok=rep["ok"], wall_s=wall,
               bounds={k: rep.get(k) for k in ("Rbar", "thU", "DW", "DWdec", "DWp", "DB",
                                                "DBdec", "LL")},
               cells=rows, created=sp.now(), env=sp.env_info(),
               code=dict(files=files, fingerprint=sp.fingerprint(files), git_head=head,
                         git_dirty_src=dirty),
               script_sha256=sp.sha256_file(os.path.abspath(__file__)))
    out = os.path.join(sp.PARTS, "experiments", "centre-tstar_fac-%g" % args.tstar_fac,
                       "%s_tau_quad_diag.json" % p["id"].split("/")[1])
    sp.write_json_atomic(out, rec)
    print("%s tstar_fac=%g ok=%s wall=%.1fs -> %s" % (p["id"], args.tstar_fac, rep["ok"],
                                                     wall, sp.rel(out)))
    for r in rows:
        print("  eps %-12s K2 off %d, R1 off %d: Ah %.3g Bh %.3g vertex %s K0R0 %.4g margin %.4g"
              % (r["eps"], r["K2_offset"], r["R1_min_offset"], r["Ah"], r["Bh"],
                 "%.1f" % r["vertex_tau"] if r["vertex_tau"] else None, r["K0R0"],
                 r["margin"]))


if __name__ == "__main__":
    main()
