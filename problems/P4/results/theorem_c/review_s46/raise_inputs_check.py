#!/usr/bin/env python3
"""S4-6 review: the inputs of the centre parts and the offset raise, on the
certified balls (not on the exact polynomials: identity_independent.py does
those).

For the certified centre background (ec_centre) it checks
  * the (a, mu) balls: midpoint + Krawczyk offset of results/a3_midpoint.json,
    re-read here from the JSON (own parsing), and that they exclude the bare
    a3 midpoint (the pre-2026-10-16 bug);
  * the certified series: every odd coefficient of n, w, v exactly 0 for
    k <= K; the tail radius nu against the largest TM domain YCUT;
  * tc_tm.centre_components: offsets (0, -1, 2, 1), A = 1 + y^2 Ahat with
    Ahat's TM coefficients overlapping the system's own Ahat series
    (centre.CentreExpansion.Ahat_series, from systems.centre_Ahat_fraction);
and, per piece, on build_pairs(AtomKP(at, centre_components(ce, hi))):
  * nB's kappa^j coefficients in all three jet slots: trimmed offset, the TM's
    exact-coefficient range ve, whether the two coefficients the raise drops
    contain 0, their radii, and the size of the first kept coefficient;
  * raise_nb's output equals dropping exactly those two coefficients;
  * the sonic derivative raise of tc_layer.LayerData ((tt dT_k)', which would
    drop coefficients): whether it fires on the centre data.

Run from the repository root (about 40 s per piece):
    PYTHONPATH=problems/P4/src python problems/P4/results/theorem_c/review_s46/raise_inputs_check.py [pieces] [out.json]
"""
import json
import os
import sys
import time

from flint import arb

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
A3 = os.path.join(ROOT, "problems", "P4", "results", "a3_midpoint.json")
PLAN = os.path.join(ROOT, "problems", "P4", "results", "theorem_c", "parts", "plan.json")


def parse_ball(s):
    mid, rad = s.strip()[1:-1].split(" +/- ")
    return arb(mid), arb(rad)


def main(pieces, out_path):
    t0 = time.time()
    from p4.validated import tc_commutator as tc
    from p4.validated import tc_energy as te
    from p4.validated import tc_layer_centre as tcc
    from p4.validated.arbseries import precision
    from p4.validated.ec_centre import certified_centre_balls, certified_centre_expansion
    from p4.validated.tc_tm import AtomKP, centre_components
    rec = dict(schema="review-s46-raise-inputs/1")
    with precision(256):
        aj = json.load(open(A3))
        am, ar = parse_ball(aj["a"])
        mm, mr = parse_ball(aj["mu"])
        oa, oar = parse_ball(aj["K_minus_m"][1])
        om, omr = parse_ball(aj["K_minus_m"][2])
        a, mu = certified_centre_balls()
        # the true (a*, mu*) lie in [mid - r, mid + r] + [off - r', off + r']
        lo_a, hi_a = am - ar + oa - oar, am + ar + oa + oar
        lo_m, hi_m = mm - mr + om - omr, mm + mr + om + omr
        rec["balls"] = dict(
            a=str(a), mu=str(mu), a_mid_a3=aj["a"], mu_mid_a3=aj["mu"],
            K_minus_m=aj["K_minus_m"],
            a_contains_box=bool(a.contains(lo_a) and a.contains(hi_a)),
            mu_contains_box=bool(mu.contains(lo_m) and mu.contains(hi_m)),
            a_excludes_bare_midpoint=not a.contains(am),
            mu_excludes_bare_midpoint=not mu.contains(mm),
            a_rad=float(a.rad()), mu_rad=float(mu.rad()))
    ce = certified_centre_expansion(256, K=30)
    K = ce.K
    odd_nonzero = [(k, i) for k in range(1, K + 1, 2) for i in range(3)
                   if not ce.balls[k][i].is_zero()]
    rec["series"] = dict(K=K, odd_coefficients_not_exactly_zero=odd_nonzero,
                         nu=str(ce.cert.nu), eps=str(ce.cert.eps), cert_ok=bool(ce.cert.ok),
                         nu_ge_ycut=bool(ce.cert.nu >= arb(tc.YCUT)),
                         ycut=tc.YCUT)
    # A = 1 + y^2 Ahat against the system's own Ahat series
    with precision(256):
        comps, (n, w, v) = centre_components(ce, tc.YCUT)
        As = ce.Ahat_series()
        Atm = comps[0].tm
        ah = [Atm.c[k + 2] for k in range(K - 1)]
        ser = [As.coeffs(K + 1)[k] if hasattr(As, "coeffs") else As[k] for k in range(K - 1)]
        mism = [k for k in range(K - 1) if not ah[k].overlaps(ser[k])]
        rec["components"] = dict(
            offsets=[cp.off for cp in comps],
            A_c0_is_1=bool(Atm.c[0] == 1 and Atm.c[0].rad() == 0),
            A_c1_is_exact_0=bool(Atm.c[1].is_zero()),
            Ahat_coefficients_compared=K - 1,
            Ahat_mismatch=mism,
            ve=[cp.tm.ve for cp in comps])
    at = te.Atoms(full=True)
    plan = json.load(open(PLAN))
    byid = {p["id"].split("/")[1]: p for p in plan["pieces"] if p["kind"] == "centre"}
    rec["pieces"] = {}
    for pid in pieces:
        p = byid[pid]
        tp = time.time()
        with precision(320):
            comps, _ = centre_components(ce, p["hi"])
            pr = tc.build_pairs(AtomKP(at, comps, sd_shift=0))
            rows = []
            for j, jet in enumerate(pr["nB"]):
                for sl, x in enumerate(jet):
                    t = x.trim()
                    c = t.tm.c
                    rows.append(dict(
                        kappa=j, slot=sl, off=t.off, ve=t.tm.ve,
                        c0_contains_0=bool(c[0].contains(0)),
                        c1_contains_0=bool(len(c) > 1 and c[1].contains(0)),
                        c0_rad=float(c[0].rad()), c1_rad=float(c[1].rad()) if len(c) > 1 else None,
                        c0_abs=float(abs(c[0]).upper()),
                        c2_abs_mid=float(abs(c[2]).mid()) if len(c) > 2 else None))
            raised = tcc.raise_nb(pr)
            same = True
            for j in range(len(pr["nB"])):
                for sl in range(3):
                    a_ = pr["nB"][j][sl].trim()
                    b_ = raised["nB"][j][sl].trim()
                    if j in range(1, 8):
                        want_c = a_.tm.c[2:]
                        ok = (b_.off == a_.off + 2 and len(b_.tm.c) == len(want_c)
                              and all(u.mid() == v.mid() and u.rad() == v.rad()
                                      for u, v in zip(b_.tm.c, want_c)))
                    else:
                        ok = b_.off == a_.off
                    same = same and ok
            # the LayerData (tt dT_k)' raise: does it fire on the centre data?
            s0 = lambda kp, sl: [x[sl].trim() for x in kp]
            dT0, dT1 = s0(raised["dT"], 0), s0(raised["dT"], 1)
            TNt = [raised["TN"][-1][k].trim() for k in (0, 1)]
            dTt = [dT0[-1], dT1[-1]]
            ti = dTt[0].inv()
            tt = (TNt[0] * ti).trim()
            ttp = ((TNt[1] * dTt[0] - TNt[0] * dTt[1]) * (ti * ti)).trim()
            fires = []
            for j in range(1, len(dT0)):
                val_off = (tt * dT0[j - 1]).trim().off
                tgt = max(val_off - 1, 0) if val_off >= 0 else val_off - 1
                der = (ttp * dT0[j - 1] + tt * dT1[j - 1]).trim()
                fires.append(dict(j=j, val_off=val_off, der_off=der.off, tgt=tgt,
                                  fires=der.off < tgt))
        r1to7 = [r for r in rows if 1 <= r["kappa"] <= 7]
        rec["pieces"][pid] = dict(
            lo=p["lo"], hi=p["hi"],
            nB_rows=rows,
            raise_ok=all(r["off"] == 8 and r["ve"] >= 2 and r["c0_contains_0"]
                         and r["c1_contains_0"] for r in r1to7),
            min_ve_kappa1to7=min(r["ve"] for r in r1to7),
            max_dropped_rad=max(max(r["c0_rad"], r["c1_rad"]) for r in r1to7),
            min_dropped_rad=min(min(r["c0_rad"], r["c1_rad"]) for r in r1to7),
            other_kappa_offsets=sorted({(r["kappa"], r["off"]) for r in rows
                                        if not 1 <= r["kappa"] <= 7}),
            raise_nb_equals_drop_of_two=same,
            ttdT_derivative_raise=fires,
            ttdT_raise_fires=any(f["fires"] for f in fires),
            seconds=time.time() - tp)
        print("%s: raise_ok=%s ve>=%d dropped rad in [%.3g, %.3g] same=%s fires=%s %.1fs"
              % (pid, rec["pieces"][pid]["raise_ok"], rec["pieces"][pid]["min_ve_kappa1to7"],
                 rec["pieces"][pid]["min_dropped_rad"], rec["pieces"][pid]["max_dropped_rad"],
                 same, rec["pieces"][pid]["ttdT_raise_fires"], time.time() - tp), flush=True)
    rec["seconds"] = time.time() - t0
    rec["ok"] = (rec["balls"]["a_contains_box"] and rec["balls"]["mu_contains_box"]
                 and rec["balls"]["a_excludes_bare_midpoint"]
                 and not odd_nonzero and rec["series"]["nu_ge_ycut"]
                 and rec["components"]["offsets"] == [0, -1, 2, 1]
                 and rec["components"]["A_c0_is_1"] and rec["components"]["A_c1_is_exact_0"]
                 and not rec["components"]["Ahat_mismatch"]
                 and all(v["raise_ok"] and v["raise_nb_equals_drop_of_two"]
                         and not v["ttdT_raise_fires"] for v in rec["pieces"].values()))
    with open(out_path, "w") as fh:
        json.dump(rec, fh, indent=1, sort_keys=True, default=str)
    print("ok=%s (%.0f s) -> %s" % (rec["ok"], rec["seconds"], out_path))


if __name__ == "__main__":
    ps = (sys.argv[1] if len(sys.argv) > 1 else "c00,c01,c15,c28,c29").split(",")
    main(ps, sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "raise_inputs_check.json"))
