"""Follow-up (2026-10-19) on review finding S1: the S closing box's bounds on |P_N(rho)| at complex rho.

At the review head, `throat.certify_closing` bounds |P_N(rho)|, rho = lam rho_1 e^{2Y(t)}, through
`throat.nh_P_bounds` (cells) and `throat.nh_P_disc` (tail).  This script, using the author's code:
  (1) evaluates |P_N| in Arb at sampled points (t, lam) of every base cell of both closing boxes (the cell
      grid of certify_closing: multiples of 1/16 on [-8, t_a] plus the design's breakpoints) and checks that
      the new bounds enclose every sample;
  (2) counts how often the first version's real-segment bounds, |nh_P_seg(nh, acb(rho_1 e^{2t}))|, fail at
      the same samples (the defect, on the actual contour);
  (3) checks the tail's disc bound on the ray t in [-40, -8];
  (4) re-runs the body of tests/test_part1.py::test_S1_closing_P_bounds_at_complex_rho with nh_P_bounds
      replaced by the first version's real-segment bound, to see whether the new test would have caught S1;
  (5) re-runs certify_closing on both sides and compares every field with the stored closing leaves.
Usage: s1_recheck.py OUT.json
"""
from __future__ import annotations

import json
import math
import sys
from fractions import Fraction
from pathlib import Path

W = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(W / "src"))

from flint import acb, arb  # noqa: E402

from p2 import throat as TH  # noqa: E402
from p2.arbutil import git_commit, precision  # noqa: E402
from p2.m3 import fa  # noqa: E402
from p2.params import Point  # noqa: E402

CERT = W / "certificates" / "part1"
WB = Fraction(99, 12500)


def rho1_of(pt, nh):
    z1 = fa(TH.ZETA_1) * pt.OmegaP
    rho1 = arb(((z1 * z1) / nh["g1"]).upper()) * arb("1.0001")
    for _ in range(40):
        Pl = TH.nh_P(nh, acb(rho1)).real - TH.nh_P_err(nh, rho1)
        if rho1 < nh["R"] and rho1 * Pl >= z1 * z1:
            return rho1
        rho1 = rho1 * arb("1.05")
    raise ArithmeticError("rho_1")


def old_bounds(nh, rho1, civ):
    seg = abs(TH.nh_P_seg(nh, acb(rho1 * (2 * civ).exp())))
    return arb(seg.lower()), arb(seg.upper())


def side_check(pt, side, x):
    D = TH.SDesign.from_key(x["design"])
    nh = TH.nh_model(pt, side, TH.D_R[side])
    rho1 = rho1_of(pt, nh)
    nodes = sorted({TH.T_LEFT, D.ta} | {v for v in D.breakpoints() if TH.T_LEFT < v < D.ta}
                   | {m * TH.H_T for m in range(math.ceil(TH.T_LEFT / TH.H_T), math.floor(D.ta / TH.H_T) + 1)
                      if TH.T_LEFT < m * TH.H_T < D.ta})
    n_samp, n_new_fail, n_old_lo_fail, n_old_hi_fail = 0, 0, 0, 0
    worst_old_lo = 0.0
    for a, b in zip(nodes[:-1], nodes[1:]):
        pb = D.piece_bend(a, b)
        civ = (fa(a) + fa(b)) / 2 + arb(0, 1) * (fa(b - a) / 2)
        rabs = arb((rho1 * (2 * civ).exp()).upper())
        z = acb(rho1) * (2 * acb(civ, D.eta_value(civ, pb))).exp()
        Plo, Phi = TH.nh_P_bounds(nh, z, rabs)
        olo, ohi = old_bounds(nh, rho1, civ)
        for frac in (Fraction(0), Fraction(1, 4), Fraction(1, 2), Fraction(3, 4), Fraction(1)):
            t = fa(a + (b - a) * frac)
            for lam in ("1", "0.5", "0.1", "0.01", "0"):
                v = abs(TH.nh_P(nh, arb(lam) * rho1 * (2 * acb(t, D.eta_value(t, pb))).exp()))
                n_samp += 1
                if not (v.lower() >= Plo and v.upper() <= Phi):
                    n_new_fail += 1
                if v.upper() < olo:
                    n_old_lo_fail += 1
                    worst_old_lo = max(worst_old_lo, float(((olo - v) / olo).upper()))
                if v.lower() > ohi:
                    n_old_hi_fail += 1
    # tail: t in [-40, -8] on the ray (eta = t), |rho| <= rho_1 e^{2 t_L}
    rL = arb((rho1 * (2 * fa(TH.T_LEFT)).exp()).upper())
    _, hiL = TH.nh_P_disc(nh, rL)
    n_tail, n_tail_fail = 0, 0
    for k in range(0, 321):
        t = fa(TH.T_LEFT) - fa(Fraction(k, 10))
        for lam in ("1", "0.3"):
            v = abs(TH.nh_P(nh, arb(lam) * rho1 * (2 * acb(t, t)).exp()))
            n_tail += 1
            if not (v.upper() <= hiL):
                n_tail_fail += 1
    return {"side": "plus" if side > 0 else "c", "cells": len(nodes) - 1, "samples": n_samp,
            "new_bound_violations": n_new_fail, "old_lower_bound_violations": n_old_lo_fail,
            "old_upper_bound_violations": n_old_hi_fail, "old_lower_worst_relative_excess": worst_old_lo,
            "tail_samples": n_tail, "tail_disc_violations": n_tail_fail}


def test_with_old_bound(pt):
    """The body of test_S1_closing_P_bounds_at_complex_rho with nh_P_bounds replaced by the old bound."""
    D = TH.SDesign(2, Fraction(1, 20), 3, 0, 1)
    cells = ((Fraction(-2), Fraction(-31, 16)), (Fraction(-13, 8), Fraction(-25, 16)), (Fraction(-1), Fraction(-15, 16)),
             (Fraction(-1, 4), Fraction(-3, 16)), (Fraction(1), Fraction(17, 16)), (Fraction(2), Fraction(41, 20)))
    fails = []
    for side, rho1 in ((1, arb("0.00031")), (-1, arb("0.0034"))):
        nh = TH.nh_model(pt, side, TH.D_R[side])
        for a, b in cells:
            pb = D.piece_bend(a, b)
            civ = (fa(a) + fa(b)) / 2 + arb(0, 1) * (fa(b - a) / 2)
            Plo, Phi = old_bounds(nh, rho1, civ)
            for frac in (Fraction(0), Fraction(1, 2), Fraction(1)):
                t = fa(a + (b - a) * frac)
                for lam in ("1", "0.5", "0.1"):
                    zt = arb(lam) * rho1 * (2 * acb(t, D.eta_value(t, pb))).exp()
                    v = abs(TH.nh_P(nh, zt))
                    if not (v.lower() >= Plo and v.upper() <= Phi):
                        fails.append([side, str(a), str(b), lam, str(frac)])
    return {"assertion_failures_with_old_bound": len(fails), "first": fails[:4]}


def main():
    pt = Point("PB", 256)
    out = {"commit": git_commit()}
    g = TH.FGeo(pt)
    with precision(TH.PREC):
        out["sides"] = []
        out["regenerated_closing_leaves"] = []
        for side, tag in ((1, "plus"), (-1, "c")):
            cert = json.loads((CERT / f"S_PB_{tag}.json").read_text())
            x = [y for y in cert["leaves"] if y.get("closing")][0]
            r = side_check(pt, side, x)
            out["sides"].append(r)
            print(json.dumps(r), flush=True)
            new = json.loads(json.dumps(TH.certify_closing(pt, g, side, WB, design={"design": x["design"], "target": x["target"]})))
            diffs = [k for k in ("design", "near_horizon", "exit", "w1_check", "target", "mu_cells", "mu_tail", "worst_cell",
                                 "cells", "n_fails", "pass", "eps_max") if new.get(k) != x.get(k)]
            out["regenerated_closing_leaves"].append({"side": tag, "stored_mu_cells": x["mu_cells"], "stored_mu_tail": x["mu_tail"],
                                                      "recomputed_mu_cells": new["mu_cells"], "diffs": diffs})
            print(tag, x["mu_cells"], new["mu_cells"], diffs, flush=True)
        out["test_body_with_old_bound"] = test_with_old_bound(pt)
        print(out["test_body_with_old_bound"], flush=True)
    Path(sys.argv[1] if len(sys.argv) > 1 else "s1_recheck.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
