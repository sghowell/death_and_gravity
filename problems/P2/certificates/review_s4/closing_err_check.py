"""The S closing box: the truncation-error term (S4 review, finding S1 on throat.certify_closing.err_cell).

The defect described below is that of the code at 0e329217; it was fixed at 68de8184 (see s1_recheck.py).
This script re-implements both versions of the error term itself, so it gives the same numbers at either commit.

throat.certify_closing bounds |P_N(rho)| (rho = rho_T e^{2Y}, complex) through nh_P_seg evaluated at the REAL point
rho_1 e^{2t} (err_cell: `PN = nh_P_seg(nh, acb(rho1 * e2t))`), and likewise |P_N| on the tail (UL).  That is not an
enclosure of P_N on the complex set the contour visits.  This script re-runs the closing-box cell check with the
author's code, unchanged except for err_cell and the tail term, where |P_N| is bounded rigorously by
  (a) the disc bounds |g_1| -/+ sum_{n>=2} |g_n| |rho|^{n-1}  (valid for every complex |rho| <= rabs), and
  (b) the complex segment enclosure nh_P_seg(nh, rho_1 * e^{2Y(cell)}),
taking for each cell the better (larger) lower bound and the smaller upper bound of the two.
It reports the stored and corrected mu_cells and mu_tail, the largest change of a cell's lower bound, and the
largest ratio between the corrected and the original error term.
"""
from __future__ import annotations

import itertools
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

W = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(W / "src"))

from flint import acb, acb_series, arb  # noqa: E402

from p2 import throat as TH  # noqa: E402
from p2.arbutil import precision  # noqa: E402
from p2.m3 import fa  # noqa: E402
from p2.params import Point  # noqa: E402


def disc_bounds(nh, rabs):
    g = nh["g"]
    tail = arb(0)
    zp = arb(1)
    for n in range(2, nh["N"] + 1):
        zp = zp * rabs
        tail = tail + abs(g[n]).upper() * zp
    g1 = abs(g[1])
    return arb((g1 - tail).lower()), arb((g1 + tail).upper())


def run(side, design, stored):
    pt = Point("PB", 256)
    W_bound = Fraction(99, 12500)
    with precision(TH.PREC):
        D = TH.SDesign.from_key(design["design"])
        kap = pt.kp if side > 0 else pt.kc
        nh = TH.nh_model(pt, side, TH.D_R[side])
        z1 = fa(TH.ZETA_1) * pt.OmegaP
        g1 = nh["g1"]
        rho1 = arb(((z1 * z1) / g1).upper()) * arb("1.0001")
        for _ in range(40):
            Pl = TH.nh_P(nh, acb(rho1)).real - TH.nh_P_err(nh, rho1)
            if rho1 < nh["R"] and rho1 * Pl >= z1 * z1:
                break
            rho1 = rho1 * arb("1.05")
        eps_max = TH.eps_max_of(pt, Fraction(0))
        errT = TH.nh_P_err(nh, rho1)
        PT = TH.nh_P_seg(nh, acb(rho1))
        PT_abs_lo = arb((abs(PT).lower() - errT).lower())
        prm = TH.Params(pt, None, None, Fraction(0), eps_max)
        k2 = kap * kap

        def model(tt, pb, pp):
            e1, _, _ = D.eta_derivs(tt, pb)
            eta = D.eta_value(tt, pb)
            e2Y = (2 * acb_series([acb(tt, eta), acb(1, e1)], prec=2)).exp()
            Ut = e2Y * TH.nh_P_seg(nh, acb(rho1) * e2Y) / PT
            return TH.parts(pt, D, tt, None, None, side, kap, pb, pp, model=Ut), e2Y

        def err_orig(civ, p):
            e2t = (2 * civ).exp()
            rabs = arb((rho1 * e2t).upper())
            PN = TH.nh_P_seg(nh, acb(rho1 * e2t))
            e1 = TH.nh_P_err(nh, rabs) / arb((abs(PN).lower()).max(arb("1e-300")))
            e2 = errT / PT_abs_lo
            delta = (e1 + e2) / (1 - e2)
            Uabs = e2t * arb((abs(PN).upper() + TH.nh_P_err(nh, rabs)).upper()) / PT_abs_lo
            q = rabs / nh["R"]
            Wabs = nh["M_W"] * q / (1 - q)
            return p["Yp2"][0] * (Uabs * delta + prm.E2.upper() * Wabs)

        def err_fixed(civ, p, e2Y):
            e2t = (2 * civ).exp()
            rabs = arb((rho1 * e2t).upper())
            lo_d, hi_d = disc_bounds(nh, rabs)
            PNc = TH.nh_P_seg(nh, acb(rho1) * e2Y.coeffs()[0])     # P_N on {lam rho_1 e^{2Y(t)}: lam in [0,1], t in cell}
            lo_c, hi_c = arb(abs(PNc).lower()), arb(abs(PNc).upper())
            P_lo = lo_d.max(lo_c)
            P_hi = hi_d.min(hi_c)
            if not (P_lo > 0):
                return None
            e1 = TH.nh_P_err(nh, rabs) / P_lo
            e2 = errT / PT_abs_lo
            delta = (e1 + e2) / (1 - e2)
            Uabs = e2t * arb((P_hi + TH.nh_P_err(nh, rabs)).upper()) / PT_abs_lo
            q = rabs / nh["R"]
            Wabs = nh["M_W"] * q / (1 - q)
            return p["Yp2"][0] * (Uabs * delta + prm.E2.upper() * Wabs)

        target = fa(Fraction(design["target"]))
        nodes = sorted({TH.T_LEFT, D.ta} | {x for x in D.breakpoints() if TH.T_LEFT < x < D.ta}
                       | {m * TH.H_T for m in range(math.ceil(TH.T_LEFT / TH.H_T), math.floor(D.ta / TH.H_T) + 1)
                          if TH.T_LEFT < m * TH.H_T < D.ta})
        mu_o, mu_f = None, None
        worst_change, worst_ratio, n = 0.0, 0.0, 0
        for a, b in itertools.pairwise(nodes):
            stack = [(a, b, 0)]
            while stack:
                a_, b_, dpt = stack.pop()
                pb, pp = D.piece_bend(a_, b_), D.piece_phase(a_, b_)
                A_, B_ = fa(a_), fa(b_)
                w = B_ - A_
                tm = (A_ + B_) / 2
                civ = tm + arb(0, 1) * (w / 2)
                pm, _ = model(tm, pb, pp)
                pc, e2Y = model(civ, pb, pp)
                qm = prm.combine(pm, 0, k2)
                dq = prm.combine(pc, 1, k2)
                base = qm.lower() - abs(dq).upper() * (w / 2)
                eo = err_orig(civ, pc)
                ef = err_fixed(civ, pc, e2Y)
                low_o = arb((base - eo.upper()).lower())
                low_f = arb((base - ef.upper()).lower()) if ef is not None else arb(-1)
                n += 1
                # the stored partition: refine exactly as production did (on the original bound)
                if low_o < target and dpt < TH.MAX_DEPTH:
                    m = (a_ + b_) / 2
                    stack.append((m, b_, dpt + 1))
                    stack.append((a_, m, dpt + 1))
                    continue
                mu_o = low_o if mu_o is None or low_o < mu_o else mu_o
                mu_f = low_f if mu_f is None or low_f < mu_f else mu_f
                worst_change = max(worst_change, float((low_o - low_f).upper()))
                if eo.upper() > 0:
                    worst_ratio = max(worst_ratio, float((ef.upper() / eo.upper()).upper()) if ef is not None else float("inf"))
        # the tail t <= t_L, original and fixed
        phinf = D.phinf()
        e2L = (2 * fa(TH.T_LEFT)).exp()
        rL = arb((rho1 * e2L).upper())
        UL_o = e2L * arb((abs(TH.nh_P_seg(nh, acb(rL))).upper() + TH.nh_P_err(nh, rL)).upper()) / PT_abs_lo
        _, hi_d = disc_bounds(nh, rL)
        UL_f = e2L * arb((hi_d + TH.nh_P_err(nh, rL)).upper()) / PT_abs_lo
        WL = nh["M_W"] * (rL / nh["R"]) / (1 - rL / nh["R"])
        base_t = (-2 * phinf.sin() + prm.EN * (-4 * phinf.cos()) + prm.E2N2 * (2 * phinf.sin()) - 2 * WL * prm.E2.upper())
        tail_o, tail_f = base_t - 2 * UL_o, base_t - 2 * UL_f
        return {"side": "plus" if side > 0 else "c", "cells": n, "stored_mu_cells": stored["mu_cells"],
                "mu_cells_original_err": TH.lohi(mu_o), "mu_cells_fixed_err": TH.lohi(mu_f),
                "largest_drop_of_a_cell_bound": worst_change, "largest_err_ratio_fixed_over_original": worst_ratio,
                "stored_mu_tail": stored["mu_tail"], "mu_tail_original": TH.lohi(tail_o), "mu_tail_fixed": TH.lohi(tail_f),
                "pass_fixed": bool(mu_f > 0 and tail_f > 0)}


def main():
    out = []
    for side, tag in ((1, "plus"), (-1, "c")):
        cert = json.loads((W / "certificates" / "part1" / f"S_PB_{tag}.json").read_text())
        x = [y for y in cert["leaves"] if y.get("closing")][0]
        r = run(side, {"design": x["design"], "target": x["target"]}, x)
        print(json.dumps(r, indent=1), flush=True)
        out.append(r)
    Path(sys.argv[1] if len(sys.argv) > 1 else "closing_err_check.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
