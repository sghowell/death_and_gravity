"""Sharper negative/positive control of the author's rigorous checker (reviewer script, 2026-10-17).

The committed A2 control (a) uses the window's ground state, designed at Im w = -0.55 kappa_-, on a box
around the PS mode at Im w ~ -1.68 kappa_-; (b) the box-specific design overflows and never reaches the
checker.  Here the ground state is designed close to the mode (Im w = -0.88 gamma0, where the float
margin is still positive) and the SAME checker (m3.certify_l with an explicit box) is run on:
  * IN:     a box of half-width gamma0/50 around the l-th PS fundamental mode (contains a QNM): must fail;
  * STRAD:  Re as IN, Im in [-1.05, -0.95] gamma0 (contains the mode): must fail;
  * OUT:    Re as IN, Im in [-0.90, -0.86] gamma0 (no QNM; the mode is ~0.1 gamma0 below): may certify.
A checker that certified IN or STRAD would be unsound; OUT certifying shows the checker resolves the
mode location to ~0.1 gamma0, i.e. the A2 failure is not merely a design artefact.
Run (repo root):  nice -n 10 .venv/bin/python problems/P2/certificates/review_s2/sharp_negative_control.py [l ...]
(default l = 100 50; the committed output used 100 50 20, 416 s of CPU)
"""
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from flint import arb  # noqa: E402

from p2 import m3  # noqa: E402
from p2.arbutil import precision  # noqa: E402
from p2.params import Point  # noqa: E402

PS = {20: (100.823703512, -1.68357264593), 50: (248.396138233, -1.6832715886), 100: (494.340406988, -1.68322053463)}
ls = [int(x) for x in sys.argv[1:]] or [100, 50]
out = {}
pt = Point("PB", 256)
with precision(m3.PREC):
    geo, rays, tails = m3.geometry_certificate(pt)
    km, g0 = float(pt.km.mid()), float(pt.gamma0.mid())
    for l in ls:
        t0 = time.time()
        re_ps, im_ps = PS[l][0] * km, PS[l][1] * km
        e = g0 / 50
        wd = complex(re_ps, -0.88 * g0)
        lam, nodes = m3.design_ground_state(pt, l, wd)
        muHO = float((2 * (l + arb(1) / 2) * pt.OmegaP * pt.gamma0).mid())
        boxes = {"IN (contains PS mode)": (re_ps - e, re_ps + e, im_ps - e, im_ps + e),
                 "STRAD Im in [-1.05,-0.95]g0 (contains mode)": (re_ps - e, re_ps + e, -1.05 * g0, -0.95 * g0),
                 "OUT Im in [-0.90,-0.86]g0 (no QNM)": (re_ps - e, re_ps + e, -0.90 * g0, -0.86 * g0)}
        row = {"design_point_Im_over_gamma0": -0.88, "float_lambda_at_design_over_muHO": lam / muHO,
               "PS_mode_Im_over_gamma0": im_ps / g0}
        for name, (a1, a2, b1, b2) in boxes.items():
            c = m3.certify_l(pt, l, rays, tails, nodes=nodes, lam_float=lam,
                             box=(arb(a1), arb(a2), arb(b1), arb(b2)))
            mu = float(arb(c["mu_cert"].split()[0].strip("[")).mid())
            row[name] = {"certified": c["A1_pass"], "mu_cert_over_muHO": mu / muHO, "cells": c["cells_evaluated"]}
        row["expected_pattern_holds"] = (not row["IN (contains PS mode)"]["certified"]
                                         and not row["STRAD Im in [-1.05,-0.95]g0 (contains mode)"]["certified"])
        row["seconds"] = round(time.time() - t0, 1)
        out[f"l{l}"] = row
        print(l, json.dumps(row), flush=True)
(HERE / "sharp_negative_control.json").write_text(json.dumps(out, indent=1) + "\n")
