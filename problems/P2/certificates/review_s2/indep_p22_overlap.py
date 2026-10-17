"""P2-2 overlap check for the S2 windows (reviewer script, 2026-10-17).

1. a_l stored in the m3 certificates = W_l_exact of the S1 lg certificates (l = 20, 50) and of
   certificates/m3/lg_PB_l100.json (l = 100).
2. Certified (Arb, params.json balls, reviewer's own window): a_l < (l+1/2) Omega_P + 6 gamma0, and the
   margin in kappa_- units.
3. Replay of the S1 criterion (*) for l = 100 at the stored (a_100, r_X) with the S1 code (lg_bounds.LG).
4. Independent float check of what (*) bounds, at l = 100: |W e^{i Theta}/(2 i w) - 1| for
   w = a_100 + i b (b = 0, -beta/2, -beta) and larger Re w, from the reviewer's own ODE integration
   along the real tortoise axis (Jost data at x = -200 and x = +120, DOP853), with
   Theta = int (sqrt(w^2 - V) - w) dx.  The certified bound is criterion - 1 = 0.881.
Run (repo root):  nice -n 10 .venv/bin/python problems/P2/certificates/review_s2/indep_p22_overlap.py
"""
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np
from flint import arb, ctx, fmpq
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
P2 = HERE.parents[1]
sys.path.insert(0, str(P2 / "src"))
ctx.prec = 128
T0 = time.time()
J = json.loads((P2 / "certificates" / "params.json").read_text())["points"]["PB"]["quantities"]
OMP, G0, KM = (arb(J[k]["ball"]) for k in ("Omega_P", "gamma0", "kappa_minus"))
out = {"overlap": {}}

# 1-2: stored thresholds and certified comparison
for l in (20, 50, 100):
    m3c = json.loads((P2 / "certificates" / "m3" / f"m3_PB_l{l}.json").read_text())
    src = P2 / "certificates" / ("lg" if l <= 50 else "m3") / f"lg_PB_l{l}.json"
    lgc = json.loads(src.read_text())
    a = Fraction(lgc["W_l_exact"])
    al = arb(fmpq(a.numerator, a.denominator))
    top = (l + arb(1) / 2) * OMP + 6 * G0
    out["overlap"][f"l{l}"] = {"source": str(src.relative_to(P2.parents[1])), "a_l_exact": lgc["W_l_exact"],
                               "same_as_m3_certificate": m3c["P2_2_overlap"]["a_l_exact"] == lgc["W_l_exact"],
                               "lg_criterion_certified_flag": lgc["criterion"]["certified"],
                               "lg_criterion_value_upper": lgc["criterion"]["value_upper"],
                               "a_l_over_kappa_minus": (al / KM).str(10), "top_over_kappa_minus": (top / KM).str(10),
                               "a_l_lt_top_certified": bool(al < top), "margin_over_kappa_minus": ((top - al) / KM).str(6)}
    print(l, out["overlap"][f"l{l}"], flush=True)

# 3: replay the S1 criterion at l = 100 (S1 code path)
from p2 import lg_bounds  # noqa: E402
from p2.params import Point  # noqa: E402

st = json.loads((P2 / "certificates" / "m3" / "lg_PB_l100.json").read_text())
lg = lg_bounds.LG(Point("PB"), 100)
a0 = fmpq(*[int(x) for x in st["W_l_exact"].split("/")])
rX = fmpq(*[int(x) for x in st["matching_point_exact"].split("/")])
res = lg.criterion(a0, rX)
out["l100_criterion_replay"] = {"value": res["value"].str(12), "ok": res["ok"], "stored": st["criterion"]["value_upper"]}
print("criterion replay l=100:", out["l100_criterion_replay"], flush=True)

# 4: independent float Wronskian at l = 100
mp.mp.dps = 40
g = lambda k: float(J[k]["approx"])
rn, rm, rp, rc, Q2 = g("r_n"), g("r_minus"), g("r_plus"), g("r_c"), g("Q2")
Amp = [mp.mpf(J["A_coefficients"][k]["approx"]) for k in ("r_n", "r_minus", "r_plus", "r_c")]
Rmp = [mp.mpf(J[k]["approx"]) for k in ("r_n", "r_minus", "r_plus", "r_c")]
beta = float(J["strip_depth_beta"]["approx"])
km = g("kappa_minus")
L = 100 * 101


def rstar(r):
    return sum(A * mp.log(abs(r - ri)) for A, ri in zip(Amp, Rmp))


def r_of_x(x):
    lo, hi = Rmp[2] + mp.mpf("1e-30"), Rmp[3] - mp.mpf("1e-30")
    for _ in range(400):
        mid = (lo + hi) / 2
        if rstar(mid) > x:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def Fr(r):
    return -(r - rn) * (r - rm) * (r - rp) * (r - rc) / (50 * r * r)


def Vr(r):
    D = -r ** 4 / 50 + r * r - 2 * r + Q2
    Dp = -4 * r ** 3 / 50 + 2 * r - 2
    return D / r ** 2 * (L / r ** 2 + (r * Dp - 2 * D) / r ** 4)


def integrate(w, x_from, x_to, sign):
    """Jost solution e^{sign i w x} at x_from, carried to x_to; also int (p - w) dx over the span."""
    r0 = float(r_of_x(x_from))
    ph = np.exp(sign * 1j * w * x_from)

    def f(x, y):
        r = y[0]
        Phi = y[1] + 1j * y[2]
        dPhi = y[3] + 1j * y[4]
        v = Vr(r)
        dd = (v - w * w) * Phi
        p = np.sqrt(w * w - v)
        return [Fr(r), dPhi.real, dPhi.imag, dd.real, dd.imag, (p - w).real, (p - w).imag]
    y0 = [r0, ph.real, ph.imag, (sign * 1j * w * ph).real, (sign * 1j * w * ph).imag, 0.0, 0.0]
    sol = solve_ivp(f, (x_from, x_to), y0, method="DOP853", rtol=1e-12, atol=1e-14)
    y = sol.y[:, -1]
    return y[1] + 1j * y[2], y[3] + 1j * y[4], y[5] + 1j * y[6], y[0]


xm = float(rstar(mp.mpf(int(rX.p)) / int(rX.q)))
a100 = float(Fraction(st["W_l_exact"]))
rows = []
for a in (a100, 1.01 * a100, 1.05 * a100, 1.2 * a100):
    for b in (0.0, -beta / 2, -beta):
        w = complex(a, b)
        P, dP, thp, _ = integrate(w, -200.0, xm, -1)
        C, dC, thc, _ = integrate(w, 120.0, xm, +1)
        Wr = P * dC - dP * C
        Theta = thp - thc                       # int_{-200}^{xm} + int_{xm}^{120}
        dev = abs(Wr * np.exp(1j * Theta) / (2j * w) - 1)
        rows.append({"Re_w_over_kappa_minus": a / km, "Im_w_over_kappa_minus": b / km, "deviation": float(dev)})
        print(rows[-1], flush=True)
out["l100_float_wronskian"] = {"rows": rows, "max_deviation": max(r["deviation"] for r in rows),
                               "certified_bound": "criterion - 1 = 0.881"}
out["seconds"] = round(time.time() - T0, 1)
(HERE / "indep_p22_overlap.json").write_text(json.dumps(out, indent=1) + "\n")
print("done", out["seconds"], "s")
