"""Independent check of the L certificates (S4 review).  No p2 code.

For sampled L leaves (near Re w = 0, near Re w = 8 kappa_-, near Im w = -d_delta, and random), the Jost solutions
R_+ and R_c are computed at the stored exits r_+ + h^2 tau_a and r_c - h^2 tau_a^c for actual integer l (and a very
large l as a proxy for the h -> 0 limit) by a Frobenius series in r itself (not in the scaled variable tau),
mpmath at 30 digits, with the recurrence of Delta^2 R'' + Delta Delta' R' + (w^2 r^4 - L Delta) R = 0 derived here
independently.  y = Phi_x/Phi = F/r + F R'/R.  Every value must satisfy the stored bound: Re y_+ >= the stored lower
bound (side +), Re y_c <= the stored upper bound (side c).  The exit rule h^-2 U >= max(4|w|^2_max, C_W) is
re-checked on a grid of h^2 in [0, h0^2].
Usage: indep_L.py OUT.json        52 boxes, 6 points each, l in {20, 21, 33, 100, 1000, 10^5}
       indep_L.py OUT.json all    all 1760 boxes, corners and centre, l in {20, 21, 100, 10^5}
"""
from __future__ import annotations

import json
import random
import sys
from fractions import Fraction
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from indep_params import PB  # noqa: E402

CERT = Path(__file__).resolve().parents[1] / "part1"
mp.mp.dps = 30
p = PB(40)
H0 = mp.mpf(2) / 41
CW = 4 * max(p.kp ** 2, mp.mpf(99) / 12500)


def delta_shift(rh):
    """Coefficients d_k of Delta(rh + t) = sum d_k t^k (k = 0..4), d_0 = 0 at a root."""
    c4 = -p.Lam / 3
    # Delta(r) = c4 r^4 + r^2 - 2 r + Q2
    d0 = p.Delta(rh)
    d1 = 4 * c4 * rh ** 3 + 2 * rh - 2
    d2 = 6 * c4 * rh ** 2 + 1
    d3 = 4 * c4 * rh
    d4 = c4
    return [mp.mpf(0), d1, d2, d3, d4], d0


def pmul(a, b):
    out = [mp.mpc(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def jost_y(side, l, w, t_exit):
    """y = Phi_x/Phi of the Jost solution at r = r_h + t_exit (t_exit > 0 at r_+, < 0 at r_c)."""
    rh = p.rp if side > 0 else p.rc
    kap = p.kp if side > 0 else p.kc
    d, _ = delta_shift(rh)
    dp = [d[k + 1] * (k + 1) for k in range(4)]
    A = pmul(d, d)                 # Delta^2
    B = pmul(d, dp)                # Delta Delta'
    rr = [rh, mp.mpf(1)]
    r2 = pmul(rr, rr)
    r4 = pmul(r2, r2)
    L = l * (l + 1)
    C = [w * w * r4[k] - L * d[k] for k in range(5)]
    rho = -1j * w / (2 * kap)
    # sum_k a_k [A_{n-k+2}(k+rho)(k+rho-1) + B_{n-k+1}(k+rho) + C_{n-k}] = 0, the k = n term is d1^2 n (n + 2 rho)
    a = [mp.mpc(1)]
    S, Sd = mp.mpc(1), rho / t_exit
    tn = mp.mpf(1)
    n = 0
    small = 0
    while True:
        n += 1
        acc = mp.mpc(0)
        for k in range(max(0, n - 6), n):
            j2, j1, j0 = n - k + 2, n - k + 1, n - k
            term = mp.mpc(0)
            if j2 < len(A):
                term += A[j2] * (k + rho) * (k + rho - 1)
            if j1 < len(B):
                term += B[j1] * (k + rho)
            if j0 < len(C):
                term += C[j0]
            acc += a[k] * term
        an = -acc / (d[1] ** 2 * n * (n + 2 * rho))
        a.append(an)
        tn *= t_exit
        S += an * tn
        Sd += an * (n + rho) * tn / t_exit
        if abs(an * tn) < mp.mpf(10) ** (-mp.mp.dps) * abs(S) and n > 20:
            small += 1
            if small > 10:
                break
        else:
            small = 0
        if n > 20000:
            raise RuntimeError("series too long")
    # R = |t|^rho S, dR/dr = |t|^rho Sd  (d/dr |t|^rho = rho |t|^rho / t for either sign of t)
    r = rh + t_exit
    F = p.F(r)
    return F / r + F * Sd / S


def main():
    out = {"rows": []}
    certp = json.loads((CERT / "L_PB_plus.json").read_text())
    certc = json.loads((CERT / "L_PB_c.json").read_text())
    lp, lc = certp["leaves_detail"], certc["leaves_detail"]
    assert [x["box"] for x in lp] == [x["box"] for x in lc]
    rng = random.Random(20261019)
    picks = []
    ALL = len(sys.argv) > 2 and sys.argv[2] == "all"
    for i, x in enumerate(lp):
        b = [Fraction(v) for v in x["box"]]
        if b[0] == 0 or b[1] == 8 or b[2] == Fraction(-11, 20):
            if rng.random() < 0.25:
                picks.append(i)
    picks += rng.sample(range(len(lp)), 12)
    picks = sorted(set(picks))
    ls = [20, 21, 33, 100, 1000, 100000]
    if ALL:
        picks = list(range(len(lp)))
        ls = [20, 21, 100, 100000]
    nbad = 0
    worst_margin_p, worst_margin_c = None, None
    for i in picks:
        bp, bc = lp[i], lc[i]
        b = [Fraction(v) for v in bp["box"]]
        ylo = mp.mpf(bp["Re_y_bound"].split("+/-")[0].strip("[ ")) - mp.mpf(bp["Re_y_bound"].split("+/-")[1].strip("] "))
        yhi = mp.mpf(bc["Re_y_bound"].split("+/-")[0].strip("[ ")) + mp.mpf(bc["Re_y_bound"].split("+/-")[1].strip("] "))
        tau_p, tau_c = Fraction(bp["tau_a"]), Fraction(bc["tau_a"])
        pts = [(b[0], b[2]), (b[1], b[2]), (b[0], b[3]), (b[1], b[3]), ((b[0] + b[1]) / 2, (b[2] + b[3]) / 2)]
        if not ALL:
            pts.append((b[0] + (b[1] - b[0]) * Fraction(rng.random()), b[2] + (b[3] - b[2]) * Fraction(rng.random())))
        # exit rule on a grid of h^2 in [0, h0^2]
        wmax2 = max(abs(mp.mpc(float(re), float(im))) ** 2 for re in (b[0], b[1]) for im in (b[2], b[3])) * p.km ** 2
        M = max(4 * wmax2, CW)
        ex_ok = True
        for j in range(0, 201):
            h2 = H0 ** 2 * j / 200
            for side, tau in ((1, tau_p), (-1, tau_c)):
                rh = p.rp if side > 0 else p.rc
                if h2 == 0:
                    d, _ = delta_shift(rh)
                    val = d[1] * side * mp.mpf(tau.numerator) / tau.denominator / rh ** 4      # D(tau)/r^4 at h = 0
                else:
                    r = rh + side * h2 * mp.mpf(tau.numerator) / tau.denominator
                    val = p.U(r) / h2
                ex_ok = ex_ok and val >= M
        row = {"box": bp["box"], "tau_plus": str(tau_p), "tau_c": str(tau_c), "stored_lo_plus": float(ylo),
               "stored_hi_c": float(yhi), "exit_rule_ok_grid": bool(ex_ok), "values": []}
        for (re_, im_) in pts:
            w = mp.mpc(mp.mpf(re_.numerator) / re_.denominator, mp.mpf(im_.numerator) / im_.denominator) * p.km
            for l in ls:
                h2 = 1 / (mp.mpf(l) + mp.mpf(1) / 2) ** 2
                yp = jost_y(1, l, w, h2 * mp.mpf(tau_p.numerator) / tau_p.denominator)
                yc = jost_y(-1, l, w, -h2 * mp.mpf(tau_c.numerator) / tau_c.denominator)
                okp = mp.re(yp) >= ylo
                okc = mp.re(yc) <= yhi
                if not (okp and okc):
                    nbad += 1
                mp_ = float(mp.re(yp) - ylo)
                mc_ = float(yhi - mp.re(yc))
                worst_margin_p = mp_ if worst_margin_p is None else min(worst_margin_p, mp_)
                worst_margin_c = mc_ if worst_margin_c is None else min(worst_margin_c, mc_)
                row["values"].append({"w_over_km": [float(re_), float(im_)], "l": l, "Re_y_plus": float(mp.re(yp)),
                                      "Re_y_c": float(mp.re(yc)), "ok": bool(okp and okc)})
        out["rows"].append(row)
        if ALL and i % 100:
            continue
        print(bp["box"], "exit ok", ex_ok, "min Re y+ - bound", min(v["Re_y_plus"] for v in row["values"]) - float(ylo),
              "bound - max Re y_c", float(yhi) - max(v["Re_y_c"] for v in row["values"]), flush=True)
    out["n_boxes"] = len(picks)
    out["n_values"] = sum(len(r["values"]) for r in out["rows"])
    out["n_violations"] = nbad
    out["all_exit_rules_ok"] = all(r["exit_rule_ok_grid"] for r in out["rows"])
    out["min_margin_plus"] = worst_margin_p
    out["min_margin_c"] = worst_margin_c
    print({k: v for k, v in out.items() if k != "rows"})
    Path(sys.argv[1] if len(sys.argv) > 1 else "indep_L.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
