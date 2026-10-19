"""Independent check of the S closing box's near-horizon model (S4 review).  No p2 code.

At each horizon r_h (side + : r_+, d_r = 3/50; side c : r_c, d_r = 1), with s = side (r - r_h) and
rho(s) = s f(s) = exp(+-2 kappa_h (x - c_h)):
  * f, |f| on the box |Re s|, |Im s| <= d_r (dense sampling), against the stored lower bound m_f; |s f(s)| on the
    circle |s| = d_r against R = m_f d_r (the Rouche condition);
  * sup |U| and sup |W~| on the box (dense sampling) against the stored M_G, M_W;
  * the Taylor coefficients g_n of G(rho) = U(r_h + side psi(rho)) by a Cauchy integral (FFT on |rho| = R/2,
    psi solved by Newton), against the stored g_1, g_2, and the Cauchy bound |g_n| <= M_G/R^n;
  * rho_1: G(rho_1) >= zeta_1^2, and rho_1 e^{2 t_a} < R.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from indep_params import PB  # noqa: E402

CERT = Path(__file__).resolve().parents[1] / "part1"
p = PB()


def num(s):
    return float(mp.mpf(s.split("+/-")[0].strip("[ ")))


def side_check(side):
    mp.mp.dps = 30
    i = 2 if side > 0 else 3
    rh = p.roots[i]
    kap = p.kp if side > 0 else p.kc
    dr = mp.mpf(3) / 50 if side > 0 else mp.mpf(1)
    others = [(a, r) for j, (a, r) in enumerate(zip(p.A, p.roots)) if j != i]

    def f(s):
        Lsum = sum(a * mp.log(1 + side * s / (rh - r)) for a, r in others)
        return mp.exp(side * 2 * kap * Lsum)

    def U(s):
        r = rh + side * s
        return p.Delta(r) / r ** 4

    def Wt(s):
        return p.Wt(rh + side * s)

    cert = json.loads((CERT / f"S_PB_{'plus' if side > 0 else 'c'}.json").read_text())
    x = [y for y in cert["leaves"] if y.get("closing")][0]
    nh = x["near_horizon"]
    # box sampling
    n = 121
    g = np.linspace(-1, 1, n)
    fmin, Umax, Wmax = None, 0, 0
    for a in g:
        for b in g:
            s = mp.mpc(a * dr, b * dr)
            fv = abs(f(s))
            fmin = fv if fmin is None or fv < fmin else fmin
            Umax = max(Umax, abs(U(s)))
            Wmax = max(Wmax, abs(Wt(s)))
    circ = min(abs(dr * mp.expjpi(2 * k / 720) * f(dr * mp.expjpi(2 * k / 720))) for k in range(720))
    m_f, R, M_G, M_W = num(nh["m_f"]), num(nh["R"]), num(nh["M_G"]), num(nh["M_W"])
    # Taylor coefficients of G by FFT on |rho| = R/2
    N = 64
    rad = mp.mpf(R) / 2
    vals = []
    s_guess = rad / f(0)
    for k in range(N):
        rho = rad * mp.expjpi(2 * mp.mpf(k) / N)
        s = s_guess if k == 0 else s
        for _ in range(60):
            fs = f(s)
            # d/ds (s f(s)) = f + s f' ; f'/f = side 2 kap sum a * side/(rh - r) / (1 + side s/(rh - r))
            dl = side * 2 * kap * sum(a * side / (rh - r) / (1 + side * s / (rh - r)) for a, r in others)
            g_ = s * fs - rho
            dg = fs * (1 + s * dl)
            ds = g_ / dg
            s = s - ds
            if abs(ds) < mp.mpf(10) ** (-28):
                break
        vals.append(U(s))
    coeffs = []
    for nn in range(16):
        c = sum(vals[k] * mp.expjpi(-2 * mp.mpf(k) * nn / N) for k in range(N)) / N / rad ** nn
        coeffs.append(c)
    g1s, g2s = num(nh["g1"]), num(nh["g2"])
    cauchy_ok = all(abs(coeffs[k]) <= M_G / R ** k * 1.000001 for k in range(1, 16))
    # rho_1
    rho1 = num(nh["rho_1"])
    z1 = mp.mpf(1) / 50 * p.OmegaP
    s = rho1 / f(0)
    for _ in range(60):
        fs = f(s)
        dl = side * 2 * kap * sum(a * side / (rh - r) / (1 + side * s / (rh - r)) for a, r in others)
        ds = (s * fs - rho1) / (fs * (1 + s * dl))
        s = s - ds
    G1 = U(s)
    ta = 2.05
    return {"side": "plus" if side > 0 else "c",
            "m_f_stored": m_f, "min_abs_f_sampled": float(fmin), "m_f_ok": bool(m_f <= fmin),
            "R_stored": R, "min_abs_sf_on_circle": float(circ), "rouche_ok": bool(circ >= R),
            "M_G_stored": M_G, "max_abs_U_sampled": float(Umax), "M_G_ok": bool(Umax <= M_G),
            "M_W_stored": M_W, "max_abs_Wt_sampled": float(Wmax), "M_W_ok": bool(Wmax <= M_W),
            "g1_stored": g1s, "g1_fft": float(mp.re(coeffs[1])), "g2_stored": g2s, "g2_fft": float(mp.re(coeffs[2])),
            "g0_fft_abs": float(abs(coeffs[0])), "cauchy_bounds_ok": bool(cauchy_ok),
            "rho_1": rho1, "G_rho1_over_zeta1sq": float(mp.re(G1) / z1 ** 2), "rho1_ok": bool(mp.re(G1) >= z1 ** 2),
            "rho_max_over_R": float(rho1 * mp.e ** (2 * ta) / R)}


if __name__ == "__main__":
    out = [side_check(1), side_check(-1)]
    print(json.dumps(out, indent=1))
    Path(sys.argv[1] if len(sys.argv) > 1 else "nh_check.json").write_text(json.dumps(out, indent=1) + "\n")
