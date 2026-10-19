"""Point checks of W_l(w) != 0 in H+ for l in [20, 2000] by an independent float method (S4 review).  No p2 code.

Riccati in the tortoise coordinate x: dr/dx = F(r), dy/dx = V_l(r) - w^2 - y^2 (y = Phi_x/Phi), integrated in complex
double precision (scipy DOP853, rtol 1e-11) from r_h +- h^2 tau0 (tau0 = 1/4, start values from a Frobenius series in
r computed in mpmath, see indep_L.jost_y) to the photon sphere x_P from both sides.  W_l = 0 iff y_+(x_P) = y_c(x_P)
(Phi_+, Phi_c != 0 there), so the statistic is  D = |y_+ - y_c| / (|y_+| + |y_c|).

Positive control: at the l = 20 and l = 50 photon-sphere fundamental modes (below the strip; from the S1 float solver as
stored in controls_PB.json N1) D must be tiny, and a secant iteration on y_+ - y_c seeded nearby must converge to them.
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from indep_L import jost_y, p  # noqa: E402

CERT = Path(__file__).resolve().parents[1] / "part1"
lam3 = float(p.Lam / 3)
roots = [float(r) for r in p.roots]
rP = float(p.rP)
km, OmP, d = float(p.km), float(p.OmegaP), float(p.d)


def Ffun(r):
    return -lam3 * (r - roots[0]) * (r - roots[1]) * (r - roots[2]) * (r - roots[3]) / r ** 2


def Vfun(r, L):
    D = -lam3 * (r - roots[0]) * (r - roots[1]) * (r - roots[2]) * (r - roots[3])
    Dp = -4 * lam3 * r ** 3 + 2 * r - 2
    F = D / r ** 2
    Fr = Dp / r ** 2 - 2 * D / r ** 3
    return F * (L / r ** 2 + Fr / r)


def y_at_rP(side, l, w, tau0=0.25):
    h2 = 1 / (l + 0.5) ** 2
    t0 = side * h2 * tau0
    rh = p.rp if side > 0 else p.rc
    y0 = complex(jost_y(side, l, mp.mpc(w.real, w.imag), mp.mpf(t0)))
    r0 = float(rh + mp.mpf(t0))
    L = l * (l + 1)
    w2 = w * w

    def rhs(x, s):
        r = s[0]
        y = s[1] + 1j * s[2]
        dy = Vfun(r, L) - w2 - y * y
        return [Ffun(r), dy.real, dy.imag]

    def hit(x, s):
        return s[0] - rP
    hit.terminal = True
    x_span = (0.0, 1e4 * (1 if side > 0 else -1))
    sol = solve_ivp(rhs, x_span, [r0, y0.real, y0.imag], method="DOP853", rtol=1e-11, atol=1e-13, events=hit)
    assert sol.status == 1, sol.message
    s = sol.y_events[0][0]
    return s[1] + 1j * s[2], sol.nfev


def D_stat(l, w):
    yp, n1 = y_at_rP(1, l, w)
    yc, n2 = y_at_rP(-1, l, w)
    return abs(yp - yc) / (abs(yp) + abs(yc)), yp, yc, n1 + n2


def secant(l, w0, steps=30):
    def f(w):
        _, yp, yc, _ = D_stat(l, w)
        return yp - yc
    x0, x1 = w0, w0 * (1 + 1e-5) + 1e-6j
    f0, f1 = f(x0), f(x1)
    for _ in range(steps):
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, f0, x1, f1 = x1, f1, x2, f(x2)
        if abs(x1 - x0) < 1e-10 * abs(x1):
            break
    return x1


def main():
    out = {}
    ctrl = json.loads((CERT / "controls_PB.json").read_text())
    pc = []
    for l in ("20", "50"):
        m = ctrl["N1_N2"]["N1"][l]["PS_mode_over_kappa_minus"]
        wm = complex(float(m[0]), float(m[1])) * km
        Dm, _, _, _ = D_stat(int(l), wm)
        seed = wm + complex(0.3, 0.2) * km
        root = secant(int(l), seed)
        Dseed, _, _, _ = D_stat(int(l), seed)
        pc.append({"l": int(l), "mode_over_km": [float(m[0]), float(m[1])], "D_at_mode": Dm, "D_at_seed": Dseed,
                   "secant_root_over_km": [root.real / km, root.imag / km],
                   "root_minus_mode_over_km": abs(root - wm) / km})
        print("control", pc[-1], flush=True)
    out["positive_control"] = pc
    rng = random.Random(20261019)
    rows = []
    for i in range(int(sys.argv[2]) if len(sys.argv) > 2 else 120):
        u = rng.random()
        l = int(round(20 * (100 ** u)))            # log-uniform in [20, 2000]
        hinv = l + 0.5
        reg = rng.choice(["L", "S", "T", "A", "Abig"])
        lo, hi = {"L": (0, 8 * km), "S": (6 * km, 0.9 * OmP * hinv), "T": (0.8 * OmP * hinv, 2 * OmP * hinv),
                  "A": (1.5 * OmP * hinv, 3 * OmP * hinv), "Abig": (3 * OmP * hinv, 10 * OmP * hinv)}[reg]
        w = complex(lo + (hi - lo) * rng.random(), -d * rng.random())
        if rng.random() < 0.15:
            w = complex(w.real, -d)                 # the strip edge
        if rng.random() < 0.1:
            w = complex(w.real, 0.0)                # the real axis
        Dv, yp, yc, nf = D_stat(l, w)
        rows.append({"l": l, "regime": reg, "w_over_km": [w.real / km, w.imag / km], "D": Dv, "nfev": nf})
        print(i, l, reg, round(w.real / km, 4), round(w.imag / km, 4), "D=%.3e" % Dv, flush=True)
    out["rows"] = rows
    out["min_D"] = min(r["D"] for r in rows)
    out["n"] = len(rows)
    Path(sys.argv[1] if len(sys.argv) > 1 else "indep_points.json").write_text(json.dumps(out, indent=1) + "\n")
    print("n", len(rows), "min D", out["min_D"])


if __name__ == "__main__":
    main()
