"""Independent check of the S2 contour geometry (reviewer script, 2026-10-17).

Does NOT use the author's ODE integrator.  Checks:
  P. parameters recomputed with mpmath (60 digits) from the exact definitions and compared with the
     S1 certified balls in certificates/params.json (which this script uses for all Arb work);
  X. x_P = r*(r_P), Re r*(0), sum A_i = 0, x_P < Re r*(0);
  C. the contour image r(s), s in [-64, 32]: an Arb Krawczyk enclosure of the root of the
     branch-tracked equation r*(r) = x_P + e^{i pi/4} s (branches from a float path), compared with
     the author's Taylor-model enclosure (m3_contour.Ray) at every step endpoint and every step
     midpoint;
  V. the vertical segments at Re x = Re r*(0) and Re x = 0 (float, scipy DOP853), plus a fan of
     vertical lines covering both closed sectors, and the location of the nearest r = 0 / r = inf
     point on the two special lines (how far the singularities are from the sector);
  T. the tail discs, recomputed in Arb with a triangle-inequality bound (not the author's box
     bound): sup |r - r_+| for Re x <= x_P - 64/sqrt2, sup |r - r_c| for Re x >= x_P + 32/sqrt2;
     float continuation of the contour to s = -400 and s = +200 stays inside the discs.
Run (repo root):  .venv/bin/python problems/P2/certificates/review_s2/indep_geometry.py
"""
import json
import math
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np
from flint import acb, arb, ctx
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
P2 = HERE.parents[1]
sys.path.insert(0, str(P2 / "src"))
ctx.prec = 128
mp.mp.dps = 60
T0 = time.time()
out = {}


# ---------------------------------------------------------------- P: parameters
def mp_params():
    cub = sorted([z.real for z in mp.polyroots([1, 0, -25, 25], maxsteps=400, extraprec=400)])
    r_ext = min(x for x in cub if x > 0)
    Qmax2 = r_ext - r_ext ** 4 / 50
    Q2 = (mp.mpf(999) / 1000) ** 2 * Qmax2
    co = [-mp.mpf(1) / 50, 0, 1, -2, Q2]                      # descending
    rts = sorted([z.real for z in mp.polyroots(co, maxsteps=400, extraprec=400)])
    dD = lambda r: -4 * r ** 3 / 50 + 2 * r - 2
    A = [r * r / dD(r) for r in rts]
    kap = [abs(dD(r)) / (2 * r * r) for r in rts]
    rP = (3 + mp.sqrt(9 - 8 * Q2)) / 2
    FP = (-rP ** 4 / 50 + rP ** 2 - 2 * rP + Q2) / rP ** 2
    OmP = mp.sqrt(FP) / rP
    lamP = mp.sqrt(FP * (2 - 3 / rP)) / rP
    return {"r_ext": r_ext, "Q2": Q2, "r_n": rts[0], "r_minus": rts[1], "r_plus": rts[2], "r_c": rts[3],
            "A": A, "kappa_minus": kap[1], "kappa_plus": kap[2], "kappa_c": kap[3], "r_P": rP,
            "Omega_P": OmP, "gamma0": lamP / 2, "beta": (mp.mpf(1) / 2 + mp.mpf(1) / 20) * kap[1]}


M = mp_params()
J = json.loads((P2 / "certificates" / "params.json").read_text())["points"]["PB"]["quantities"]
B = {k: arb(J[k]["ball"]) for k in ("r_ext", "Q2", "r_n", "r_minus", "r_plus", "r_c", "kappa_minus",
                                    "kappa_plus", "kappa_c", "r_P", "Omega_P", "gamma0")}
B["beta"] = arb(J["strip_depth_beta"]["ball"])
AB = [arb(J["A_coefficients"][k]["ball"]) for k in ("r_n", "r_minus", "r_plus", "r_c")]
ROOTS = [B["r_n"], B["r_minus"], B["r_plus"], B["r_c"]]
contain = {}
for k in B:
    xm = arb(mp.nstr(M[k], 55))
    contain[k] = bool(B[k].overlaps(xm)) and abs(float((B[k] - xm).mid())) < 1e-35
contain["A_i"] = all(bool(a.overlaps(arb(mp.nstr(m, 55)))) for a, m in zip(AB, M["A"]))
out["P_params_mp_vs_params_json"] = contain
print("P params (mpmath 60 digits inside params.json balls):", contain, flush=True)

# ---------------------------------------------------------------- X: x_P, Re r*(0)
def re_rstar_arb(r):
    return sum((a * abs(r - ri).log() for a, ri in zip(AB, ROOTS)), arb(0))


xP = re_rstar_arb(acb(B["r_P"]))
x0 = sum((a * abs(ri).log() for a, ri in zip(AB, ROOTS)), arb(0))
sumA = sum(AB, arb(0))
out["X"] = {"x_P": xP.str(15), "Re_rstar_0": x0.str(15), "sum_A": sumA.str(5),
            "margin_Re_rstar0_minus_xP": (x0 - xP).str(10), "xP_lt_Re_rstar0": bool(xP < x0),
            "author_x_P": "[-2.43796603456595 +/- 7.67e-16]", "author_Re_rstar_0": "[-2.29663689319 +/- 4.87e-12]"}
print("X", out["X"], flush=True)

# float model (numpy) for paths and branch tracking
rn, rm, rp, rc = (float(M[k]) for k in ("r_n", "r_minus", "r_plus", "r_c"))
RTS = np.array([rn, rm, rp, rc])
Af = np.array([float(a) for a in M["A"]])
xPf = float(xP.mid())
E1 = np.exp(1j * np.pi / 4)


def Ff(r):
    return -(r - rn) * (r - rm) * (r - rp) * (r - rc) / (50 * r * r)


def ray_float(s_end, n_out):
    """r(s) along x = x_P + e^{i pi/4} s from r_P, plus continuous args of (r - r_i) (c: r_c - r)."""
    f = lambda s, y: [ (E1 * Ff(y[0] + 1j * y[1])).real, (E1 * Ff(y[0] + 1j * y[1])).imag ]
    ss = np.linspace(0.0, s_end, n_out)
    sol = solve_ivp(f, (0.0, s_end), [float(M["r_P"]), 0.0], t_eval=ss, method="DOP853", rtol=1e-13, atol=1e-15)
    r = sol.y[0] + 1j * sol.y[1]
    z = np.stack([r - rn, r - rm, r - rp, rc - r])
    args = np.unwrap(np.angle(z), axis=1)
    return ss, r, args


# ---------------------------------------------------------------- C: Krawczyk enclosures of r(s)
EI4 = acb(1, 1) / arb(2).sqrt()


def krawczyk_r(s_ball, r_guess, alphas, rho0=1e-9):
    """Unique root of f(r) = sum A_i [Log(sigma_i (r - r_i) e^{-i a_i}) + i a_i] - x(s) in a box."""
    X = acb(xP) + EI4 * s_ball
    sig = [1, 1, 1, -1]
    rot = [acb(arb(0), arb(-a)).exp() for a in alphas]           # e^{-i a}
    m = acb(r_guess.real, r_guess.imag)

    def f(r):
        return sum((AB[i] * ((sig[i] * (r - ROOTS[i])) * rot[i]).log() + acb(0, 1) * AB[i] * arb(alphas[i])
                    for i in range(4)), acb(0)) - X

    def fp(r):
        return sum((AB[i] / (r - ROOTS[i]) for i in range(4)), acb(0))
    Y = 1 / acb(complex(fp(m).mid()))
    rho = rho0
    for _ in range(30):
        R = acb(m.real + arb(0, 1) * rho, m.imag + arb(0, 1) * rho)
        ok_branch = all(bool(((sig[i] * (R - ROOTS[i])) * rot[i]).real > 0) for i in range(4))
        if ok_branch:
            K = m - Y * f(m) + (1 - Y * fp(R)) * (R - m)
            if (K.real.lower() > R.real.lower() and K.real.upper() < R.real.upper()
                    and K.imag.lower() > R.imag.lower() and K.imag.upper() < R.imag.upper()):
                return K
        rho *= 4
    return None


from p2 import m3 as A3                     # author's code: ONLY to obtain the enclosures compared against
from p2.params import Point

pt = Point("PB", 256)
t_auth = time.time()
left = A3.Ray(pt, A3.theta_unit(), arb(-64), arb(0.0625))
right = A3.Ray(pt, A3.theta_unit(), arb(32), arb(0.0625))
t_auth = time.time() - t_auth
res_C = {}
for name, ray, S in (("left", left, -64.0), ("right", right, 32.0)):
    n = int(abs(S) * 32) + 1                                    # grid of spacing 1/32: step ends + midpoints
    ss, rr, args = ray_float(S, n)
    worst_sep, n_overlap, n_fail, max_rad_author, max_rad_mine, max_middiff = 0.0, 0, 0, 0.0, 0.0, 0.0
    for k in range(n):
        s = arb(float(ss[k]))                                    # exact dyadic (k/32)
        K = krawczyk_r(s, complex(rr[k]), [float(a) for a in args[:, k]])
        if K is None:
            n_fail += 1
            continue
        Ra = ray.at(s)
        ovl = Ra.real.overlaps(K.real) and Ra.imag.overlaps(K.imag)
        n_overlap += int(ovl)
        max_rad_author = max(max_rad_author, float(Ra.real.rad()), float(Ra.imag.rad()))
        max_rad_mine = max(max_rad_mine, float(K.real.rad()), float(K.imag.rad()))
        max_middiff = max(max_middiff, float(abs(acb(Ra.mid()) - acb(K.mid())).upper()))   # Arb, not doubles
    res_C[name] = {"points": n, "krawczyk_failed": n_fail, "author_ball_overlaps_krawczyk_ball": n_overlap,
                   "max_radius_author": max_rad_author, "max_radius_krawczyk": max_rad_mine,
                   "max_midpoint_difference": max_middiff,
                   "float_r_at_end": [rr[-1].real, rr[-1].imag], "winding_args_at_end": [float(a) for a in args[:, -1]]}
    print("C", name, res_C[name], flush=True)
out["C_contour"] = res_C
out["C_contour"]["author_Ray_build_seconds"] = round(t_auth, 2)


# ---------------------------------------------------------------- V: vertical segments and fan
def real_point(X):
    g = lambda r: sum(float(a) * mp.log(abs(r - ri)) for a, ri in zip(M["A"], (M["r_n"], M["r_minus"], M["r_plus"], M["r_c"]))) - X
    return float(mp.findroot(g, (M["r_plus"] + mp.mpf("1e-14"), M["r_c"] - mp.mpf("1e-14")), solver="anderson"))


def vertical(X, Y_end, n=2001, stop_small=1e-3, stop_big=1e3):
    r0 = real_point(X)
    sgn = 1 if Y_end > 0 else -1

    def f(Y, y):
        z = 1j * Ff(y[0] + 1j * y[1])
        return [z.real, z.imag]

    def ev_small(Y, y):
        return abs(y[0] + 1j * y[1]) - stop_small
    ev_small.terminal = True

    def ev_big(Y, y):
        return abs(y[0] + 1j * y[1]) - stop_big
    ev_big.terminal = True
    sol = solve_ivp(f, (0.0, Y_end), [r0, 0.0], method="DOP853", rtol=1e-12, atol=1e-14,
                    t_eval=np.linspace(0.0, Y_end, n), events=(ev_small, ev_big))
    r = sol.y[0] + 1j * sol.y[1]
    hit = None
    if sol.t_events[0].size:
        hit = ("r->0", float(sol.t_events[0][0]))
    elif sol.t_events[1].size:
        hit = ("r->inf", float(sol.t_events[1][0]))
    return r, hit, sol.t


x0f = float(x0.mid())
res_V = {}
for name, X, H in (("Re_x=Re_rstar(0)", x0f, x0f - xPf), ("Re_x=0", 0.0, -xPf)):
    r, hit, t = vertical(X, H)
    res_V[name] = {"height": H, "min_abs_r": float(np.min(np.abs(r))), "max_abs_r": float(np.max(np.abs(r))),
                   "singularity_within_segment": hit}
    for direction in (+1, -1):                                  # nearest singularity on the whole line
        r2, hit2, t2 = vertical(X, direction * 200.0, n=20001)
        res_V[name][f"nearest_singularity_{'up' if direction > 0 else 'down'}"] = hit2
print("V special segments:", res_V, flush=True)
fan = {"upper_right": [], "lower_left": []}
for X in np.concatenate([np.linspace(xPf + 1e-6, 0.5, 60), np.linspace(0.5, 40.0, 40)]):
    r, hit, _ = vertical(float(X), float(X - xPf), n=801)
    fan["upper_right"].append((float(X), float(np.min(np.abs(r))), float(np.max(np.abs(r))), hit))
for X in np.linspace(-80.0, xPf - 1e-6, 80):
    r, hit, _ = vertical(float(X), float(X - xPf), n=801)
    fan["lower_left"].append((float(X), float(np.min(np.abs(r))), float(np.max(np.abs(r))), hit))
for k, v in fan.items():
    hits = [e for e in v if e[3] is not None]
    res_V[f"fan_{k}"] = {"lines": len(v), "min_abs_r": min(e[1] for e in v), "max_abs_r": max(e[2] for e in v),
                         "lines_hitting_singularity": len(hits)}
print("V fan:", {k: res_V[k] for k in res_V if k.startswith("fan")}, flush=True)
out["V_vertical"] = res_V


# ---------------------------------------------------------------- T: tail discs (triangle inequality)
def tail_disc(h, d0, S):
    ih = 2 if h == "plus" else 3
    rh, Ah = ROOTS[ih], AB[ih]
    smin, smax = arb(0), arb(0)
    for i in range(4):
        if i == ih:
            continue
        dist = abs(rh - ROOTS[i])
        lo, hi = dist - d0, dist + d0
        assert lo > 0
        t1, t2 = AB[i] * lo.log(), AB[i] * hi.log()
        smin += arb(min(t1.lower(), t2.lower()))       # exact arb comparisons, no float rounding
        smax += arb(max(t1.upper(), t2.upper()))
    xS = xP + arb(S) * arb(2).sqrt() / 2
    d = ((xS - smin) / Ah).exp() if h == "plus" else ((xS - smax) / Ah).exp()
    return arb(d.upper())


dL = tail_disc("plus", arb("0.06"), -64)
dR = tail_disc("c", arb("0.6"), 32)
# start points inside d0 (Krawczyk enclosure at s = -64 and s = 32)
ssL, rrL, argsL = ray_float(-64.0, 2049)          # dense path: np.unwrap needs small arg increments
KL = krawczyk_r(arb(-64), complex(rrL[-1]), [float(a) for a in argsL[:, -1]])
ssR, rrR, argsR = ray_float(32.0, 1025)
KR = krawczyk_r(arb(32), complex(rrR[-1]), [float(a) for a in argsR[:, -1]])
startL = abs(KL - acb(ROOTS[2]))
startR = abs(KR - acb(ROOTS[3]))
# float continuation beyond the validated range
ssx, rrx, _ = ray_float(-400.0, 40001)
mL = np.abs(rrx - rp)[ssx <= -64.0]
ssy, rry, _ = ray_float(200.0, 20001)
mR = np.abs(rry - rc)[ssy >= 32.0]
out["T_tails"] = {"d_left_reviewer": dL.str(8), "d_left_author": "[0.00955711 +/- 1.94e-9]",
                  "d_right_reviewer": dR.str(8), "d_right_author": "[0.351053 +/- 4.63e-8]",
                  "start_left_|r(-64)-r_+|": startL.str(8), "start_left_inside_0.06": bool(startL < arb("0.06")),
                  "start_right_|r(32)-r_c|": startR.str(8), "start_right_inside_0.6": bool(startR < arb("0.6")),
                  "d_left_lt_0.06": bool(dL < arb("0.06")), "d_right_lt_0.6": bool(dR < arb("0.6")),
                  "float_sup_|r-r_+|_s_in[-400,-64]": float(mL.max()), "float_sup_|r-r_c|_s_in[32,200]": float(mR.max())}
print("T", out["T_tails"], flush=True)
out["seconds"] = round(time.time() - T0, 1)
(HERE / "indep_geometry.json").write_text(json.dumps(out, indent=1, default=str) + "\n")
print("done", out["seconds"], "s")
