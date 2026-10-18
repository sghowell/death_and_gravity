"""S4: P2-1 (part i) at P_B uniformly in l >= l0 = 20 -- driver, bulk facts, cover, controls, replay.

Pre-registration: notes/s4-preregistration.md (with Amendments).  Proof: notes/s4-proof.md.
Regimes: L (lowfreq.py), S (throat.py), T and A (part1_T.py).  Certificates: certificates/part1/.

    python -m p2.part1 bulk                       bulk_PB.json (U monotone, W~ bound, C_W)
    python -m p2.part1 launch T|A|S+|Sc|L [n]     detached workers (nice 10) + a launch record
    python -m p2.part1 run REGIME w n             one worker (resumable JSONL under results/part1/)
    python -m p2.part1 assemble REGIME            certificate from the worker files
    python -m p2.part1 cover                      cover_PB.json
    python -m p2.part1 controls                   controls_PB.json (N1-N5, K1-K4)
    python -m p2.part1 verify                     replay_PB.json
"""
from __future__ import annotations

import itertools
import json
import os
import platform
import subprocess
import sys
import time
from fractions import Fraction
from pathlib import Path

from flint import arb, fmpq

from .arbutil import Timer, git_commit, precision
from .m3 import fa
from .params import Point

PREC = 128
ROOT = Path(__file__).resolve().parents[2]
SRC = Path(__file__).resolve().parents[1]
CERT_DIR = ROOT / "certificates" / "part1"
RES_DIR = ROOT / "results" / "part1"

# regime boundaries (pre-registration section 2/5); kappa_- and Omega_P multiples, exact
H0 = Fraction(2, 41)
OMEGA_L_KM = Fraction(8)
OMEGA_S_KM = Fraction(6)
ZETA_SU, ZETA_TL, ZETA_TU, ZETA_A = Fraction(9, 10), Fraction(4, 5), Fraction(2), Fraction(3, 2)


# ---------------------------------------------------------------------------------------------
# bulk facts
# ---------------------------------------------------------------------------------------------
def Wt_real(pt: Point, r):
    """W~ = F F_r/r - U/4 for real r (arb ball or interval)."""
    D = pt.Delta(r)
    Dp = pt.dDelta(r)
    F = D / (r * r)
    Fr = (Dp * r - 2 * D) / r ** 3
    return F * Fr / r - F / (4 * r * r)


def U_real(pt: Point, r):
    return pt.Delta(r) / r ** 4


def bulk(pt: Point, n_cells=6000):
    """U monotone on (r_+, r_P) and (r_P, r_c); a certified bound max(-W~) <= W_bound on [r_+, r_c];
    C_W = 4 max(kappa_+^2, W_bound)."""
    import sympy as sp
    r, Lam, Q2 = sp.symbols("r Lambda Q2")
    D = -Lam / 3 * r ** 4 + r ** 2 - 2 * r + Q2
    ident = sp.expand(sp.diff(D, r) * r - 4 * D + 2 * (r ** 2 - 3 * r + 2 * Q2)) == 0
    # V_l = U/h^2 + W~ exactly (h = 1/(l + 1/2), L = l(l+1) = 1/h^2 - 1/4), with V_l = F(L/r^2 + F'/r)
    l = sp.symbols("l")
    Fs = D / r ** 2
    V = Fs * (l * (l + 1) / r ** 2 + sp.diff(Fs, r) / r)
    Ws = Fs * sp.diff(Fs, r) / r - Fs / r ** 2 / 4
    split = sp.simplify(V - ((l + sp.Rational(1, 2)) ** 2 * Fs / r ** 2 + Ws)) == 0
    with Timer() as tm:
        rP_inner = 3 - pt.rP                     # the other root of r^2 - 3r + 2Q^2 (sum of roots 3)
        mono = bool(rP_inner < pt.rp) and bool(pt.rP < pt.rc)
        a, b = pt.rp, pt.rc
        lo_a, hi_b = a.lower(), b.upper()
        lo_min = None
        argmin = None
        for k in range(n_cells):
            x0 = lo_a + (hi_b - lo_a) * k / n_cells
            x1 = lo_a + (hi_b - lo_a) * (k + 1) / n_cells
            cell = (x0 + x1) / 2 + arb(0, 1) * ((x1 - x0) / 2)
            w = Wt_real(pt, cell)
            if lo_min is None or w.lower() < lo_min:
                lo_min, argmin = w.lower(), float(((x0 + x1) / 2).mid())
        maxneg = -arb(lo_min)
        W_bound = Fraction(int(float(maxneg.upper()) * 1e5) + 1, 10 ** 5)
        assert fa(W_bound) > maxneg
        CW = 4 * arb(fa(W_bound)).max(pt.kp * pt.kp)
    return {"identity_Delta_prime_r_minus_4Delta": "Delta' r - 4 Delta = -2 (r^2 - 3 r + 2 Q^2)",
            "identity_sympy_ok": bool(ident), "splitting_V_eq_U_over_h2_plus_Wt_sympy_ok": bool(split),
            "inner_photon_sphere_3_minus_rP": rP_inner.str(15), "r_plus": pt.rp.str(15),
            "U_increasing_on_rplus_rP_and_decreasing_on_rP_rc": mono,
            "Wt_cells": n_cells, "Wt_min_lower": arb(lo_min).str(10), "Wt_min_near_r": argmin,
            "W_bound": str(W_bound), "W_bound_float": float(W_bound),
            "kappa_plus_sq": (pt.kp * pt.kp).str(10), "C_W": CW.str(12),
            "C_W_upper_exact": str(Fraction(4) * max(W_bound, Fraction(int(float((pt.kp * pt.kp).upper()) * 1e9) + 1, 10 ** 9))),
            "pass": bool(ident and split and mono),
            "timing_s": {"wall": round(tm.wall, 1), "cpu": round(tm.cpu, 1)}}


def load_bulk():
    st = json.loads((CERT_DIR / "bulk_PB.json").read_text())
    return Fraction(st["W_bound"]), st


# ---------------------------------------------------------------------------------------------
# negative controls N1-N5 (pre-registration section 4, Amendment 2)
# ---------------------------------------------------------------------------------------------
def _ps_mode(l: int):
    import mpmath as mp

    from .qnm_float import Geom, Solver
    g = Geom.frozen("PB")
    with mp.workdps(30):
        top = (l + mp.mpf(1) / 2) * g.OmegaP
        return Solver(g, K=6).find(l, mp.mpc(top, -g.gamma0))


def control_N1_N2(pt: Point):
    """N1: the T checker on boxes holding the l-th PS fundamental mode (l = 20, 50, 100), at h0 = 1/(l+1/2);
    IN (+-gamma0/50) and STRAD (Im in [-1.05, -0.95] gamma0).  N2: h0 = 1/20.5, zeta in
    [Omega_P - 0.001 kappa_-, Omega_P], nu in [-1.69, -1.68] kappa_-.  Each must NOT certify."""
    import mpmath as mp

    from . import part1_T as TT
    _, rays, tails = TT.setup(pt)
    ch = TT.TChecker(pt, rays, tails, TT.DESIGN_T)
    out = {"N1": {}, "N2": None}
    g0 = pt.gamma0
    e = g0 / 50
    for l in (20, 50, 100):
        ps = _ps_mode(l)
        re_ps, im_ps = arb(mp.nstr(mp.re(ps), 40)), arb(mp.nstr(mp.im(ps), 40))
        h = Fraction(2, 2 * l + 1)
        hh = fa(h)
        row = {"PS_mode_over_kappa_minus": [mp.nstr(mp.re(ps) / float(pt.km.mid()), 12), mp.nstr(mp.im(ps) / float(pt.km.mid()), 12)]}
        for name, (nlo, nhi) in (("IN", (im_ps - e, im_ps + e)), ("STRAD", (-arb("1.05") * g0, -arb("0.95") * g0))):
            r = TT.certify_box_abs(ch, pt, hh * (re_ps - e), hh * (re_ps + e), nlo, TT.DESIGN_T, h0=h)
            r["contains_mode"] = bool(nlo < im_ps < nhi)
            row[name] = r
        row["pass"] = (not row["IN"]["certified"]) and (not row["STRAD"]["certified"]) and row["IN"]["contains_mode"] \
            and row["STRAD"]["contains_mode"]
        out["N1"][str(l)] = row
    km = pt.km
    r = TT.certify_box_abs(ch, pt, pt.OmegaP - arb("0.001") * km, pt.OmegaP, -arb("1.69") * km, TT.DESIGN_T, h0=H0)
    r["pass"] = not r["certified"]
    out["N2"] = r
    return out


def _n3_box(pt: Point, LF, l: int, re_lo, re_hi, im_lo, im_hi, CW):
    """The L checker's continuation path at integer l on the box (kappa_- units): exit rule, then
    Re y_+(x_a) > 0 and Re y_c(x_a^c) < 0."""
    from flint import acb
    km = pt.km
    w = acb(LF._hull(fa(re_lo) * km, fa(re_hi) * km), LF._hull(fa(im_lo) * km, fa(im_hi) * km))
    h2 = Fraction(1) / (Fraction(l) + Fraction(1, 2)) ** 2
    M = arb((4 * arb(abs(w).upper()) ** 2).max(CW).upper())
    row = {"box_re_over_km": [str(re_lo), str(re_hi)], "box_im_over_km": [str(im_lo), str(im_hi)]}
    im_ = LF._hull(fa(im_lo) * km, fa(im_hi) * km)
    row["box_avoids_E_plus_E_c"] = all(not (im_ + n * kap).contains(0) and bool(abs(im_ + n * kap) > 0)
                                       for kap in (pt.kp, pt.kc) for n in range(1, 8))
    ok = True
    for side, name in ((1, "plus"), (-1, "c")):
        tau = LF.exit_tau_point(pt, side, M, h2)
        exit_rule = bool(LF.scaled_U(pt, side, fmpq(tau.numerator, tau.denominator), fa(h2)) >= M)
        radius = (pt.rp - pt.rm) / fa(h2) if side > 0 else (pt.rc - pt.rp) / fa(h2)
        tau0 = min(tau, Fraction(int(float(radius.lower()) / 4 * 1024), 1024))
        y = None
        while y is None:
            try:
                y, nst = LF.y_continued_tau(pt, side, w, h2, tau0, tau)
            except ArithmeticError:
                tau0 = tau0 / 2
                if tau0 < Fraction(1, 2 ** 20):
                    raise
        good = bool(y.re_lower() > 0) if side > 0 else bool(y.re_upper() < 0)
        row[name] = {"tau_a": str(tau), "exit_rule_ok": exit_rule, "tau0_frobenius": str(tau0), "taylor_steps_in_tau": nst,
                     "Re_y_lo_hi": [float(y.re_lower()), float(y.re_upper())], "sign_certified": good}
        ok = ok and good and exit_rule
    row["certified"] = ok
    return row


def control_N3(pt: Point, W_bound: Fraction):
    """N3: the L checker through its continuation path at l = 1, 2 on boxes holding the NE modes
    -1.911925 i kappa_- and -2.891324 i kappa_- (Re in [0, 1/20000], Im = mode +- 1/40000, kappa_- units:
    the continuation toward the exit runs against the recessive direction, so a box's deviation is
    amplified about 30-fold and the boxes are kept small).  Must NOT certify.  Reported alongside (not
    required): the same checker on a mode-free box at half the mode's depth, to show that the failure at
    the mode is a sign failure, not wide enclosures."""
    from . import lowfreq as LF
    CW = 4 * arb(fa(W_bound)).max(pt.kp * pt.kp)
    out = {}
    for l, mode in ((1, Fraction(-1911925, 10 ** 6)), (2, Fraction(-2891324, 10 ** 6))):
        d = Fraction(1, 40000)
        row = {"l": l, "mode_over_km": str(mode)}
        row["IN"] = _n3_box(pt, LF, l, Fraction(0), 2 * d, mode - d, mode + d, CW)
        row["OUT_half_depth"] = _n3_box(pt, LF, l, Fraction(0), 2 * d, mode / 2 - d, mode / 2 + d, CW)
        row["pass"] = (not row["IN"]["certified"]) and row["IN"]["box_avoids_E_plus_E_c"]
        out[str(l)] = row
    return out


def control_N4(pt: Point, W_bound: Fraction):
    """N4: the S checker with the exit placed before the turning point (allowed region): the design selected
    for zeta in [0.3, 0.32] Omega_P with X_a moved to x_T - 1/(2k) (side +) / x_T + 1/(2k) (side c).
    Must NOT certify."""
    from . import throat as TH
    g = TH.FGeo(pt)
    out = {}
    zlo, zhi = Fraction(3, 10), Fraction(8, 25)
    zl_f, zh_f = float(zlo) * g.OmP, float(zhi) * g.OmP
    for side, name in ((1, "plus"), (-1, "c")):
        D, _ = TH.select_design(g, side, zl_f, zh_f, float(TH.eps_max_of(pt, zlo)), float(W_bound))
        import math
        zref = math.sqrt(zl_f * zh_f)
        rT = g.turning(zref, side)
        k = abs(g.Ux(rT)) / (2 * zref ** 2)
        Xa = TH.dyadic(g.rstar(rT) - side * 0.5 / k)
        # the cells near the misplaced exit have q^ < 0 on a whole interval, so deeper refinement cannot change
        # the verdict: this control runs with eps_levels = 0 and cell depth 2 (production: 4 and 8)
        rec = TH.certify_zbox(pt, g, side, zlo, zhi, W_bound,
                              design={"design": D.key(), "X_a": str(Xa), "k": str(TH.dyadic(k))},
                              eps_levels=0, max_depth=2)
        out[name] = {"design": D.key(), "X_a": str(Xa), "refinement": {"eps_levels": 0, "cell_depth": 2}, "exit": rec.get("exit"),
                     "eps_boxes": [{"mu_cells": e["mu_cells"], "n_fails": e["n_fails"]} for e in rec.get("eps_boxes", [])],
                     "error": rec.get("error"), "certified": rec["pass"], "pass": not rec["pass"]}
    return out


def control_N5(pt: Point, W_bound: Fraction):
    """N5 (Amendment 2): the L checker on its production path (h^2 in [0, h0^2], no continuation) with the
    exit in the allowed region, D/r^4 >= (Re w)_min^2/4 instead of the exit rule.  Must NOT certify."""
    from . import lowfreq as LF
    ch = LF.LChecker(pt, W_bound)
    out = {}
    for b in ((Fraction(159, 20), Fraction(8), Fraction(-11, 20), Fraction(-1, 2)),
              (Fraction(159, 20), Fraction(8), Fraction(-3, 10), Fraction(-1, 4)),
              (Fraction(159, 20), Fraction(8), Fraction(-1, 20), Fraction(0)),
              (Fraction(6), Fraction(121, 20), Fraction(-1, 20), Fraction(0))):
        rec = ch.box(b, exit_scale=Fraction(1, 4))
        rec["certified"] = rec["pass"]
        rec["pass"] = not rec["certified"]
        out["/".join(str(x) for x in b)] = rec
    return out


# ---------------------------------------------------------------------------------------------
# acceptance checks K1-K4
# ---------------------------------------------------------------------------------------------
K2_POINTS = [(0, 0), (2, -0.2), (5, -0.55), (8, -0.3)]


def check_K2(ls=(20, 30, 50), points=K2_POINTS):
    """K2: W_l from L's scaled-Frobenius Jost values (continued by Taylor steps) against zsol and arbw
    (S3's C0 wrapper) at exact dyadic points; both comparison balls need relative radius < 1e-10."""
    from flint import acb

    from . import part3, zsol
    from .lowfreq import wronskian_L
    pt = Point("PB", 256)
    out = []
    for l in ls:
        for re_k, im_k in points:
            w_q = (part3._kq(pt, re_k), part3._kq(pt, im_k))
            with precision(256):
                w = acb(arb(w_q[0]), arb(w_q[1]))
                WL = wronskian_L(pt, l, w)
                Wz = zsol.wronskian_z(pt, l, w, theta=part3.THETA_POINT)[0]
                rl = float(abs(WL).rad() / abs(WL).mid())
                rz = float(abs(Wz).rad() / abs(Wz).mid())
            Wa, pa, ra = part3._arbw_reference("PB", l, w_q)
            with precision(256):
                ov = bool(WL.overlaps(Wz) and WL.overlaps(Wa) and Wz.overlaps(Wa))
            out.append({"l": l, "w_kappa": [re_k, im_k], "w_exact": [str(w_q[0]), str(w_q[1])],
                        "W_L": WL.str(15), "relrad_L": rl, "relrad_zsol": rz, "relrad_arbw": ra, "arbw_prec": pa,
                        "overlap_all": ov, "pass": bool(ov and rl < 1e-10 and rz < 1e-10 and ra < 1e-10)})
    return {"rows": out, "pass": all(r["pass"] for r in out)}


def _k4_points(n=200, seed=20261018):
    """Deterministic sample of the cover: l in [20, 400], a regime, Re(omega) in its range, nu in [-d, 0]."""
    import random
    rng = random.Random(seed)
    pts = []
    for i in range(n):
        l = rng.randint(20, 400)
        reg = rng.choice(["L", "S", "T", "A"])
        u = rng.random()
        v = rng.random()
        pts.append({"i": i, "l": l, "regime": reg, "u": u, "v": v})
    return pts


def run_K4(worker: int, nworkers: int):
    """K4 worker: float W_l at the sample points and a secant QNM search seeded there (float evidence)."""
    import mpmath as mp

    from .qnm_float import Geom, Solver
    RES_DIR.mkdir(parents=True, exist_ok=True)
    out = RES_DIR / f"K4_PB.w{worker}.jsonl"
    done = set()
    if out.exists():
        done = {json.loads(x)["i"] for x in out.read_text().splitlines() if x.strip()}
    g = Geom.frozen("PB")
    with mp.workdps(20):
        sol = Solver(g, K=6)
        km, OmP = g.km, g.OmegaP
        d = mp.mpf("0.55") * km
        for p in _k4_points():
            if p["i"] % nworkers != worker or p["i"] in done:
                continue
            l = p["l"]
            hinv = l + mp.mpf(1) / 2
            lo, hi = {"L": (0, 8 * km), "S": (6 * km, mp.mpf("0.9") * OmP * hinv),
                      "T": (mp.mpf("0.8") * OmP * hinv, 2 * OmP * hinv),
                      "A": (mp.mpf("1.5") * OmP * hinv, 3 * OmP * hinv)}[p["regime"]]
            w = mp.mpc(lo + (hi - lo) * p["u"], -d * p["v"])
            t0 = time.time()
            W = sol.wronskian(w, l)
            row = {**p, "w_over_km": [mp.nstr(mp.re(w) / km, 10), mp.nstr(mp.im(w) / km, 10)],
                   "abs_W": mp.nstr(abs(W), 6), "W_nonzero": bool(abs(W) > 0)}
            try:
                root = sol.find(l, w, maxsteps=12)
                inH = bool(mp.re(root) >= 0 and -d <= mp.im(root) <= 0)
                row["secant_root_over_km"] = [mp.nstr(mp.re(root) / km, 10), mp.nstr(mp.im(root) / km, 10)]
                row["root_in_H_plus"] = inH
            except (ValueError, ZeroDivisionError) as ex:
                row["secant"] = str(ex)[:80]
                row["root_in_H_plus"] = False
            row["seconds"] = round(time.time() - t0, 1)
            with out.open("a") as f:
                f.write(json.dumps(row) + "\n")


def check_K3():
    """K3: at 3 sample (zeta, nu, h) points in each overlap, the adjacent certificates both hold there
    (looked up in the assembled certificates)."""
    pt = Point("PB", 256)
    km, OmP = float(pt.km.mid()), float(pt.OmegaP.mid())
    certs = {n: json.loads((CERT_DIR / f"{n}.json").read_text()) for n in
             ("T_PB", "A_PB", "S_PB_plus", "S_PB_c", "L_PB_plus", "L_PB_c")}

    def in_T(z):
        return any(float(Fraction(b["zeta_over_OmegaP"][0])) <= z / OmP <= float(Fraction(b["zeta_over_OmegaP"][1]))
                   and b["pass"] for b in certs["T_PB"]["boxes"])

    def in_S(side, z, eps):
        for x in certs[f"S_PB_{side}"]["leaves"]:
            a, b = (float(Fraction(v)) for v in x["zeta_over_OmegaP"])
            if a <= z / OmP <= b and x["pass"]:
                if x.get("closing"):
                    return eps <= float(Fraction(x["eps_max"]))
                for e in x["eps_boxes"]:
                    if float(Fraction(e["eps"][0])) <= eps <= float(Fraction(e["eps"][1])) and e["pass"]:
                        return True
        return False

    def in_L(side, re_k, im_k):
        for x in certs[f"L_PB_{side}"]["leaves_detail"]:
            b = [float(Fraction(v)) for v in x["box"]]
            if b[0] <= re_k <= b[1] and b[2] <= im_k <= b[3] and x["y_sign_ok"] and x["exit_rule_ok"]:
                return True
        return False

    rows = []
    for (h, re_k, im_k) in ((1 / 20.5, 7.0, -0.3), (1 / 100.5, 6.5, -0.1), (1 / 1000.5, 7.9, -0.5)):
        re = re_k * km
        rows.append({"overlap": "L/S", "h": h, "re_over_km": re_k, "im_over_km": im_k,
                     "L": in_L("plus", re_k, im_k) and in_L("c", re_k, im_k),
                     "S": in_S("plus", h * re, 1 / re) and in_S("c", h * re, 1 / re)})
    for (h, zf, im_k) in ((1 / 20.5, 0.85, -0.3), (1 / 50.5, 0.82, -0.55), (1 / 200.5, 0.89, 0.0)):
        z = zf * OmP
        rows.append({"overlap": "S/T", "h": h, "zeta_over_OmegaP": zf, "im_over_km": im_k,
                     "S": in_S("plus", z, h / z) and in_S("c", z, h / z), "T": in_T(z)})
    for (h, zf, im_k) in ((1 / 20.5, 1.6, -0.2), (1 / 40.5, 1.9, -0.55), (1 / 400.5, 1.55, 0.0)):
        z = zf * OmP
        rows.append({"overlap": "T/A", "h": h, "zeta_over_OmegaP": zf, "im_over_km": im_k, "T": in_T(z),
                     "A": bool(certs["A_PB"]["pass"] and zf >= 1.5 and certs["A_PB"]["monotonicity_W7"]["pass"])})
    for r in rows:
        r["pass"] = all(v for k, v in r.items() if k in ("L", "S", "T", "A"))
    return {"rows": rows, "pass": all(r["pass"] for r in rows)}


def cover():
    """cover_PB.json: the exact inequalities of the regime cover, the tiling of each certificate's
    parameter domain, and every regime certificate's verdict."""
    pt = Point("PB", 256)
    out = {"constants": {"h0": str(H0), "Omega_L_over_km": str(OMEGA_L_KM), "Omega_S_over_km": str(OMEGA_S_KM),
                         "zeta_SU": str(ZETA_SU), "zeta_TL": str(ZETA_TL), "zeta_TU": str(ZETA_TU), "zeta_A": str(ZETA_A)}}
    with precision(PREC):
        km, OmP, d = pt.km, pt.OmegaP, pt.beta
        ineq = {
            "Omega_S < Omega_L (L/S overlap)": OMEGA_S_KM < OMEGA_L_KM,
            "zeta_TL < zeta_SU (S/T overlap)": ZETA_TL < ZETA_SU,
            "zeta_A < zeta_TU (T/A overlap)": ZETA_A < ZETA_TU,
            "zeta_SU Omega_P/h0 > Omega_L (S reaches above L for every h <= h0)":
                bool(fa(ZETA_SU) * OmP / fa(H0) > fa(OMEGA_L_KM) * km),
            "d_delta < kappa_+ < kappa_c (E_+, E_c outside the strip)": bool(d < pt.kp and pt.kp < pt.kc),
            "Omega_S > d_delta (Jost decay on S's tails: Re w + Im w > 0)": bool(fa(OMEGA_S_KM) * km > d),
            "zeta_TL Omega_P/h0 > d_delta (T)": bool(fa(ZETA_TL) * OmP / fa(H0) > d),
            "zeta_A Omega_P/h0 > d_delta (A)": bool(fa(ZETA_A) * OmP / fa(H0) > d),
        }
    out["inequalities"] = {k: bool(v) for k, v in ineq.items()}
    certs = {}
    for n in ("bulk_PB", "T_PB", "A_PB", "S_PB_plus", "S_PB_c", "L_PB_plus", "L_PB_c"):
        f = CERT_DIR / f"{n}.json"
        certs[n] = json.loads(f.read_text()) if f.exists() else None
    out["certificates_pass"] = {k: (v is not None and bool(v.get("pass"))) for k, v in certs.items()}
    tiles = {}
    T = certs["T_PB"]
    if T:
        iv = sorted((Fraction(b["zeta_over_OmegaP"][0]), Fraction(b["zeta_over_OmegaP"][1])) for b in T["boxes"])
        tiles["T zeta-boxes tile [zeta_TL, zeta_TU]"] = bool(iv and iv[0][0] == ZETA_TL and iv[-1][1] == ZETA_TU
                                                              and all(a[1] == b[0] for a, b in itertools.pairwise(iv)))
        tiles["T nu range [-d_delta, 0], h0"] = all(b["h0"] == str(H0) for b in T["boxes"])
    A = certs["A_PB"]
    if A:
        tiles["A at zeta_A, monotone above (W7)"] = bool(A["boxes"][0]["zeta_over_OmegaP"][0] == str(ZETA_A)
                                                        and A["monotonicity_W7"]["pass"])
    for side in ("plus", "c"):
        Sx = certs[f"S_PB_{side}"]
        if Sx:
            ok = Sx["zeta_leaves_tile_0_to_zeta_SU"]
            from .throat import eps_max_of
            for x in Sx["leaves"]:
                zl = Fraction(x["zeta_over_OmegaP"][0])
                em = Fraction(x["eps_max"])
                ok = ok and em == eps_max_of(pt, zl)
                if not x.get("closing"):
                    ivs = sorted((Fraction(e["eps"][0]), Fraction(e["eps"][1])) for e in x["eps_boxes"])
                    ok = ok and ivs[0][0] == 0 and ivs[-1][1] == em and all(a[1] == b[0] for a, b in itertools.pairwise(ivs))
            tiles[f"S {side}: zeta-leaves tile [0, zeta_SU], eps-boxes tile [0, eps_max(zeta_lo)]"] = bool(ok)
        Lx = certs[f"L_PB_{side}"]
        if Lx:
            area = sum((Fraction(b["box"][1]) - Fraction(b["box"][0])) * (Fraction(b["box"][3]) - Fraction(b["box"][2]))
                       for b in Lx["leaves_detail"])
            inside = all(0 <= Fraction(b["box"][0]) and Fraction(b["box"][1]) <= OMEGA_L_KM
                         and Fraction(-11, 20) <= Fraction(b["box"][2]) and Fraction(b["box"][3]) <= 0 for b in Lx["leaves_detail"])
            tiles[f"L {side}: leaves inside [0, 8] x [-11/20, 0] with total area 8 * 11/20 (quadrisected boxes)"] = \
                bool(inside and area == OMEGA_L_KM * Fraction(11, 20))
    out["tiling"] = tiles
    out["cited"] = {"Im w > 0": "P2-0 (at P2-0's level)", "Re w < 0": "W_l(-conj w) = conj W_l(w) (FORMULATION 1.2; zeq_identity Z4)"}
    out["pass"] = bool(all(out["inequalities"].values()) and all(out["certificates_pass"].values()) and tiles
                       and all(tiles.values()))
    out["commit"] = git_commit()
    return out


def k4_residual(row):
    """A secant 'root' counts only if the Wronskian is small there: |W(root)| < 1e-8 |W(seed)| (float, the
    run's settings).  Each flagged root is also evaluated with the rigorous zsol code (256 bits)."""
    import mpmath as mp
    from flint import acb

    from . import part3, zsol
    from .qnm_float import Geom, Solver
    g = Geom.frozen("PB")
    with mp.workdps(20):
        sol = Solver(g, K=6)
        km = g.km
        root = mp.mpc(mp.mpf(row["secant_root_over_km"][0]), mp.mpf(row["secant_root_over_km"][1])) * km
        seed = mp.mpc(mp.mpf(row["w_over_km"][0]), mp.mpf(row["w_over_km"][1])) * km
        Wr, Ws = abs(sol.wronskian(root, row["l"])), abs(sol.wronskian(seed, row["l"]))
        ratio = float(Wr / Ws)
    pt = Point("PB", 256)
    with precision(256):
        w = acb(arb(mp.nstr(mp.re(root), 30)), arb(mp.nstr(mp.im(root), 30)))
        Wz = zsol.wronskian_z(pt, row["l"], w, theta=part3.THETA_POINT)[0]
    return {"i": row["i"], "l": row["l"], "root_over_km": row["secant_root_over_km"], "abs_W_root_float": mp.nstr(Wr, 8),
            "abs_W_seed_float": mp.nstr(Ws, 8), "residual_ratio": ratio, "zsol_abs_W_root": abs(Wz).str(10),
            "zsol_excludes_0": bool(not Wz.contains(0)), "is_zero": bool(ratio < 1e-8 and Wz.contains(0))}


def assemble_K4():
    rows = {}
    for f in sorted(RES_DIR.glob("K4_PB.w*.jsonl")):
        for line in f.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                rows[r["i"]] = r
    n = len(_k4_points())
    missing = [i for i in range(n) if i not in rows]
    vals = [rows[i] for i in sorted(rows)]
    flagged = [r for r in vals if r["root_in_H_plus"]]
    checks = [k4_residual(r) for r in flagged]
    return {"points": n, "missing": missing, "all_W_nonzero": all(r["W_nonzero"] for r in vals),
            "secant_converged": sum(1 for r in vals if "secant_root_over_km" in r),
            "secant_stopped_in_H_plus_raw": len(flagged),
            "residual_checks": checks,
            "zeros_in_H_plus": [c for c in checks if c["is_zero"]],
            "by_regime": {g: sum(1 for r in vals if r["regime"] == g) for g in "LSTA"},
            "seconds_total": round(sum(r["seconds"] for r in vals), 1),
            "pass": (not missing) and all(r["W_nonzero"] for r in vals) and not any(c["is_zero"] for c in checks)
                    and all(c["zsol_excludes_0"] for c in checks),
            "rows": vals}


def controls():
    """controls_PB.json: N1-N5 (must not certify) and K1-K4."""
    from . import part1_T as TT
    pt = Point("PB", 256)
    W_bound, _ = load_bulk()
    out = {"commit": git_commit()}
    with precision(PREC):
        with Timer() as tm:
            out["N1_N2"] = control_N1_N2(pt)
            out["N3"] = control_N3(pt, W_bound)
            out["N4"] = control_N4(pt, W_bound)
            out["N5"] = control_N5(pt, W_bound)
            geo, _, _ = TT.setup(pt)
            out["K1_geometry_replay"] = {"A0_pass": geo["A0_pass"], "pass": bool(geo["A0_pass"])}
        out["timing_N_K1_s"] = {"wall": round(tm.wall, 1), "cpu": round(tm.cpu, 1)}
    with Timer() as tm:
        out["K2"] = check_K2()
    out["K2"]["timing_s"] = {"wall": round(tm.wall, 1), "cpu": round(tm.cpu, 1)}
    out["K3"] = check_K3()
    out["K4"] = assemble_K4()
    n1 = all(r["pass"] for r in out["N1_N2"]["N1"].values()) and out["N1_N2"]["N2"]["pass"]
    n3 = all(r["pass"] for r in out["N3"].values())
    n4 = all(r["pass"] for r in out["N4"].values())
    n5 = all(r["pass"] for r in out["N5"].values())
    out["summary"] = {"N1_N2": n1, "N3": n3, "N4": n4, "N5": n5, "K1": out["K1_geometry_replay"]["pass"],
                      "K2": out["K2"]["pass"], "K3": out["K3"]["pass"], "K4": out["K4"]["pass"]}
    out["kill_K_unsound"] = not (n1 and n3)
    out["pass"] = all(out["summary"].values())
    return out


def verify(frac_stride=20):
    """replay_PB.json: K1; bulk; T/A (5% of cells from the stored partitions, incl. each box's worst cell,
    and the tails exactly); S (every 20th zeta-leaf recomputed from its stored design and target); L
    (every 20th omega-leaf); the cover file."""
    from . import lowfreq as LF
    from . import part1_T as TT
    from . import throat as TH
    pt = Point("PB", 256)
    probs = []
    stats = {}
    with precision(PREC):
        geo, rays, tails = TT.setup(pt)
        if not geo["A0_pass"]:
            probs.append("K1: geometry replay failed")
        b = bulk(pt)
        st = json.loads((CERT_DIR / "bulk_PB.json").read_text())
        if b["W_bound"] != st["W_bound"] or b["pass"] != st["pass"]:
            probs.append("bulk: replay differs")
        for name, design in (("T", TT.DESIGN_T), ("A", TT.DESIGN_A)):
            cert = json.loads((CERT_DIR / f"{name}_PB.json").read_text())
            ch = TT.TChecker(pt, rays, tails, design)
            nc = 0
            for bx in cert["boxes"]:
                zb = tuple(Fraction(x) for x in bx["zeta_over_OmegaP"])
                cells = ch.cells_from_leaves(bx["leaves"])
                idx = set(range(0, len(cells), frac_stride))
                wc = bx["worst_cell"]
                if wc is not None:
                    idx |= {i for i, (a_, b_, _) in enumerate(cells) if float(a_) == wc[0] and float(b_) == wc[1]}
                r = TT.certify_box(ch, pt, zb, design, leaves=bx["leaves"], sample=idx)
                nc += r["cells"]
                mu_st = arb(bx["mu_cert"])
                if not (arb(r["mu_cert"]).lower() >= mu_st.lower()) or r["tail_left"] != bx["tail_left"] \
                        or r["tail_right"] != bx["tail_right"] or r["n_fails"] != 0:
                    probs.append(f"{name} box {bx['zeta_over_OmegaP']}: sample replay inconsistent")
            stats[name] = {"boxes": len(cert["boxes"]), "cells_recomputed": nc}
        g = TH.FGeo(pt)
        W_bound, _ = load_bulk()
        for side, tag in ((1, "plus"), (-1, "c")):
            cert = json.loads((CERT_DIR / f"S_PB_{tag}.json").read_text())
            leaves = cert["leaves"]
            pick = list(range(0, len(leaves), frac_stride))
            for i in pick:
                x = leaves[i]
                dsg = {"design": x["design"], "target": x["target"]}
                if x.get("closing"):
                    r = TH.certify_closing(pt, g, side, W_bound, design=dsg)
                    same = r["pass"] == x["pass"] and r["mu_cells"] == x["mu_cells"]
                else:
                    dsg.update({"X_a": x["X_a"], "k": x["k"]})
                    zl, zh = (Fraction(v) for v in x["zeta_over_OmegaP"])
                    r = TH.certify_zbox(pt, g, side, zl, zh, W_bound, design=dsg)
                    same = r["pass"] == x["pass"] and [e["mu_cells"] for e in r["eps_boxes"]] == \
                        [e["mu_cells"] for e in x["eps_boxes"]]
                if not same:
                    probs.append(f"S {tag} leaf {x['zeta_over_OmegaP']}: replay differs")
            stats[f"S_{tag}"] = {"leaves": len(leaves), "recomputed": len(pick)}
        ch = LF.LChecker(pt, W_bound)
        cert = json.loads((CERT_DIR / "L_PB_plus.json").read_text())
        certc = json.loads((CERT_DIR / "L_PB_c.json").read_text())
        lv, lvc = cert["leaves_detail"], certc["leaves_detail"]
        pick = list(range(0, len(lv), frac_stride))
        for i in pick:
            bx = tuple(Fraction(v) for v in lv[i]["box"])
            r = ch.box(bx)
            if r["plus"].get("Re_y_bound") != lv[i]["Re_y_bound"] or r["c"].get("Re_y_bound") != lvc[i]["Re_y_bound"]:
                probs.append(f"L leaf {lv[i]['box']}: replay differs")
        stats["L"] = {"leaves": len(lv), "recomputed": len(pick)}
    cv = cover()
    st = json.loads((CERT_DIR / "cover_PB.json").read_text())
    if cv["pass"] != st["pass"] or cv["tiling"] != st["tiling"] or cv["inequalities"] != st["inequalities"]:
        probs.append("cover: replay differs")
    return {"problems": probs, "pass": not probs, "stats": stats, "commit": git_commit()}


# ---------------------------------------------------------------------------------------------
# detached, resumable runs
# ---------------------------------------------------------------------------------------------
def launch(regime: str, nworkers: int = 1):
    RES_DIR.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%dT%H%M%S")
    rec = {"regime": regime, "workers": nworkers, "commit": git_commit(), "started": stamp,
           "host": platform.node(), "python": sys.executable, "nice": 10, "pids": [], "logs": []}
    for w in range(nworkers):
        log = RES_DIR / f"{regime}_w{w}_{stamp}.log"
        cmd = ["nice", "-n", "10", sys.executable, "-m", "p2.part1", "run", regime, str(w), str(nworkers)]
        with log.open("w") as fh:
            p = subprocess.Popen(cmd, cwd=SRC, stdout=fh, stderr=subprocess.STDOUT, start_new_session=True)
        rec["pids"].append(p.pid)
        rec["logs"].append(log.name)
        rec["cmd"] = " ".join(cmd[3:])
    (RES_DIR / f"launch_{regime}_{stamp}.json").write_text(json.dumps(rec, indent=1) + "\n")
    return rec


def run(regime: str, worker: int, nworkers: int):
    if regime in ("T", "A"):
        from . import part1_T
        part1_T.run_regime(regime, worker, nworkers)
    elif regime in ("S+", "Sc"):
        from . import throat
        throat.run_side(regime[1], worker, nworkers)
    elif regime == "L":
        from . import lowfreq
        lowfreq.run_L(worker, nworkers)
    elif regime == "K4":
        run_K4(worker, nworkers)
    else:
        raise SystemExit(f"unknown regime {regime}")


def assemble(regime: str):
    if regime in ("T", "A"):
        from . import part1_T
        return part1_T.assemble(regime)
    if regime in ("S+", "Sc"):
        from . import throat
        return throat.assemble(regime[1])
    if regime == "L":
        from . import lowfreq
        return lowfreq.assemble()
    raise SystemExit(f"unknown regime {regime}")


def main(argv):
    cmd = argv[1] if len(argv) > 1 else "help"
    if cmd == "bulk":
        CERT_DIR.mkdir(parents=True, exist_ok=True)
        pt = Point("PB", 256)
        with precision(PREC):
            c = bulk(pt)
        c["commit"] = git_commit()
        c["precision_bits"] = PREC
        (CERT_DIR / "bulk_PB.json").write_text(json.dumps(c, indent=1) + "\n")
        print(json.dumps(c, indent=1))
    elif cmd == "launch":
        print(json.dumps(launch(argv[2], int(argv[3]) if len(argv) > 3 else 1), indent=1))
    elif cmd == "run":
        run(argv[2], int(argv[3]), int(argv[4]))
    elif cmd == "assemble":
        c = assemble(argv[2])
        print(argv[2], "pass:", c.get("pass"))
    elif cmd == "cover":
        c = cover()
        (CERT_DIR / "cover_PB.json").write_text(json.dumps(c, indent=1) + "\n")
        print("cover pass:", c["pass"])
    elif cmd == "controls":
        c = controls()
        (CERT_DIR / "controls_PB.json").write_text(json.dumps(c, indent=1) + "\n")
        print("controls pass:", c["pass"])
    elif cmd == "verify":
        c = verify()
        (CERT_DIR / "replay_PB.json").write_text(json.dumps(c, indent=1) + "\n")
        print("replay:", "OK" if c["pass"] else c["problems"])
    else:
        print(__doc__)
        print(f"(pid {os.getpid()})")


if __name__ == "__main__":
    main(sys.argv)
