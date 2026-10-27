"""S6 part (i) at depth e kappa_- (e = 23/25), uniformly in l >= l0 = 20 (h <= h0 = 2/41).

Pre-registration O4 and section 3 (c98ca16d; Amendments f82114bf).  S4's four regimes on the band,
with the S6 changes:
  T  S4's contour and splitting, design T6 (s1 = 12, phi_inf = -1, taper (3, 6)) on zeta in [0.8, 1.3] Omega_P;
  A  zeta_A = 6/5 Omega_P, and Lemma W7;
  S  S4's procedure (both sides, the closing box and the zeta-boxes of [0, 0.9] Omega_P) with the design grid
     extended (beta_common.S_GRID_EXT), nu in [-e kappa_-, 0];
  L  the band rows [0, 8] x [-19/20, -11/20] (kappa_- units), box side 1/20 (quadrisection to 1/80), with the
     pole-free series R~ = (1 + 2 rho) R in every box (W6'; beta_common.make_lchecker).
S, T and A cover nu in [-e kappa_-, 0] at once; L covers the band rows (the rows above are S4's L, inside P2-4).

Commands (cwd problems/P2/src or PYTHONPATH=src):
  python -m p2.beta_part1 T | A | cover | controls
  python -m p2.beta_part1 launch {S+|Sc|L|B7} [workers]     (detached, nice 10, resumable JSONL)
  python -m p2.beta_part1 run {S+|Sc|L|B7} worker nworkers
  python -m p2.beta_part1 assemble {S+|Sc|L|B7}
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

from flint import arb

from .arbutil import A, Timer, precision
from .beta_common import (A_DESIGN, CERT, EDGE, H0, L_BOX, L_IM, L_MIN_BOX, L_RE, PREC, RES, T_DESIGN, T_FALLBACK,
                          FloatTA, commit_tag, deep_point, make_lchecker, use_extended_s_grid, write_json)
from .m3 import fa

RES1 = RES / "part1"
CERT1 = CERT / "part1"
SRC = Path(__file__).resolve().parents[1]
ZETA_SU = Fraction(9, 10)
OMEGA_S_KM = Fraction(6)
OMEGA_L_KM = Fraction(8)


# ---------------------------------------------------------------------------------------------
# T and A (fast, foreground)
# ---------------------------------------------------------------------------------------------
def t_boxes():
    lo_, hi_ = (Fraction(x) for x in T_DESIGN["zeta_over_OmegaP"])
    w = Fraction(T_DESIGN["zeta_box_over_OmegaP"])
    n = int((hi_ - lo_) / w)
    return [(lo_ + k * w, lo_ + (k + 1) * w) for k in range(n)]


def run_T():
    from . import part1_T as P
    pt = deep_point()
    f = FloatTA(pt)
    rows = []
    with precision(PREC), Timer() as tm:
        geo, rays, tails = P.setup(pt)
        ch = P.TChecker(pt, rays, tails, T_DESIGN)
        for zb in t_boxes():
            r = P.certify_box(ch, pt, zb, T_DESIGN, h0=H0)
            design = T_DESIGN
            if not r["pass"]:
                # the registered per-box fallback: float maximin over T_FALLBACK, then Arb
                best = None
                for s1, phinf in T_FALLBACK:
                    d = dict(T_DESIGN, s1=s1, phi_inf=phinf)
                    m = f.box(d, zb, EDGE, H0)["F_min_over_astar"]
                    if best is None or m > best[0]:
                        best = (m, d)
                design = best[1]
                ch2 = P.TChecker(pt, rays, tails, design)
                r = P.certify_box(ch2, pt, zb, design, h0=H0)
                r["fallback_design"] = {"s1": design["s1"], "phi_inf": design["phi_inf"]}
            r["design_name"] = design["name"]
            rows.append(r)
            print("T", zb, r["pass"], r["mu_over_a_star"], r["cells"], flush=True)
    ok = bool(geo["A0_pass"] and all(r["pass"] for r in rows))
    return {"regime": "T", "point": "PB", "depth_over_kappa_minus": str(EDGE), "h0": str(H0), "precision_bits": PREC,
            "design": T_DESIGN, "fallback_grid": [{"s1": a, "phi_inf": b} for a, b in T_FALLBACK],
            "geometry_A0_pass": geo["A0_pass"], "boxes": rows, "pass": ok,
            "cells_total": sum(r["cells"] for r in rows), "timing_s": {"wall": round(tm.wall, 1), "cpu": round(tm.cpu, 1)},
            "commit": commit_tag()}


def run_A():
    from . import part1_T as P
    pt = deep_point()
    with precision(PREC), Timer() as tm:
        geo, rays, tails = P.setup(pt)
        ch = P.TChecker(pt, rays, tails, A_DESIGN)
        zA = Fraction(A_DESIGN["zeta_over_OmegaP"][0])
        r = P.certify_box(ch, pt, (zA, zA), A_DESIGN, h0=H0)
        mono = P.monotonicity_A(pt, A_DESIGN, H0)
    ok = bool(geo["A0_pass"] and r["pass"] and mono["pass"])
    return {"regime": "A", "point": "PB", "depth_over_kappa_minus": str(EDGE), "h0": str(H0), "precision_bits": PREC,
            "design": A_DESIGN, "geometry_A0_pass": geo["A0_pass"], "boxes": [r], "monotonicity_W7": mono, "pass": ok,
            "timing_s": {"wall": round(tm.wall, 1), "cpu": round(tm.cpu, 1)}, "commit": commit_tag()}


# ---------------------------------------------------------------------------------------------
# S (resumable)
# ---------------------------------------------------------------------------------------------
def s_jobs():
    from . import throat as T
    return [("closing", "0", str(T.ZETA_1))] + [("box", str(a), str(b)) for a, b in T.zboxes()]


def run_S(side_name: str, worker: int = 0, nworkers: int = 1):
    from . import throat as T
    from .part1 import load_bulk
    grid = use_extended_s_grid()
    side = 1 if side_name == "+" else -1
    tag = "plus" if side > 0 else "c"
    W_bound, _ = load_bulk()
    RES1.mkdir(parents=True, exist_ok=True)
    out = RES1 / f"S_{tag}.w{worker}.jsonl"
    done = set()
    if out.exists():
        done = {tuple(json.loads(x)["key"]) for x in out.read_text().splitlines() if x.strip()}
    pt = deep_point()
    g = T.FGeo(pt)
    with precision(T.PREC):
        for j, key in enumerate(s_jobs()):
            if j % nworkers != worker or key in done:
                continue
            with Timer() as tm:
                if key[0] == "closing":
                    leaves = [T.certify_closing(pt, g, side, W_bound)]
                else:
                    leaves = T.certify_zbox_refined(pt, g, side, Fraction(key[1]), Fraction(key[2]), W_bound)
            row = {"key": list(key), "pass": all(x["pass"] for x in leaves), "leaves": leaves, "commit": commit_tag(),
                   "grid_sizes": {k: len(v) for k, v in grid["used"].items()},
                   "timing_s": {"wall": round(tm.wall, 1), "cpu": round(tm.cpu, 1)},
                   "finished": time.strftime("%Y-%m-%dT%H:%M:%S")}
            with out.open("a") as fh:
                fh.write(json.dumps(row) + "\n")
            print("S", side_name, key, row["pass"], row["timing_s"], flush=True)


def assemble_S(side_name: str):
    from . import throat as T
    grid = use_extended_s_grid()
    side = 1 if side_name == "+" else -1
    tag = "plus" if side > 0 else "c"
    rows = {}
    for f in sorted(RES1.glob(f"S_{tag}.w*.jsonl")):
        for line in f.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                rows[tuple(r["key"])] = r
    jobs = s_jobs()
    missing = [j for j in jobs if j not in rows]
    leaves = [x for j in jobs if j in rows for x in rows[j]["leaves"]]
    iv = sorted((Fraction(x["zeta_over_OmegaP"][0]), Fraction(x["zeta_over_OmegaP"][1])) for x in leaves)
    tiles = bool(iv) and iv[0][0] == 0 and iv[-1][1] == ZETA_SU and all(a[1] == b[0] for a, b in itertools.pairwise(iv))
    pt = deep_point()
    with precision(PREC):
        x0re = sum((Ai * abs(ri).log() for Ai, ri in zip(pt.Acoef, pt.roots)), arb(0))
        tail_strip, eps_tiles = True, True
        tail_disc_max = 0.0
        for x in leaves:
            if x.get("closing"):
                continue
            Xa, k = fa(Fraction(x["X_a"])), fa(Fraction(x["k"]))
            ta = Fraction(x["design"][0]) + Fraction(x["design"][1])
            XL = Xa + side * fa(T.T_LEFT - ta) / k
            tail_strip = tail_strip and (bool(XL < x0re) if side > 0 else bool(XL > 0))
            em = Fraction(x["eps_max"])
            eps_tiles = eps_tiles and em == T.eps_max_of(pt, Fraction(x["zeta_over_OmegaP"][0]))
            ivs = sorted((Fraction(e["eps"][0]), Fraction(e["eps"][1])) for e in x["eps_boxes"])
            eps_tiles = eps_tiles and ivs[0][0] == 0 and ivs[-1][1] == em and \
                all(a[1] == b[0] for a, b in itertools.pairwise(ivs))
            for e in x["eps_boxes"]:
                if e.get("tail_disc"):
                    tail_disc_max = max(tail_disc_max, float(arb(e["tail_disc"]).upper()))
    ok = (not missing) and tiles and eps_tiles and all(x["pass"] for x in leaves) and tail_strip
    commits = sorted({rows[j]["commit"] for j in jobs if j in rows})
    cert = {"regime": "S", "side": tag, "point": "PB", "depth_over_kappa_minus": str(EDGE), "h0": str(H0),
            "precision_bits": T.PREC, "design_grid": grid,
            "constants": {"Omega_S_over_km": str(OMEGA_S_KM), "zeta_SU_over_OmegaP": str(ZETA_SU),
                          "zeta_1_over_OmegaP": str(T.ZETA_1), "zeta_box": str(T.ZBOX), "zeta_box_min": str(T.ZBOX_MIN),
                          "eps_levels": T.EPS_LEVELS, "t_L": str(T.T_LEFT), "cell_depth": T.MAX_DEPTH},
            "missing": [list(j) for j in missing], "zeta_leaves_tile_0_to_zeta_SU": tiles,
            "eps_boxes_tile_0_to_eps_max": eps_tiles, "tail_region_in_strip": tail_strip,
            "tail_disc_max": tail_disc_max, "n_zeta_leaves": len(leaves),
            "cpu_s_total": round(sum(rows[j]["timing_s"]["cpu"] for j in jobs if j in rows), 1),
            "commits": commits, "pass": bool(ok), "leaves": leaves}
    write_json(CERT1 / f"S_PB_{tag}.json", cert)
    return cert


# ---------------------------------------------------------------------------------------------
# L (resumable)
# ---------------------------------------------------------------------------------------------
def l_boxes():
    out = []
    n_re = int((L_RE[1] - L_RE[0]) / L_BOX)
    n_im = int((L_IM[1] - L_IM[0]) / L_BOX)
    for i in range(n_re):
        for j in range(n_im):
            out.append((L_RE[0] + i * L_BOX, L_RE[0] + (i + 1) * L_BOX, L_IM[1] - (j + 1) * L_BOX, L_IM[1] - j * L_BOX))
    return out


def run_L(worker: int = 0, nworkers: int = 1):
    from . import lowfreq as LF
    from .part1 import load_bulk
    W_bound, _ = load_bulk()
    RES1.mkdir(parents=True, exist_ok=True)
    out = RES1 / f"L.w{worker}.jsonl"
    done = set()
    if out.exists():
        done = {tuple(json.loads(x)["grid_box"]) for x in out.read_text().splitlines() if x.strip()}
    pt = deep_point()
    with precision(LF.PREC):
        ch = make_lchecker(pt, W_bound)
        for k, b in enumerate(l_boxes()):
            key = tuple(str(x) for x in b)
            if k % nworkers != worker or key in done:
                continue
            with Timer() as tm:
                leaves = LF.certify_box_refined(ch, b, L_MIN_BOX)
            row = {"grid_box": list(key), "pass": all(x["pass"] for x in leaves), "leaves": leaves, "commit": commit_tag(),
                   "timing_s": {"wall": round(tm.wall, 2), "cpu": round(tm.cpu, 2)}}
            with out.open("a") as fh:
                fh.write(json.dumps(row) + "\n")


def assemble_L():
    rows = {}
    for f in sorted(RES1.glob("L.w*.jsonl")):
        for line in f.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                rows[tuple(r["grid_box"])] = r
    grid = [tuple(str(x) for x in b) for b in l_boxes()]
    missing = [g for g in grid if g not in rows]
    tiles_ok = True
    leaves = []
    for g in grid:
        if g not in rows:
            continue
        G = tuple(Fraction(x) for x in g)
        lv = [tuple(Fraction(v) for v in x["box"]) for x in rows[g]["leaves"]]
        inside = all(G[0] <= b[0] and b[1] <= G[1] and G[2] <= b[2] and b[3] <= G[3] for b in lv)
        disjoint = all(min(a[1], b[1]) <= max(a[0], b[0]) or min(a[3], b[3]) <= max(a[2], b[2])
                       for i, a in enumerate(lv) for b in lv[i + 1:])
        area = sum((b[1] - b[0]) * (b[3] - b[2]) for b in lv) == (G[1] - G[0]) * (G[3] - G[2])
        tiles_ok = tiles_ok and inside and disjoint and area
        leaves += rows[g]["leaves"]
    worst_p = min((float(arb(x["plus"]["Re_y_bound"]).lower()) for x in leaves if "Re_y_bound" in x.get("plus", {})),
                  default=None)
    worst_c = max((float(arb(x["c"]["Re_y_bound"]).upper()) for x in leaves if "Re_y_bound" in x.get("c", {})),
                  default=None)
    ok = (not missing) and tiles_ok and all(x["pass"] for x in leaves)
    from . import lowfreq as LF
    cert = {"regime": "L", "point": "PB", "depth_over_kappa_minus": str(EDGE), "h0": str(H0), "precision_bits": LF.PREC,
            "design": {"Re_range_over_km": [str(x) for x in L_RE], "Im_range_over_km": [str(x) for x in L_IM],
                       "box_side_over_km": str(L_BOX), "min_box_side_over_km": str(L_MIN_BOX),
                       "h2_range": ["0", str(H0 * H0)], "h2_subintervals": LF.K_H, "tau_grid": str(LF.TAU_STEP),
                       "exit_rule": "D/r^4 >= max(4|w|^2_max, C_W)", "series": "R~ = (1 + 2 rho) R in every box (W6')"},
            "grid_boxes": len(grid), "missing": [list(m) for m in missing], "leaves": len(leaves),
            "refined_grid_boxes": sum(1 for g in grid if g in rows and len(rows[g]["leaves"]) > 1),
            "tiles_ok": tiles_ok, "min_Re_y_plus_lower": worst_p, "max_Re_y_c_upper": worst_c,
            "cpu_s_total": round(sum(r["timing_s"]["cpu"] for r in rows.values()), 1),
            "commits": sorted({r["commit"] for r in rows.values()}), "pass": bool(ok),
            "leaves_detail": leaves}
    write_json(CERT1 / "L_PB.json", cert)
    return cert


# ---------------------------------------------------------------------------------------------
# the cover (section 3) and the part-(i) controls NB3-NB5
# ---------------------------------------------------------------------------------------------
def cover():
    pt = deep_point()
    certs = {}
    for n in ("T_PB", "A_PB", "S_PB_plus", "S_PB_c", "L_PB"):
        f = CERT1 / f"{n}.json"
        certs[n] = json.loads(f.read_text()) if f.exists() else None
    zTL, zTU = (Fraction(x) for x in T_DESIGN["zeta_over_OmegaP"])
    zA = Fraction(A_DESIGN["zeta_over_OmegaP"][0])
    with precision(PREC):
        km, OmP, d = pt.km, pt.OmegaP, pt.beta
        ineq = {
            "Omega_S < Omega_L (L/S overlap)": OMEGA_S_KM < OMEGA_L_KM,
            "zeta_TL < zeta_SU (S/T overlap)": zTL < ZETA_SU,
            "zeta_A < zeta_TU (T/A overlap)": zA < zTU,
            "zeta_SU Omega_P/h0 > Omega_L (S reaches above L for every h <= h0)":
                bool(fa(ZETA_SU) * OmP / fa(H0) > fa(OMEGA_L_KM) * km),
            "Omega_S kappa_- > d (Re w + Im w > 0 in S)": bool(fa(OMEGA_S_KM) * km > d),
            "zeta_TL Omega_P/h0 > d (T)": bool(fa(zTL) * OmP / fa(H0) > d),
            "zeta_A Omega_P/h0 > d (A)": bool(fa(zA) * OmP / fa(H0) > d),
            "d < 2 kappa_+ and d < kappa_c (only omega0 = -i kappa_+ of E_+ u E_c in S_e)":
                bool(d < 2 * pt.kp and d < pt.kc),
            "kappa_+ < d (omega0 lies in S_e)": bool(pt.kp < d),
            "L's band rows reach -e kappa_- and overlap S_delta": bool(A(L_IM[0]) <= -A(EDGE) and L_IM[1] >= Fraction(-11, 20)),
            "h0 = 2/41 gives l0 = 20": H0 == Fraction(2, 41),
        }
    tiles = {}
    if certs["T_PB"]:
        iv = sorted((Fraction(b["zeta_over_OmegaP"][0]), Fraction(b["zeta_over_OmegaP"][1])) for b in certs["T_PB"]["boxes"])
        tiles["T zeta-boxes tile [zeta_TL, zeta_TU]"] = bool(iv and iv[0][0] == zTL and iv[-1][1] == zTU and
                                                              all(a[1] == b[0] for a, b in itertools.pairwise(iv)))
        tiles["T nu_low encloses -e kappa_- and h0 = 2/41"] = all(
            arb(b["nu_low"]).contains(-A(EDGE) * pt.km) and b["h0"] == str(H0) for b in certs["T_PB"]["boxes"])
    if certs["A_PB"]:
        b = certs["A_PB"]["boxes"][0]
        tiles["A at zeta_A, nu_low encloses -e kappa_-, W7"] = bool(
            Fraction(b["zeta_over_OmegaP"][0]) == zA and arb(b["nu_low"]).contains(-A(EDGE) * pt.km)
            and certs["A_PB"]["monotonicity_W7"]["pass"])
    for tag in ("plus", "c"):
        c = certs[f"S_PB_{tag}"]
        if c:
            tiles[f"S {tag}: zeta-leaves tile [0, zeta_SU], eps-boxes tile [0, eps_max]"] = bool(
                c["zeta_leaves_tile_0_to_zeta_SU"] and c["eps_boxes_tile_0_to_eps_max"])
    if certs["L_PB"]:
        tiles["L: in every grid box of the band rows, the leaves are disjoint and fill it"] = bool(certs["L_PB"]["tiles_ok"])
    passes = {k: (v is not None and bool(v.get("pass"))) for k, v in certs.items()}
    ok = all(ineq.values()) and all(passes.values()) and tiles and all(tiles.values()) and len(tiles) == 6
    return {"point": "PB", "depth_over_kappa_minus": str(EDGE), "h0": str(H0),
            "constants": {"Omega_L_over_km": str(OMEGA_L_KM), "Omega_S_over_km": str(OMEGA_S_KM),
                          "zeta_SU": str(ZETA_SU), "zeta_TL": str(zTL), "zeta_TU": str(zTU), "zeta_A": str(zA)},
            "inequalities": {k: bool(v) for k, v in ineq.items()}, "certificates_pass": passes, "tiling": tiles,
            "cited": {"Im w >= -0.55 kappa_- (S_delta)": "P2-4 (S5)",
                      "Re w < 0": "W_l(-conj w) = conj W_l(w) (FORMULATION 1.2; zeq_identity Z4)",
                      "omega0 = -i kappa_+": "Lemma E (certificates/beta/epole_PB.json; W6' in L)"},
            "pass": bool(ok), "commit": commit_tag()}


def controls():
    """NB3: L's box containing omega0 with S4's (not pole-free) series must not certify.  NB4: T with S4's
    design at depth e on [0.8, 0.85] Omega_P must not certify.  NB5: A at h0 = 2 with zeta_A = 5 (< 21/4)
    must not certify."""
    from . import lowfreq as LF
    from . import part1_T as P
    from .part1 import load_bulk
    pt = deep_point()
    out = {}
    W_bound, _ = load_bulk()
    with precision(PREC):
        ch = LF.LChecker(pt, W_bound)
        b = (Fraction(0), Fraction(1, 20), Fraction(-17, 20), Fraction(-4, 5))
        rec = ch.box(b)
        ok_y = all(rec.get(s, {}).get("y_sign_ok") for s in ("plus", "c"))
        out["NB3_L_pole_box_S4_series"] = {"box": [str(x) for x in b], "plus": rec.get("plus"),
                                           "certified_y_signs": bool(ok_y), "fails_as_required": not ok_y}
        _, rays, tails = P.setup(pt)
        chT = P.TChecker(pt, rays, tails, P.DESIGN_T)
        r = P.certify_box(chT, pt, (Fraction(4, 5), Fraction(17, 20)), P.DESIGN_T, h0=H0)
        out["NB4_T_S4_design_at_e"] = {"zeta": ["4/5", "17/20"], "mu_over_a_star": r["mu_over_a_star"],
                                       "pass": r["pass"], "fails_as_required": not r["pass"]}
        chA = P.TChecker(pt, rays, tails, P.DESIGN_A)
        r = P.certify_box(chA, pt, (Fraction(5), Fraction(5)), P.DESIGN_A, h0=Fraction(2))
        out["NB5_A_h0_2_zeta_5"] = {"zeta_A": "5", "h0": "2", "mu_over_a_star": r["mu_over_a_star"],
                                    "pass": r["pass"], "fails_as_required": not r["pass"]}
    out["pass"] = all(v["fails_as_required"] for v in out.values())
    out["commit"] = commit_tag()
    return out


# ---------------------------------------------------------------------------------------------
# B7: float cross-check of the band for l >= 20 (A6; not used by the proof)
# ---------------------------------------------------------------------------------------------
def b7_points(n=100, seed=20261027):
    import math
    import random
    rng = random.Random(seed)
    pts = []
    for i in range(n):
        l = int(round(math.exp(rng.uniform(math.log(20), math.log(1000)))))
        reg = rng.choice(["L", "S", "T", "A"])
        pts.append({"i": i, "l": l, "regime": reg, "u": rng.random(), "v": rng.random()})
    return pts


def run_B7(worker: int, nworkers: int):
    import mpmath as mp

    from .qnm_float import Geom, Solver
    RES1.mkdir(parents=True, exist_ok=True)
    out = RES1 / f"B7.w{worker}.jsonl"
    done = set()
    if out.exists():
        done = {json.loads(x)["i"] for x in out.read_text().splitlines() if x.strip()}
    g = Geom.frozen("PB")
    with mp.workdps(20):
        sol = Solver(g, K=6)
        km, OmP = g.km, g.OmegaP
        e = mp.mpf(EDGE.numerator) / EDGE.denominator
        for p in b7_points():
            if p["i"] % nworkers != worker or p["i"] in done:
                continue
            l = p["l"]
            hinv = l + mp.mpf(1) / 2
            lo, hi = {"L": (0, 8 * km), "S": (6 * km, mp.mpf("0.9") * OmP * hinv),
                      "T": (mp.mpf("0.8") * OmP * hinv, mp.mpf("1.3") * OmP * hinv),
                      "A": (mp.mpf("1.2") * OmP * hinv, 3 * OmP * hinv)}[p["regime"]]
            im = -(mp.mpf("0.55") + (e - mp.mpf("0.55")) * p["v"]) * km
            w = mp.mpc(lo + (hi - lo) * p["u"], im)
            t0 = time.time()
            W = sol.wronskian(w, l)
            row = {**p, "w_over_km": [mp.nstr(mp.re(w) / km, 10), mp.nstr(mp.im(w) / km, 10)], "abs_W": mp.nstr(abs(W), 6)}
            try:
                root = sol.find(l, w, maxsteps=12)
                inband = bool(mp.re(root) >= 0 and -e * km <= mp.im(root) <= -mp.mpf("0.55") * km)
                Wr = abs(sol.wronskian(root, l))
                row.update({"secant_root_over_km": [mp.nstr(mp.re(root) / km, 10), mp.nstr(mp.im(root) / km, 10)],
                            "root_in_band": inband, "residual_ratio": float(Wr / abs(W)),
                            "is_zero_by_A6_rule": bool(inband and Wr < mp.mpf("1e-8") * abs(W))})
            except (ValueError, ZeroDivisionError) as ex:
                row.update({"secant": str(ex)[:80], "root_in_band": False, "is_zero_by_A6_rule": False})
            row["seconds"] = round(time.time() - t0, 1)
            with out.open("a") as fh:
                fh.write(json.dumps(row) + "\n")


def assemble_B7():
    rows = {}
    for f in sorted(RES1.glob("B7.w*.jsonl")):
        for line in f.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                rows[r["i"]] = r
    pts = b7_points()
    missing = [p["i"] for p in pts if p["i"] not in rows]
    flagged = [r for r in rows.values() if r.get("is_zero_by_A6_rule")]
    # a flagged stop would also be evaluated with zsol at 128 bits (A6); none expected
    cert = {"check": "B7 (A6): float W_l and a secant search at 100 points of the band, l log-uniform in [20, 1000]; "
                     "float evidence, not used by the proof", "points": len(pts), "missing": missing,
            "secant_stops_in_band": sum(1 for r in rows.values() if r.get("root_in_band")),
            "zeros_by_A6_rule": len(flagged), "flagged": flagged,
            "min_abs_W": min((float(r["abs_W"]) for r in rows.values()), default=None),
            "cpu_s_total": round(sum(r["seconds"] for r in rows.values()), 1),
            "pass": bool(not missing and not flagged), "rows": [rows[i] for i in sorted(rows)], "commit": commit_tag()}
    write_json(CERT1 / "B7_PB.json", cert)
    return cert


# ---------------------------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------------------------
def launch(job: str, nworkers: int = 1):
    RES1.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%dT%H%M%S")
    rec = {"job": job, "workers": nworkers, "commit": commit_tag(), "started": stamp, "host": platform.node(),
           "python": sys.executable, "nice": 10, "pids": [], "logs": []}
    for w in range(nworkers):
        log = RES1 / f"{job}_w{w}_{stamp}.log"
        cmd = ["nice", "-n", "10", sys.executable, "-m", "p2.beta_part1", "run", job, str(w), str(nworkers)]
        with log.open("w") as fh:
            p = subprocess.Popen(cmd, cwd=SRC, stdout=fh, stderr=subprocess.STDOUT, start_new_session=True)
        rec["pids"].append(p.pid)
        rec["logs"].append(log.name)
    (RES1 / f"launch_{job}_{stamp}.json").write_text(json.dumps(rec, indent=1) + "\n")
    return rec


def main(argv):
    cmd = argv[1]
    if cmd in ("T", "A"):
        c = run_T() if cmd == "T" else run_A()
        write_json(CERT1 / f"{cmd}_PB.json", c)
        print(cmd, "pass:", c["pass"], c["timing_s"], c["commit"])
    elif cmd == "cover":
        c = cover()
        write_json(CERT1 / "cover_PB.json", c)
        print("cover pass:", c["pass"])
    elif cmd == "controls":
        c = controls()
        write_json(CERT1 / "controls_PB.json", c)
        print("controls pass:", c["pass"])
    elif cmd == "launch":
        print(json.dumps(launch(argv[2], int(argv[3]) if len(argv) > 3 else 1), indent=1))
    elif cmd == "run":
        job, w, n = argv[2], int(argv[3]), int(argv[4])
        if job in ("S+", "Sc"):
            run_S(job[1], w, n)
        elif job == "L":
            run_L(w, n)
        elif job == "B7":
            run_B7(w, n)
        else:
            raise SystemExit(f"unknown job {job}")
    elif cmd == "assemble":
        job = argv[2]
        c = assemble_S(job[1]) if job in ("S+", "Sc") else (assemble_L() if job == "L" else assemble_B7())
        print(job, "pass:", c.get("pass"))
    else:
        print(__doc__)
        print(f"(pid {os.getpid()})")


if __name__ == "__main__":
    main(sys.argv)
