#!/usr/bin/env python
"""S5 review: an independent re-check of the P2-4 assembly bookkeeping at P_B.

Shares no code with the p2 package.  Exact rationals (fractions.Fraction) and rigorous interval
arithmetic (mpmath.iv) only; the parameters are recomputed from their definitions (FORMULATION 1.3)
by exact bisection, not read from params.json (params.json is only compared against).

What it checks (each item is a named boolean in the output):
  A  the point P_B from first principles: r_ext, Q_max^2, Q^2, the four roots of Delta, the surface
     gravities, r_P, Omega_P, gamma0; overlap with params.json; the photon-sphere uniqueness that
     HV17 2.2 uses ((r - r_P)(r^-2 F)' < 0 on (r_+, r_c));
  B  the seam inequalities: b - d_delta > 0 (and b is the least 2^-28 dyadic >= d_delta), 0 < e <= eta,
     b < kappa_+, d_delta < kappa_+, kappa_c, gamma0; kappa_- < 2 kappa_+ (and that it is implied by
     d_delta < kappa_+); exponents 21/20 and 11/10;
  C  for l = 0..19, from the raw JSON: a_l identity (lg vs winding), lg metadata, criterion upper < 2,
     strip depth, the right-half path (start, end, contiguity, edges, monotone), every piece inside
     its disc of radius rho/2, pole clearance of every disc box, exact re-sum of the increments with
     an interval for pi, winding values;
  D  l >= 20: h0 = 2/41 <=> l >= 20; the regime Re-ranges union to [0, oo) for every h <= h0;
  E  a sampled cover of S_delta x {l}: seam points and random points, each assigned to an exclusion;
  F  beta: the Krawczyk enclosure of omega_* (P2-K6) against 11/20 and 0.90938098, off E_+ u E_c;
     D4's W_0'(0) = i (r_+^2 + r_c^2)/(r_+ r_c).

Run from anywhere:  python indep_assembly.py   (writes indep_assembly_out.json next to this file)
"""
from __future__ import annotations

import json
import random
import re
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

import mpmath
from mpmath import iv

HERE = Path(__file__).resolve().parent
P2 = HERE.parents[1]                      # problems/P2
CERT = P2 / "certificates"
OUT = HERE / "indep_assembly_out.json"
iv.dps = 80
mpmath.mp.dps = 50          # display only; every decision uses iv or exact rationals

checks: dict[str, dict] = {}


def rec(name, ok, **info):
    checks[name] = {"ok": bool(ok), **{k: (str(v) if not isinstance(v, (bool, int, float, str, list, dict)) else v)
                                        for k, v in info.items()}}
    return bool(ok)


# ------------------------------------------------------------------------------------------------
# exact / interval helpers (independent of p2)

def ivq(q: Fr):
    """An mpmath interval containing the rational q (outward-rounded division)."""
    return iv.mpf(q.numerator) / iv.mpf(q.denominator)


def ivhull(lo: Fr, hi: Fr):
    a, b = ivq(lo), ivq(hi)
    return iv.mpf([a.a, b.b])


def _tuple_to_fr(t) -> Fr:
    """Exact value of an mpmath mpf tuple (sign, man, exp, bc); no rounding (mpmath.mpf(...) would
    round to the mp context precision, so it is not used for endpoints)."""
    sign, man, exp, _bc = t
    if man == 0:
        if exp != 0:
            raise ValueError("special value (inf/nan) in an interval endpoint")
        return Fr(0)
    v = Fr(int(man)) * (Fr(2) ** exp) if exp >= 0 else Fr(int(man), 2 ** (-exp))
    return -v if sign else v


def iv_lo(x) -> Fr:
    return _tuple_to_fr(x._mpi_[0])


def iv_hi(x) -> Fr:
    return _tuple_to_fr(x._mpi_[1])


BALL = re.compile(r"^\[?\s*([+-]?[0-9.eE+-]*)\s*(?:\+/-\s*([0-9.eE+-]+))?\s*\]?$")


def parse_ball(s: str) -> tuple[Fr, Fr]:
    """'[m +/- r]', '[+/- r]' or 'm' -> (m - r, m + r) as exact rationals (the printed Arb ball contains
    the true ball, so this interval is a valid enclosure)."""
    s = s.strip()
    if s.startswith("[") and "+/-" in s:
        body = s[1:-1]
        m_str, r_str = body.split("+/-")
        m = Fr(m_str.strip()) if m_str.strip() else Fr(0)
        r = Fr(r_str.strip())
        return m - r, m + r
    v = Fr(s.strip("[]"))
    return v, v


def ball_iv(s: str):
    lo, hi = parse_ball(s)
    return ivhull(lo, hi)


# ------------------------------------------------------------------------------------------------
# A: the point from first principles

def bisect_root(f, lo: Fr, hi: Fr, bits: int = 190):
    """f(x) returns (lo, hi) bounds of the value at a rational x; requires a certified sign change."""
    flo, fhi = f(lo), f(hi)
    s_lo = 1 if flo[0] > 0 else (-1 if flo[1] < 0 else 0)
    s_hi = 1 if fhi[0] > 0 else (-1 if fhi[1] < 0 else 0)
    assert s_lo * s_hi == -1, "no certified sign change"
    w = Fr(1, 2 ** bits)
    while hi - lo > w:
        mid = (lo + hi) / 2
        fm = f(mid)
        s = 1 if fm[0] > 0 else (-1 if fm[1] < 0 else 0)
        if s == 0:
            raise RuntimeError("sign undetermined; increase precision of the coefficients")
        if s == s_lo:
            lo = mid
        else:
            hi = mid
    return lo, hi


def point_PB():
    t0 = time.time()
    # r_ext: smallest positive root of r^3 - 25 r + 25 (exact polynomial)
    p = lambda r: (r ** 3 - 25 * r + 25,) * 2
    re_lo, re_hi = bisect_root(p, Fr(1), Fr(105, 100), bits=260)
    # Q_max^2 = r - r^4/50, increasing on the bracket (1 - 4 r^3/50 > 0 there)
    g = lambda r: r - r ** 4 / 50
    assert 1 - 4 * re_hi ** 3 / 50 > 0
    qm2 = (g(re_lo), g(re_hi))
    q2 = (Fr(999, 1000) ** 2 * qm2[0], Fr(999, 1000) ** 2 * qm2[1])
    D0 = lambda r: -r ** 4 / 50 + r * r - 2 * r
    Delta = lambda r: (D0(r) + q2[0], D0(r) + q2[1])
    brackets = {"r_n": (Fr(-8), Fr(-79, 10)), "r_minus": (Fr(99, 100), Fr(1)),
                "r_plus": (Fr(109, 100), Fr(11, 10)), "r_c": (Fr(58, 10), Fr(59, 10))}
    roots = {k: bisect_root(Delta, *v) for k, v in brackets.items()}
    # the four brackets are disjoint and Delta has degree 4: exactly one simple root in each
    order = [roots["r_n"], roots["r_minus"], roots["r_plus"], roots["r_c"]]
    disjoint = all(order[i][1] < order[i + 1][0] for i in range(3))
    R = {k: ivhull(*v) for k, v in roots.items()}
    Q2 = ivhull(*q2)

    def dDelta(r):
        return -iv.mpf(4) / 50 * r ** 3 + 2 * r - 2

    kap = {k: abs(dDelta(R[k])) / (2 * R[k] ** 2) for k in R}
    rP = (3 + iv.sqrt(9 - 8 * Q2)) / 2
    rPm = (3 - iv.sqrt(9 - 8 * Q2)) / 2
    F = lambda r: 1 - 2 / r + Q2 / r ** 2 - r ** 2 / 50
    FrP = F(rP)
    gamma0 = iv.sqrt(FrP * (2 - 3 / rP)) / (2 * rP)
    OmegaP = iv.sqrt(FrP) / rP
    return {"r_ext": ivhull(re_lo, re_hi), "Qmax2": ivhull(*qm2), "Q2": Q2, "roots": R, "kappa": kap,
            "r_P": rP, "r_P_minus": rPm, "F_rP": FrP, "gamma0": gamma0, "Omega_P": OmegaP,
            "roots_disjoint": disjoint, "time_s": round(time.time() - t0, 2)}


def mid(x, n=12) -> str:
    """Midpoint of an interval, for display only."""
    return mpmath.nstr((mpmath.mpf(x.a) + mpmath.mpf(x.b)) / 2, n)


def overlaps(x, y) -> bool:
    return not (x.b < y.a or y.b < x.a)


def main():
    t_start = time.time()
    out = {"tool": "indep_assembly.py (S5 review; shares no code with p2)", "mpmath": mpmath.__version__,
           "iv_dps": iv.dps}
    pt = point_PB()
    km, kp, kc = pt["kappa"]["r_minus"], pt["kappa"]["r_plus"], pt["kappa"]["r_c"]
    rp, rc, rm = pt["roots"]["r_plus"], pt["roots"]["r_c"], pt["roots"]["r_minus"]
    out["point"] = {k: mid(v, 20) for k, v in
                    [("kappa_minus", km), ("kappa_plus", kp), ("kappa_c", kc), ("gamma0", pt["gamma0"]),
                     ("Omega_P", pt["Omega_P"]), ("r_minus", rm), ("r_plus", rp), ("r_c", rc),
                     ("r_n", pt["roots"]["r_n"]), ("r_P", pt["r_P"])]}
    out["point"]["widths"] = {k: mpmath.nstr(mpmath.mpf(v.delta), 3) for k, v in [("kappa_minus", km), ("gamma0", pt["gamma0"])]}
    rec("A.roots_disjoint_four_simple", pt["roots_disjoint"])
    params = json.loads((CERT / "params.json").read_text())["points"]["PB"]
    Q = params["quantities"]
    for key, mine in [("kappa_minus", km), ("kappa_plus", kp), ("kappa_c", kc), ("gamma0", pt["gamma0"]),
                      ("Omega_P", pt["Omega_P"]), ("r_minus", rm), ("r_plus", rp), ("r_c", rc),
                      ("r_n", pt["roots"]["r_n"]), ("r_P", pt["r_P"]), ("Q2", pt["Q2"]), ("Qmax2", pt["Qmax2"]),
                      ("r_ext", pt["r_ext"])]:
        rec(f"A.params_overlap.{key}", overlaps(mine, ball_iv(Q[key]["ball"])))
    # HV17 2.2: (r^-2 F)' = -2 r^-5 (r^2 - 3 r + 2 Q^2) has exactly one zero in (r_+, r_c)
    rec("A.photon_sphere_unique_in_exterior", (pt["r_P_minus"].b < rp.a) and (rp.b < pt["r_P"].a)
        and (pt["r_P"].b < rc.a), r_P_minus=mid(pt["r_P_minus"], 12))
    rec("A.F_rP_pos_and_trapping", (pt["F_rP"].a > 0) and ((2 - 3 / pt["r_P"]).a > 0))

    # ---------------- B: seams
    dfrac = Fr(11, 20)
    d = ivq(dfrac) * km
    b = Fr(6280879, 2 ** 28)
    e = Fr(4567911, 2 ** 30)
    marg = ivq(b) - d
    rec("B.b_minus_d_delta_pos", marg.a > 0, margin=mid(marg, 6),
        margin_lower=mpmath.nstr(mpmath.mpf(marg.a), 6), width_of_d_delta=mpmath.nstr(mpmath.mpf(d.delta), 3))
    k_floor = Fr(6280878, 2 ** 28)
    rec("B.b_is_least_2^-28_dyadic_above_d", (ivq(k_floor) - d).b < 0)
    rec("B.e_pos_le_eta", e > 0 and (ivq(e) - km / 10).b <= 0, e_over_km=mid((ivq(e) / km), 10))
    rec("B.b_lt_kappa_plus", (ivq(b) - kp).b < 0, margin=mid((kp - ivq(b)), 6))
    rec("B.d_lt_kappa_plus", (d - kp).b < 0, margin=mid((kp - d), 6))
    rec("B.d_lt_kappa_c", (d - kc).b < 0, margin=mid((kc - d), 6))
    rec("B.d_lt_gamma0", (d - pt["gamma0"]).b < 0, margin=mid((pt["gamma0"] - d), 6),
        d_over_gamma0=mid((d / pt["gamma0"]), 6))
    rec("B.km_lt_2kp", (km - 2 * kp).b < 0, margin=mid((2 * kp - km), 6))
    # (1/2 + delta) kappa_- < kappa_+  ==>  2 kappa_+ > (1 + 2 delta) kappa_- > kappa_- : CF17's 2 kappa_+ > kappa_-
    # carries nothing beyond d_delta < kappa_+ (already needed for E_+ to miss S_delta)
    rec("B.km_lt_2kp_is_implied_by_d_lt_kp", (2 * dfrac) > 1)
    rec("B.HV17_exponent_21_20", Fr(1, 2) + dfrac == Fr(21, 20) and Fr(21, 20) > 1)
    rec("B.CF17_2min_over_km_11_10", 2 * dfrac == Fr(11, 10))

    # ---------------- C: l <= 19
    pi_iv = iv.pi
    per_l = {}
    b_set, e_set = set(), set()
    all_present = True
    for l in range(20):
        fl = CERT / "lg" / f"lg_PB_l{l}.json"
        fw = CERT / "part3" / f"winding_PB_l{l}.json"
        if not (fl.exists() and fw.exists()):
            all_present = False
            rec(f"C.l{l}.present", False)
            continue
        lg, wd = json.loads(fl.read_text()), json.loads(fw.read_text())
        reg = wd["region"]
        a_w, b_w, e_w = Fr(reg["a"]), Fr(reg["b"]), Fr(reg["e"])
        a_lg = Fr(lg["W_l_exact"])
        b_set.add(b_w); e_set.add(e_w)
        rec(f"C.l{l}.a_identical", a_w == a_lg)
        rec(f"C.l{l}.b_e_frozen", b_w == b and e_w == e)
        rec(f"C.l{l}.lg_meta", lg["point"] == "PB" and lg["l"] == l and Fr(lg["qfrac"]) == Fr(999, 1000)
            and Fr(lg["delta"]) == Fr(1, 20))
        crit_hi = parse_ball(lg["criterion"]["value_upper"])[1]
        rec(f"C.l{l}.lg_criterion_lt_2", lg["criterion"]["certified"] is True and crit_hi < 2,
            criterion_upper=float(crit_hi))
        rec(f"C.l{l}.lg_strip_depth_is_d", overlaps(ball_iv(lg["strip_depth_beta"]), d))
        rec(f"C.l{l}.wind_meta", wd["point"] == "PB" and wd["l"] == l and wd["certified"] is True
            and reg["a_source"] == f"certificates/lg/lg_PB_l{l}.json")
        pcs = wd["pieces"]
        # path: (0,-b) -> (a,-b) -> (a,e) -> (0,e), counterclockwise, contiguous, on the right edges
        P = lambda pnt: (Fr(pnt[0]), Fr(pnt[1]))
        ok = P(pcs[0]["p0"]) == (0, -b_w) and P(pcs[-1]["p1"]) == (0, e_w)
        corner1 = corner2 = False
        for i, s in enumerate(pcs):
            (x0, y0), (x1, y1) = P(s["p0"]), P(s["p1"])
            if i and P(pcs[i - 1]["p1"]) != (x0, y0):
                ok = False
            on_bottom = y0 == y1 == -b_w and 0 <= x0 < x1 <= a_w
            on_right = x0 == x1 == a_w and -b_w <= y0 < y1 <= e_w
            on_top = y0 == y1 == e_w and a_w >= x0 > x1 >= 0
            if not (on_bottom or on_right or on_top):
                ok = False
            corner1 |= (x1, y1) == (a_w, -b_w)
            corner2 |= (x1, y1) == (a_w, e_w)
        # a contiguous path of axis-parallel pieces from (0,-b) to (0,e) that stays on the three
        # edges and passes both corners traces the right half exactly once (monotone on each edge)
        rec(f"C.l{l}.path_right_half", ok and corner1 and corner2, pieces=len(pcs))
        # discs: both endpoints within rho/2 of the centre (the segment then lies in the disc)
        in_disc = all(all((Fr(p[0]) - Fr(s["c"][0])) ** 2 + (Fr(p[1]) - Fr(s["c"][1])) ** 2
                          <= (Fr(s["rho"]) / 2) ** 2 for p in (s["p0"], s["p1"])) for s in pcs)
        rec(f"C.l{l}.pieces_in_discs", in_disc)
        # pole clearance of the box [c +- rho]^2 (corner distance sqrt2 rho) from -i n kappa, n >= 1,
        # kappa in {kappa_+, kappa_c}: Im c + kappa > 0 (checked), so Im c + n kappa > 0 increases
        # with n and |c + i n kappa| >= |c + i kappa|; the exact n = 1 distance therefore covers all n.
        clear = True
        worst = None
        for s in pcs:
            cx, cy, rho = ivq(Fr(s["c"][0])), ivq(Fr(s["c"][1])), ivq(Fr(s["rho"]))
            for kap in (kp, kc):
                if not (cy + kap).a > 0:
                    clear = False
                dist2 = cx ** 2 + (cy + kap) ** 2
                ratio = dist2 / (2 * rho ** 2)
                if not ratio.a > 1:
                    clear = False
                if worst is None or ratio.a < worst:
                    worst = ratio.a
        rec(f"C.l{l}.discs_clear_E", clear, min_dist2_over_2rho2=mpmath.nstr(mpmath.mpf(worst), 6))
        # exact re-sum of the increments
        lo = sum((parse_ball(s["increment"])[0] for s in pcs), Fr(0))
        hi = sum((parse_ball(s["increment"])[1] for s in pcs), Fr(0))
        wn = ivhull(lo, hi) / pi_iv
        k = wd["winding"]
        rec(f"C.l{l}.resum_within_quarter", isinstance(k, int) and (abs(wn - k)).b < Fr(1, 4).__float__(),
            resum_over_pi=f"[{mpmath.nstr(mpmath.mpf(wn.a), 8)}, {mpmath.nstr(mpmath.mpf(wn.b), 8)}]")
        rec(f"C.l{l}.winding_value", k == (1 if l == 0 else 0), stored=k)
        stored_sum = ball_iv(wd["sum_increments_over_H"])
        rec(f"C.l{l}.resum_matches_stored_sum", overlaps(ivhull(lo, hi), stored_sum))
        per_l[l] = {"a_l": str(a_w), "a_l_over_km": mid((ivq(a_w) / km), 10), "pieces": len(pcs),
                    "winding": k}
    rec("C.all_l_0_19_present", all_present)
    rec("C.b_e_common", len(b_set) == 1 and len(e_set) == 1)
    out["l_le_19"] = per_l
    out["l_le_19_pieces_total"] = sum(v["pieces"] for v in per_l.values())

    # ---------------- D: l >= 20
    cov = json.loads((CERT / "part1" / "cover_PB.json").read_text())
    C = {k: Fr(v) for k, v in cov["constants"].items()}
    h0 = C["h0"]
    rec("D.h0_2_41", h0 == Fr(2, 41))
    # 1/(l + 1/2) <= h0  <=>  l >= 20, at the integers around the seam (and monotone in l)
    rec("D.l19_outside_l20_inside", (Fr(1) / (19 + Fr(1, 2)) > h0) and (Fr(1) / (20 + Fr(1, 2)) <= h0))
    OmP = pt["Omega_P"]
    # union of the Re-ranges for every h in (0, h0]: [0, 8 km] u [6 km, zSU OmP/h] u [zTL, zTU] OmP/h u [zA OmP/h, oo)
    # Overlaps: 6 < 8; zSU OmP/h >= zSU OmP/h0 > 8 km (worst case h = h0, since OmP/h decreases in h); zTL < zSU; zA < zTU
    top_S_h0 = ivq(C["zeta_SU"]) * OmP / ivq(h0)
    rec("D.regimes_union_is_half_line", C["Omega_S_over_km"] < C["Omega_L_over_km"]
        and (top_S_h0 - ivq(C["Omega_L_over_km"]) * km).a > 0 and C["zeta_TL"] < C["zeta_SU"]
        and C["zeta_A"] < C["zeta_TU"] and (ivq(C["Omega_S_over_km"]) * km - top_S_h0).b < 0,
        S_top_at_h0_over_km=mid((top_S_h0 / km), 8))
    rec("D.cover_flags", cov["pass"] is True and all(cov["certificates_pass"].values())
        and all(cov["inequalities"].values()) and all(cov["tiling"].values()))
    rp1 = json.loads((CERT / "part1" / "replay_PB.json").read_text())
    rec("D.part1_full_replay_recorded", rp1.get("pass") is True and rp1.get("stride") == 1)

    # ---------------- E: sampled cover of S_delta for every l
    km_lo, km_hi = iv_lo(km), iv_hi(km)
    d_lo, d_hi = iv_lo(d), iv_hi(d)
    om_lo, om_hi = iv_lo(OmP), iv_hi(OmP)
    a_tab = {l: Fr(per_l[l]["a_l"]) for l in per_l}

    def covers(l: int, x: Fr, y: Fr) -> str | None:
        """Which certified exclusion contains (l, w = x + i y); y >= -d_lo (inside S_delta for sure).
        Robust to the interval uncertainty of kappa_-, Omega_P: a regime counts only if it contains the
        point for every value in the balls."""
        if y > 0:
            return "P2-0" if not (l == 0 and x == 0 and y == 0) else None
        ax = abs(x)                                 # Re w < 0 by the conjugation symmetry
        if l <= 19:
            if ax <= a_tab[l] and -b <= y <= e:     # R'_l (count 1 at l = 0: the zero at 0; else 0)
                return "P2-3B" if not (l == 0 and ax == 0 and y == 0) else "zero-at-0"
            if ax >= a_tab[l] and 0 >= y:            # y >= -d_lo >= -d_delta by construction
                return "P2-2"
            return None
        h = Fr(1) / (l + Fr(1, 2))
        if not h <= h0:
            return None
        if ax <= C["Omega_L_over_km"] * km_lo:
            return "P2-1:L"
        if C["Omega_S_over_km"] * km_hi <= ax <= C["zeta_SU"] * om_lo / h:
            return "P2-1:S"
        if C["zeta_TL"] * om_hi / h <= ax <= C["zeta_TU"] * om_lo / h:
            return "P2-1:T"
        if ax >= C["zeta_A"] * om_hi / h:
            return "P2-1:A"
        return None

    rng = random.Random(20261026)
    pts = []
    # seam points
    for l in (0, 1, 18, 19, 20, 21, 22):
        for y in (Fr(0), -d_lo, -d_lo / 2, e, Fr(1, 10 ** 12)):
            xs = [Fr(0), km_lo * 8, km_hi * 6, km_hi * 8 + Fr(1, 10 ** 9)]
            if l <= 19:
                xs += [a_tab[l], -a_tab[l], a_tab[l] + Fr(1, 2 ** 40), a_tab[l] - Fr(1, 2 ** 40)]
            else:
                h = Fr(1) / (l + Fr(1, 2))
                xs += [C["zeta_SU"] * om_lo / h, C["zeta_TL"] * om_hi / h, C["zeta_TU"] * om_lo / h,
                       C["zeta_A"] * om_hi / h, 10 ** 6 * km_hi]
            for x in xs:
                pts.append((l, x, y))
    # random points: small and large l, Re w up to 3 (l+1/2) Omega_P, Im w in [-d, e]
    for _ in range(20000):
        l = rng.choice([rng.randint(0, 25), rng.randint(0, 200), rng.randint(20, 10 ** 6)])
        xmax = 3 * (l + Fr(1, 2)) * om_hi + 10 * km_hi
        x = Fr(rng.randint(-10 ** 9, 10 ** 9), 10 ** 9) * xmax
        y = Fr(rng.randint(-10 ** 9, 10 ** 9), 10 ** 9) * d_lo
        if rng.random() < 0.1:
            y = abs(y) * Fr(2, 10)                   # some points in the band 0 < Im w <= e and above
        pts.append((l, x, y))
    tally: dict[str, int] = {}
    unc = []
    for (l, x, y) in pts:
        c_ = covers(l, x, y)
        if c_ is None:
            unc.append((l, str(x), str(y)))
        tally[c_ or "UNCOVERED"] = tally.get(c_ or "UNCOVERED", 0) + 1
    rec("E.every_sampled_point_covered", not unc, n_points=len(pts), tally=tally, uncovered=unc[:10])

    # ---------------- F: beta, D4
    kz = json.loads((CERT / "part3" / "krawczyk_PB_l0.json").read_text())
    im = ball_iv(kz["enclosure_over_kappa_minus"]["im"])
    rec("F.krawczyk_flags", kz.get("certified_unique_zero") is True and kz.get("l") == 0 and kz.get("point") == "PB")
    rec("F.beta_upper_0.90938098", (-im - ivq(Fr(90938098, 10 ** 8))).b <= 0,
        minus_im_over_km=mid((-im), 20))
    rec("F.omega_star_below_strip", (-im - ivq(dfrac)).a > 0)
    # off E_+ u E_c: |Im w_*| / km vs n kappa / km for n = 1, 2 (kappa_+), n = 1 (kappa_c); Re part ~ 0
    offE = all((abs(-im - n * kap / km) - ivq(Fr(1, 100))).a > 0 for kap, ns in ((kp, (1, 2)), (kc, (1,))) for n in ns)
    rec("F.omega_star_off_E", offE, kplus_over_km=mid((kp / km), 10),
        kc_over_km=mid((kc / km), 10))
    dW = (rp ** 2 + rc ** 2) / (rp * rc)
    rec("F.D4_W0prime_nonzero", dW.a > 0, W0prime_over_i=mid(dW, 12))

    out["checks"] = checks
    out["failed"] = [k for k, v in checks.items() if not v["ok"]]
    out["pass"] = not out["failed"]
    out["n_checks"] = len(checks)
    out["timing_s"] = round(time.time() - t_start, 2)
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(f"independent assembly check: pass = {out['pass']}  ({len(checks)} checks, failed: {out['failed']})")
    print(f"  b - d_delta = {checks['B.b_minus_d_delta_pos']['margin']} (lower {checks['B.b_minus_d_delta_pos']['margin_lower']},"
          f" d_delta width {checks['B.b_minus_d_delta_pos']['width_of_d_delta']})")
    print(f"  cover sample: {checks['E.every_sampled_point_covered']['tally']}")
    print(f"  wrote {OUT.name} in {out['timing_s']} s")
    return 0 if out["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
