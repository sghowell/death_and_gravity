"""S3 review (independent): static consistency of the stored part-(iii) certificates.

Written by the reviewer; imports nothing from p2.  Parameters are recomputed from the exact
definitions of FORMULATION section 1.3 with mpmath (80 digits); all geometry is exact rational
arithmetic (fractions.Fraction); the increment balls are summed exactly as rational intervals
(Arb's printed ball [m +/- r] contains the computed ball, so [m - r, m + r] is a valid interval).

Checks, for every l in 0..19 (certificates/part3/winding_PB_l{l}.json):
  P1  the region: a equals W_l_exact of certificates/lg/lg_PB_l{l}.json; b >= (1/2 + 1/20) kappa_-;
      0 < e <= kappa_-/10; b < kappa_+;
  P2  the pieces tile H in order: start (0, -b), corners (a, -b) and (a, e), end (0, e); every
      piece axis-parallel, non-degenerate, oriented along H, with exact shared endpoints;
  P3  every piece lies within rho/2 of its centre, and the centre lies on the piece's line;
  P4  the box [c +- rho]^2 avoids E_+ u E_c (sqrt(2) rho < distance from c to every -i n kappa);
  P5  every segment increment ball has |inc| < (number of pieces) * pi (each piece is certified
      below pi when computed; the pieces themselves are not stored), and the exact rational sum / pi is within 1/4 of an
      integer (reported with its distance), equal to the stored winding, and equal to the
      expected count (1 at l = 0, else 0);
  P6  the stored totals (segments, K_max, sum, winding ball) agree with the pieces.
Also: the Krawczyk certificate (C6) and the controls file (C0 rel. differences, C2/C3 counts).
Output: check_certs.json next to this file.
"""
from __future__ import annotations

import json
import re
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
CERT = HERE.parent / "part3"
LG = HERE.parent / "lg"
mp.mp.dps = 80

# pi as an exact rational interval
PI_LO = Fr("3.14159265358979323846264338327950288419716939937510")
PI_HI = PI_LO + Fr(1, 10 ** 49)


def params(qfrac=Fr(999, 1000)):
    """Roots and surface gravities of Delta = -(Lam/3) r^4 + r^2 - 2 r + Q^2 at Lam = 3/50."""
    # r_ext: smallest positive root of r^3 - 25 r + 25 = 0
    cands = [mp.findroot(lambda r: r ** 3 - 25 * r + 25, x0) for x0 in (1.0, 1.2, 4.0)]
    rext = min(c for c in cands if c > 0)
    qmax2 = rext - rext ** 4 / 50
    Q2 = mp.mpf(qfrac.numerator) ** 2 / mp.mpf(qfrac.denominator) ** 2 * qmax2
    co = [-mp.mpf(1) / 50, 0, 1, -2, Q2]                     # descending powers
    roots = sorted(mp.re(x) for x in mp.polyroots(co, maxsteps=200, extraprec=200))
    rn, rm, rp, rc = roots
    dD = lambda r: -4 * r ** 3 / 50 + 2 * r - 2
    kap = {k: abs(dD(r)) / (2 * r * r) for k, r in (("-", rm), ("+", rp), ("c", rc))}
    return {"rn": rn, "rm": rm, "rp": rp, "rc": rc, "km": kap["-"], "kp": kap["+"], "kc": kap["c"],
            "Q2": Q2, "rext": rext}


BALL = re.compile(r"^\[\s*(?:([-+0-9.eE]+)\s*)?\+/-\s*([0-9.eE+-]+)\s*\]$")


def ball(s: str):
    """Rational interval [lo, hi] from an Arb ball string."""
    s = s.strip()
    m = BALL.match(s)
    if m:
        mid = Fr(m.group(1)) if m.group(1) else Fr(0)
        rad = Fr(m.group(2))
        return mid - rad, mid + rad
    if s.startswith("[") and s.endswith("]"):
        v = Fr(s[1:-1])
        return v, v
    v = Fr(s)
    return v, v


def fq(s: str) -> Fr:
    return Fr(s)


def check_l(l: int, P: dict) -> dict:
    st = json.loads((CERT / f"winding_PB_l{l}.json").read_text())
    lg = json.loads((LG / f"lg_PB_l{l}.json").read_text())
    probs = []
    reg = st["region"]
    a, b, e = fq(reg["a"]), fq(reg["b"]), fq(reg["e"])
    km, kp, kc = P["km"], P["kp"], P["kc"]
    if a != fq(lg["W_l_exact"]):
        probs.append("P1: a != lg W_l_exact")
    if lg.get("point") != "PB" or lg.get("delta") != "1/20" or lg.get("qfrac") != "999/1000":
        probs.append("P1: lg certificate parameters")
    if not (mp.mpf(b.numerator) / b.denominator >= mp.mpf("0.55") * km):
        probs.append("P1: b < (1/2+delta) kappa_-")
    if not (0 < e and mp.mpf(e.numerator) / e.denominator <= km / 10):
        probs.append("P1: e not in (0, eta]")
    if not (mp.mpf(b.numerator) / b.denominator < kp):
        probs.append("P1: b >= kappa_+")
    if st.get("point") != "PB" or st.get("l") != l:
        probs.append("P1: point/l label")
    pieces = st["pieces"]
    # P2 tiling of H
    corners = [(Fr(0), -b), (a, -b), (a, e), (Fr(0), e)]
    cur = corners[0]
    edge = 0
    for r in pieces:
        p0 = tuple(fq(x) for x in r["p0"])
        p1 = tuple(fq(x) for x in r["p1"])
        if p0 != cur:
            probs.append(f"P2: gap before {r['key']}")
            break
        z0, z1 = corners[edge], corners[edge + 1]
        # the piece must lie on edge `edge`, in its direction, and not past its end
        if z0[1] == z1[1]:          # horizontal edge
            ok = p0[1] == p1[1] == z0[1] and (p1[0] - p0[0]) * (z1[0] - z0[0]) > 0 \
                and min(z0[0], z1[0]) <= p1[0] <= max(z0[0], z1[0])
        else:                       # vertical edge
            ok = p0[0] == p1[0] == z0[0] and (p1[1] - p0[1]) * (z1[1] - z0[1]) > 0 \
                and min(z0[1], z1[1]) <= p1[1] <= max(z0[1], z1[1])
        if not ok or r["edge"] != edge:
            probs.append(f"P2: piece {r['key']} not on edge {edge} in order")
            break
        cur = p1
        if cur == z1:
            edge += 1
            if edge == 3:
                break
    if cur != corners[3] or edge != 3 or r is not pieces[-1]:
        probs.append("P2: H not closed at (0, e) by the last piece")
    # P3, P4, P5
    lo_sum, hi_sum = Fr(0), Fr(0)
    max_inc = Fr(0)
    min_pole_margin = None
    for r in pieces:
        p0 = tuple(fq(x) for x in r["p0"])
        p1 = tuple(fq(x) for x in r["p1"])
        c = tuple(fq(x) for x in r["c"])
        rho = fq(r["rho"])
        for p in (p0, p1):
            if (p[0] - c[0]) ** 2 + (p[1] - c[1]) ** 2 > (rho / 2) ** 2:
                probs.append(f"P3: {r['key']} endpoint farther than rho/2 from c")
        if not ((p0[0] == p1[0] == c[0]) or (p0[1] == p1[1] == c[1])):
            probs.append(f"P3: centre of {r['key']} off the piece line")
        cx = mp.mpf(c[0].numerator) / c[0].denominator
        cy = mp.mpf(c[1].numerator) / c[1].denominator
        dist = min(mp.sqrt(cx ** 2 + (cy + n * k) ** 2) for k in (kp, kc) for n in range(1, 200))
        s = mp.sqrt(2) * mp.mpf(rho.numerator) / rho.denominator
        if not dist > s:
            probs.append(f"P4: box of {r['key']} meets E")
        marg = dist / s
        min_pole_margin = marg if min_pole_margin is None else min(min_pole_margin, marg)
        lo, hi = ball(r["increment"])
        if not (max(abs(lo), abs(hi)) < r["pieces"] * PI_LO):
            probs.append(f"P5: |increment| not < pieces * pi at {r['key']}")
        if r["pieces"] < 1:
            probs.append(f"P5: no pieces at {r['key']}")
        max_inc = max(max_inc, abs(lo), abs(hi))
        lo_sum += lo
        hi_sum += hi
    # winding interval: sum / pi
    cands = [lo_sum / PI_LO, lo_sum / PI_HI, hi_sum / PI_LO, hi_sum / PI_HI]
    wlo, whi = min(cands), max(cands)
    k = round((wlo + whi) / 2)
    dist_int = max(abs(wlo - k), abs(whi - k))
    expected = 1 if l == 0 else 0
    if not dist_int < Fr(1, 4):
        probs.append("P5: winding interval not within 1/4 of an integer")
    if k != st["winding"] or k != expected or not st["certified"]:
        probs.append(f"P5: winding {k} vs stored {st['winding']} vs expected {expected}")
    # P6 stored aggregates
    if st["segments"] != len(pieces):
        probs.append("P6: segment count")
    if st["K_max_used"] != max(r["K"] for r in pieces):
        probs.append("P6: K_max")
    slo, shi = ball(st["sum_increments_over_H"])
    if not (slo <= hi_sum and lo_sum <= shi):
        probs.append("P6: stored sum does not overlap the re-sum")
    return {"l": l, "a_over_km": float(mp.mpf(a.numerator) / a.denominator / km),
            "b_over_km": float(mp.mpf(b.numerator) / b.denominator / km),
            "e_over_km": float(mp.mpf(e.numerator) / e.denominator / km),
            "segments": len(pieces), "K_max": max(r["K"] for r in pieces),
            "max_M_over_absWc": max(r["M_over_absWc"] for r in pieces),
            "winding_interval": [float(wlo), float(whi)], "winding": k,
            "max_distance_from_integer": float(dist_int), "max_abs_increment": float(max_inc),
            "min_pole_margin_ratio": float(min_pole_margin),
            "min_abs_W_lower_stored": min(r["min_abs_enclosure_lower"] for r in pieces),
            "computed_at": st.get("computed_at"), "problems": probs}


def main():
    P = params()
    out = {"params_independent": {k: mp.nstr(v, 30) for k, v in P.items()},
           "params_over_km": {"kp/km": mp.nstr(P["kp"] / P["km"], 15), "kc/km": mp.nstr(P["kc"] / P["km"], 15)},
           "per_l": []}
    allp = []
    for l in range(20):
        r = check_l(l, P)
        out["per_l"].append(r)
        allp += [f"l={l}: {p}" for p in r["problems"]]
    out["total_segments"] = sum(r["segments"] for r in out["per_l"])
    out["max_distance_from_integer"] = max(r["max_distance_from_integer"] for r in out["per_l"])
    # Krawczyk certificate and controls
    kr = json.loads((CERT / "krawczyk_PB_l0.json").read_text())
    ilo, ihi = ball(kr["enclosure_over_kappa_minus"]["im"])
    rlo, rhi = ball(kr["enclosure_over_kappa_minus"]["re"])
    k6lo, k6hi = ball(kr["KAT6_enclosure_im"])
    out["krawczyk"] = {"certified_unique_zero": kr["certified_unique_zero"], "contraction": kr["contraction"],
                       "im_interval": [str(ilo), str(ihi)], "re_interval": [str(rlo), str(rhi)],
                       "beta_upper_from_enclosure": float(-ilo), "beta_le_0.90938098": bool(-ilo <= Fr("0.90938098")),
                       "overlaps_KAT6": bool(ilo <= k6hi and k6lo <= ihi),
                       "below_kappa_plus": bool(mp.mpf(-ilo.numerator) / ilo.denominator > P["kp"] / P["km"]),
                       "above_2kappa_plus": bool(mp.mpf(-ihi.numerator) / ihi.denominator < 2 * P["kp"] / P["km"])}
    if not (kr["certified_unique_zero"] and kr["contraction"] < 1):
        allp.append("Krawczyk not certified")
    ctl = json.loads((CERT / "controls_PB.json").read_text())
    out["controls"] = {"commit": ctl["commit"],
                       "C0_max_rel_diff": max(x["rel_diff_mid"] for x in ctl["C0_zsol_vs_arbw"]),
                       "C0_all_pass": all(x["pass"] for x in ctl["C0_zsol_vs_arbw"]),
                       "C0_arbw_prec": sorted({x["arbw_prec"] for x in ctl["C0_zsol_vs_arbw"]}),
                       "C2": [(x["point"], x["l"], x["box_kappa"], x["winding"], x["expected"]) for x in ctl["C2"]],
                       "C3": [(x["point"], x["l"], x["box_kappa"], x["winding"], x["expected"]) for x in ctl["C3"]],
                       "C1": ctl["C1_D4"]["pass"], "C6": ctl["C6_pass"], "pass": ctl["pass"]}
    out["problems"] = allp
    out["verdict"] = "OK" if not allp else "PROBLEMS"
    (HERE / "check_certs.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: out[k] for k in ("total_segments", "max_distance_from_integer", "krawczyk", "controls",
                                          "problems", "verdict")}, indent=1))
    for r in out["per_l"]:
        print(r["l"], round(r["a_over_km"], 4), r["segments"], r["K_max"], f"{r['max_M_over_absWc']:.1e}",
              r["winding"], f"{r['max_distance_from_integer']:.2e}", round(r["min_pole_margin_ratio"], 3),
              f"{r['min_abs_W_lower_stored']:.4g}", r["computed_at"])


if __name__ == "__main__":
    main()
