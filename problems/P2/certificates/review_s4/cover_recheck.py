"""Independent re-check of the S4 cover (W9) from the stored certificates (S4 review).  No p2 code.

Exact rationals for every stored endpoint; P_B constants from indep_params (mpmath, 50 digits).
Checks, for every l >= 20 (h = 1/(l+1/2) in (0, 2/41]) and every w in H+ = {Re w >= 0, -d <= Im w <= 0}:
  * T: the 24 zeta-boxes tile [4/5, 2] (Omega_P units), nu range [-d, 0], h0 = 2/41, all pass, decay flags;
  * A: the box at zeta_A = 3/2 passes, W7 flags; phi0 = -1/5 < -atan(d h0/zeta_A) (independent);
  * S (each side): the zeta-leaves tile [0, 9/10]; in each leaf the eps-boxes tile [0, eps_max] and
    eps_max >= min(1/(6 kappa_-), h0/zeta_lo) (independent constants); exit, order, eta/phi-at-exit and W1 flags;
    side +: X_a < X0 and X_L < X0; side c: X_L > 0 and each crossing height inside S2's window;
  * L: the 1760 leaves tile [0, 8] x [-11/20, 0] exactly (area and pairwise disjointness), each passes both sides;
  * the regime overlaps: 6 < 8; 4/5 < 9/10; 3/2 < 2; (9/10) Omega_P/h0 > 8 kappa_- (S non-empty and above L);
    d < kappa_+ < kappa_c; 6 kappa_- > d, (4/5) Omega_P/h0 > d, (3/2) Omega_P/h0 > d.
Then, for 2000 random (l, w) with l in [20, 10^6] and Re w up to 50 Omega_P/h, it locates a covering certificate
piece by exact lookup (the K3 idea, done independently and for many points).
"""
from __future__ import annotations

import json
import random
import sys
from fractions import Fraction as Fr
from itertools import pairwise
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from indep_params import PB  # noqa: E402

CERT = Path(__file__).resolve().parents[1] / "part1"
p = PB()
h0 = Fr(2, 41)


def tiles(iv, lo, hi):
    iv = sorted(iv)
    return bool(iv) and iv[0][0] == lo and iv[-1][1] == hi and all(a[1] == b[0] for a, b in pairwise(iv))


def main():
    out = {}
    T = json.loads((CERT / "T_PB.json").read_text())
    A = json.loads((CERT / "A_PB.json").read_text())
    out["T_tiles_[4/5,2]"] = tiles([(Fr(b["zeta_over_OmegaP"][0]), Fr(b["zeta_over_OmegaP"][1])) for b in T["boxes"]], Fr(4, 5), Fr(2))
    out["T_all_pass_h0_decay"] = all(b["pass"] and b["h0"] == "2/41" and b["jost_decay"] for b in T["boxes"]) and T["pass"]
    out["T_nu_low_is_minus_d"] = all(abs(mp.mpf(b["nu_low"].split("+/-")[0].strip("[ ")) + p.d) < 1e-9 for b in T["boxes"])
    out["A_pass"] = A["pass"] and A["boxes"][0]["zeta_over_OmegaP"] == ["3/2", "3/2"] and A["monotonicity_W7"]["pass"]
    out["A_phi0_lt_minus_atan"] = bool(mp.mpf(-1) / 5 < -mp.atan(p.d * mp.mpf(2) / 41 / (mp.mpf(3) / 2 * p.OmegaP)))
    Sres = {}
    crossing_max = {"X0": 0.0, "0": 0.0}
    dist_exit_max = 0.0
    for side, tag in ((1, "plus"), (-1, "c")):
        S = json.loads((CERT / f"S_PB_{tag}.json").read_text())
        leaves = S["leaves"]
        ok_tiles = tiles([(Fr(x["zeta_over_OmegaP"][0]), Fr(x["zeta_over_OmegaP"][1])) for x in leaves], Fr(0), Fr(9, 10))
        ok_eps, ok_flags, ok_w1 = True, True, True
        for x in leaves:
            zl = Fr(x["zeta_over_OmegaP"][0])
            em = Fr(x["eps_max"])
            need = 1 / (6 * p.km) if zl == 0 else min(1 / (6 * p.km), mp.mpf(2) / 41 / (mp.mpf(zl.numerator) / zl.denominator * p.OmegaP))
            ok_eps = ok_eps and mp.mpf(em.numerator) / em.denominator >= need
            ex = x["exit"]
            ok_flags = ok_flags and x["pass"] and ex["exit_forbidden_ok"] and ex["exit_order_ok"] and ex["phi_zero_at_exit"] \
                and ex["eta_zero_at_exit"] and x["w1_check"]["pass"]
            if x.get("closing"):
                continue
            ivs = [(Fr(e["eps"][0]), Fr(e["eps"][1])) for e in x["eps_boxes"]]
            ok_eps = ok_eps and tiles(ivs, Fr(0), em) and all(e["pass"] for e in x["eps_boxes"])
            w, ya = Fr(x["design"][0]), Fr(x["design"][1])
            ta = w + ya
            Xa, k = Fr(x["X_a"]), Fr(x["k"])
            Xa_m, k_m = mp.mpf(Xa.numerator) / Xa.denominator, mp.mpf(k.numerator) / k.denominator
            XL = Xa_m + side * (-8 - mp.mpf(ta.numerator) / ta.denominator) / k_m
            if side > 0:
                ok_w1 = ok_w1 and Xa_m < p.X0 and XL < p.X0
            else:
                ok_w1 = ok_w1 and XL > 0
                for X, lim, name in ((p.X0, mp.mpf("0.1413"), "X0"), (mp.mpf(0), mp.mpf("2.438"), "0")):
                    if X <= Xa_m:
                        continue
                    ts = mp.mpf(ta.numerator) / ta.denominator - k_m * (X - Xa_m)
                    wm = mp.mpf(w.numerator) / w.denominator
                    if ts <= -wm:
                        eta = ts
                    elif ts >= wm:
                        eta = mp.mpf(0)
                    else:
                        u = -ts / wm
                        eta = -wm * (mp.mpf(5) / 32 + u / 2 + mp.mpf(15) / 32 * u ** 2 - mp.mpf(5) / 32 * u ** 4 + u ** 6 / 32)
                    hgt = abs(eta / k_m)
                    ok_w1 = ok_w1 and hgt < lim
                    crossing_max[name] = max(crossing_max[name], float(hgt))
                    dist_exit_max = max(dist_exit_max, float(k_m * (X - Xa_m)))
        Sres[tag] = {"tiles_[0,9/10]": ok_tiles, "eps_tiles_and_eps_max": bool(ok_eps), "flags": bool(ok_flags),
                     "w1_and_tail_strip": bool(ok_w1), "n_leaves": len(leaves), "cert_pass": S["pass"]}
    out["S"] = Sres
    out["S_c_max_crossing_height"] = crossing_max
    out["S_c_max_t_distance_of_crossing_from_exit"] = dist_exit_max
    Lp = json.loads((CERT / "L_PB_plus.json").read_text())["leaves_detail"]
    Lc = json.loads((CERT / "L_PB_c.json").read_text())["leaves_detail"]
    boxes = [tuple(Fr(v) for v in x["box"]) for x in Lp]
    area = sum((b[1] - b[0]) * (b[3] - b[2]) for b in boxes)
    inside = all(b[0] >= 0 and b[1] <= 8 and b[2] >= Fr(-11, 20) and b[3] <= 0 for b in boxes)
    # disjoint: all boxes are on the 1/20 grid (none refined) -> distinct keys suffice; check that
    grid = all(b[1] - b[0] == Fr(1, 20) and b[3] - b[2] == Fr(1, 20) and (b[0] * 20).denominator == 1 and (b[2] * 20).denominator == 1 for b in boxes)
    distinct = len(set(boxes)) == len(boxes)
    out["L"] = {"n": len(boxes), "area_eq_8x11/20": area == 8 * Fr(11, 20), "inside": inside, "on_1/20_grid": grid,
                "distinct": distinct, "same_boxes_both_sides": [x["box"] for x in Lp] == [x["box"] for x in Lc],
                "all_pass": all(x["y_sign_ok"] and x["exit_rule_ok"] for x in Lp + Lc),
                "poles_avoided_im_ge_minus_d_gt_minus_kplus": bool(inside and p.d < p.kp)}
    km, OmP, d = p.km, p.OmegaP, p.d
    h0m = mp.mpf(2) / 41
    out["overlaps"] = {"6<8": True, "4/5<9/10": True, "3/2<2": True,
                       "0.9 OmP/h0 > 8 km": bool(mp.mpf("0.9") * OmP / h0m > 8 * km),
                       "d < kp < kc": bool(d < p.kp < p.kc),
                       "6 km > d": bool(6 * km > d), "0.8 OmP/h0 > d": bool(mp.mpf("0.8") * OmP / h0m > d),
                       "1.5 OmP/h0 > d": bool(mp.mpf("1.5") * OmP / h0m > d)}
    # random lookups
    rng = random.Random(7)
    Sl = {tag: json.loads((CERT / f"S_PB_{tag}.json").read_text())["leaves"] for tag in ("plus", "c")}
    Lset = set(boxes)
    miss = []
    counts = {"L": 0, "S": 0, "T": 0, "A": 0}
    for _ in range(2000):
        l = int(20 * 10 ** (rng.random() * 4.7))
        h = Fr(2, 2 * l + 1)
        hm = mp.mpf(2) / (2 * l + 1)
        re = mp.mpf(rng.random()) ** 3 * 50 * OmP / hm
        if rng.random() < 0.2:
            re = mp.mpf(rng.random()) * 9 * km
        im = -d * mp.mpf(rng.random())
        zeta = hm * re
        covered = []
        if re <= 8 * km:
            rk, ik = re / km, im / km
            i0, j0 = int(mp.floor(rk * 20)), int(mp.floor(-ik * 20))
            cands = [(Fr(i, 20), Fr(i + 1, 20), Fr(-(j + 1), 20), Fr(-j, 20)) for i in (i0 - 1, i0) for j in (j0 - 1, j0)]
            if any(c in Lset and c[0] * km <= re <= c[1] * km and c[2] * km <= im <= c[3] * km for c in cands):
                covered.append("L")
        if 6 * km <= re and zeta <= mp.mpf("0.9") * OmP:
            eps = 1 / re
            okS = True
            for tag in ("plus", "c"):
                hit = False
                for x in Sl[tag]:
                    a, b = (mp.mpf(Fr(v).numerator) / Fr(v).denominator * OmP for v in x["zeta_over_OmegaP"])
                    if a <= zeta <= b and x["pass"]:
                        if x.get("closing"):
                            hit = hit or eps <= mp.mpf(Fr(x["eps_max"]).numerator) / Fr(x["eps_max"]).denominator
                        else:
                            hit = hit or any(mp.mpf(Fr(e["eps"][0]).numerator) / Fr(e["eps"][0]).denominator <= eps <=
                                             mp.mpf(Fr(e["eps"][1]).numerator) / Fr(e["eps"][1]).denominator for e in x["eps_boxes"])
                okS = okS and hit
            if okS:
                covered.append("S")
        if mp.mpf("0.8") * OmP <= zeta <= 2 * OmP:
            covered.append("T")
        if zeta >= mp.mpf("1.5") * OmP:
            covered.append("A")
        for c in covered:
            counts[c] += 1
        if not covered:
            miss.append([l, float(re / km), float(im / km)])
    out["random_lookups"] = {"n": 2000, "uncovered": miss[:10], "n_uncovered": len(miss), "by_regime": counts}
    out["pass"] = bool(out["T_tiles_[4/5,2]"] and out["T_all_pass_h0_decay"] and out["T_nu_low_is_minus_d"] and out["A_pass"]
                       and out["A_phi0_lt_minus_atan"] and all(all(v for k, v in s.items() if k != "n_leaves") for s in Sres.values())
                       and all(v for k, v in out["L"].items() if k != "n") and all(out["overlaps"].values()) and not miss)
    print(json.dumps(out, indent=1))
    Path(sys.argv[1] if len(sys.argv) > 1 else "cover_recheck.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
