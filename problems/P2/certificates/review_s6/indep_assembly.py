"""S6 review, checks C8-C11: the seams and the cover of S_e, re-checked with the reviewer's own code.

P_B comes from its definitions (indep_rnds.Point, Arb 256 bits); nothing from p2 is imported.  The
certificates are read as data.  Checks:
  S1  the E-points in S_e: only omega0 = -i kappa_+ (n kappa_+ <= e kappa_- only for n = 1; kappa_c > e kappa_-);
  S2  per band certificate (l <= 19): b' = ceil(e kappa_- 2^30)/2^30 (least dyadic >= e kappa_-), b0 = floor(kappa_-/2 2^30)/2^30,
      b0 < 0.55 kappa_- (overlap with S_delta); a'_l = A_fixed's a_prime, a'_l >= zeta_A Omega_P/h_l and a'_l is the
      least 2^-30 dyadic above it; omega0 interior; the stored pieces trace the right half of dB_l exactly; each piece
      within rho/2 of its centre; each box [c +- rho]^2 clear of every point of E_+ u E_c (all n); the exact re-sum of the
      increments over pi lies within 1/4 of the stored winding, which is 0 (l = 0) or -1;
  S3  omega_* (P2-K6's ball) inside B_0, and its distance to the bottom edge;
  S4  part (ii'): W7 at fixed h_l (arctan(e kappa_- h_l/(zeta_A Omega_P)) < 1/5), Jost decay, nu_low <= -e kappa_-;
  S5  part (i) at depth e: the cover inequalities (recomputed), the regime tilings read from the certificates;
  S6  a sampled cover: 40 000 random (l, omega) in S_e (l log-uniform up to 10^6, Re omega of both signs, the seams
      over-sampled), each assigned to a certified region by the decision rule of s6-proof.md section 6.

Run:  python indep_assembly.py      -> indep_assembly_out.json
"""
from __future__ import annotations

import json
import math
import random
from fractions import Fraction
from pathlib import Path

from flint import arb, ctx

import indep_rnds as R

HERE = Path(__file__).resolve().parent
CERT = HERE.parent
E = Fraction(23, 25)
DD = Fraction(11, 20)
H0 = Fraction(2, 41)
DEN = 2 ** 30


def A(x):
    return R.A(x)


def nu_ok(s, e) -> bool:
    """the stored nu_low (a 10-digit ball string of the lower endpoint used by the checker) reaches -e kappa_-:
    the ball contains -e kappa_- or lies below it"""
    b = arb(s)
    return bool(b.contains(-e)) or bool(b < -e)


class Rec:
    def __init__(self):
        self.items = {}

    def add(self, name, ok, note=None):
        self.items[name] = {"ok": bool(ok), "note": note}
        return bool(ok)


def main():
    ctx.prec = 256
    pt = R.Point()
    km, kp, kc, OmP = pt.km, pt.kp, pt.kc, pt.OmegaP
    e = A(E) * km
    Rr = Rec()
    out = {"point": pt.summary()}
    # S1
    n_in = [n for n in range(1, 10) if not bool(n * kp > e)]
    Rr.add("S1.E_plus_in_S_e_only_n1", n_in == [1] and bool(kp < e) and bool(2 * kp > e),
           note=f"kappa_+ = {(kp / km).str(12)} kappa_-, 2 kappa_+ = {(2 * kp / km).str(12)}")
    Rr.add("S1.E_c_below_S_e", bool(kc > e), note=f"kappa_c - e kappa_- = {((kc - e) / km).str(10)} kappa_-")
    # S2
    afix = {r["l"]: r for r in json.loads((CERT / "beta" / "part2" / "A_fixed_PB.json").read_text())["rows"]}
    bp_ok = Fraction(math.ceil(float((e * DEN).mid())), DEN)
    while not bool(A(bp_ok) >= e):
        bp_ok += Fraction(1, DEN)
    while bool(A(bp_ok - Fraction(1, DEN)) >= e):
        bp_ok -= Fraction(1, DEN)
    b0_ok = Fraction(math.floor(float((km / 2 * DEN).mid())), DEN)
    while not bool(A(b0_ok) <= km / 2):
        b0_ok -= Fraction(1, DEN)
    while bool(A(b0_ok + Fraction(1, DEN)) <= km / 2):
        b0_ok += Fraction(1, DEN)
    out["b_prime"], out["b0"] = str(bp_ok), str(b0_ok)
    out["b_prime_minus_e_kappa"] = (A(bp_ok) - e).str(5)
    out["e_kappa_enclosure_radius"] = float(e.rad())
    sqrt2 = arb(2).sqrt()
    per_l = {}
    for l in range(20):
        bd = json.loads((CERT / "beta" / "part3" / f"band_PB_l{l}.json").read_text())
        reg = bd["region"]
        a, bp, b0 = Fraction(reg["a_prime"]), Fraction(reg["b_prime"]), Fraction(reg["b0"])
        Rr.add(f"S2.l{l}.b_prime_least_dyadic_ge_e", bp == bp_ok)
        Rr.add(f"S2.l{l}.b0_greatest_dyadic_le_half", b0 == b0_ok and b0 < DD * 1 and bool(A(b0) < A(DD) * km))
        za = Fraction(afix[l]["zeta_A_over_OmegaP"])
        h = Fraction(2, 2 * l + 1)
        thr = A(za) * OmP / A(h)
        Rr.add(f"S2.l{l}.a_equals_A_fixed", a == Fraction(afix[l]["a_prime"]))
        Rr.add(f"S2.l{l}.a_ge_zetaA_OmegaP_over_h", bool(A(a) >= thr) and bool(A(a - Fraction(1, DEN)) < thr),
               note=f"a'_l/kappa_- = {(A(a) / km).str(8)}")
        Rr.add(f"S2.l{l}.omega0_interior", bool(A(b0) < kp) and bool(kp < A(bp)) and a > 0)
        pcs = bd["pieces"]
        # the path
        pts = [(Fraction(p["p0"][0]), Fraction(p["p0"][1])) for p in pcs] + \
              [(Fraction(pcs[-1]["p1"][0]), Fraction(pcs[-1]["p1"][1]))]
        cont = all((Fraction(pcs[i]["p1"][0]), Fraction(pcs[i]["p1"][1])) ==
                   (Fraction(pcs[i + 1]["p0"][0]), Fraction(pcs[i + 1]["p0"][1])) for i in range(len(pcs) - 1))
        on_edges = True
        stage = 0
        for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
            if y0 == y1 == -bp and 0 <= x0 < x1 <= a and stage == 0:
                pass
            elif x0 == x1 == a and -bp <= y0 < y1 <= -b0 and stage <= 1:
                stage = 1
            elif y0 == y1 == -b0 and a >= x0 > x1 >= 0 and stage <= 2:
                stage = 2
            else:
                on_edges = False
        Rr.add(f"S2.l{l}.path", cont and on_edges and pts[0] == (0, -bp) and pts[-1] == (0, -b0) and stage == 2
               and any(y == -bp for _, y in pts) and any(x == a for x, _ in pts))
        # discs and pole clearance
        disc_ok, pole_ok, worst = True, True, None
        for p in pcs:
            cx, cy = Fraction(p["c"][0]), Fraction(p["c"][1])
            rho = Fraction(p["rho"])
            for q in ("p0", "p1"):
                x, y = Fraction(p[q][0]), Fraction(p[q][1])
                if (x - cx) ** 2 + (y - cy) ** 2 > (rho / 2) ** 2:
                    disc_ok = False
            rr = sqrt2 * A(rho)
            for kap in (kp, kc):
                nmax = int(float(((abs(A(cy)) + 2 * rr) / kap).upper())) + 2
                for n in range(1, nmax + 1):
                    d = ((A(cx)) ** 2 + (A(cy) + n * kap) ** 2).sqrt()
                    ratio = d / rr
                    if not bool(d > rr):
                        pole_ok = False
                    if worst is None or float(ratio.lower()) < worst:
                        worst = float(ratio.lower())
        Rr.add(f"S2.l{l}.discs", disc_ok)
        Rr.add(f"S2.l{l}.pole_clearance_all_n", pole_ok, note=f"min dist/(sqrt2 rho) = {worst:.4f}")
        ctx.prec = 256
        tot = sum((arb(p["increment"]) for p in pcs), arb(0))
        wn = tot / arb.pi()
        exp_w = 0 if l == 0 else -1
        Rr.add(f"S2.l{l}.resum", bool(abs(wn - exp_w) < arb("0.25")) and bd["winding"] == exp_w,
               note=f"{wn.str(6)}")
        per_l[l] = {"a_over_km": float((A(a) / km).mid()), "pieces": len(pcs), "resum_over_pi": wn.str(6),
                    "min_pole_clearance_ratio": worst, "K_max": max(p["K"] for p in pcs)}
    out["per_l"] = per_l
    # S3 omega_*
    kz = json.loads((CERT / "part3" / "krawczyk_PB_l0.json").read_text())["enclosure_over_kappa_minus"]
    im, re_ = arb(kz["im"]), arb(kz["re"])
    a0 = Fraction(json.loads((CERT / "beta" / "part3" / "band_PB_l0.json").read_text())["region"]["a_prime"])
    Rr.add("S3.omega_star_in_B0", bool(-im * km < A(bp_ok)) and bool(-im * km > A(b0_ok)) and bool(abs(re_) * km < A(a0)),
           note=f"distance to the bottom edge {((A(bp_ok) + im * km) / km).str(8)} kappa_-")
    Rr.add("S3.omega_star_off_E", bool(-im * km > kp) and bool(-im * km < 2 * kp) and bool(-im * km < kc))
    Rr.add("S3.beta_width", bool(2 * im.rad() < 1e-8), note=f"width {2 * float(im.rad()):.3e}")
    out["beta"] = (-im).str(32, radius=True)
    # S4 part (ii')
    for l in range(20):
        r = afix[l]
        h = Fraction(2, 2 * l + 1)
        za = Fraction(r["zeta_A_over_OmegaP"])
        bound = ((e * A(h)) / (A(za) * OmP)).atan()
        Rr.add(f"S4.l{l}.W7_and_decay", bool(bound < A(Fraction(1, 5))) and bool(A(za) * OmP / A(h) > e)
               and r["box"]["h0"] == str(h) and nu_ok(r["box"]["nu_low"], e)
               and r["box"]["jost_decay"] and r["box"]["pass"] and r["W7"]["pass"],
               note=f"arctan bound {bound.str(6)}")
    # S5 part (i)
    cv = json.loads((CERT / "beta" / "part1" / "cover_PB.json").read_text())
    C = {k: Fraction(v) for k, v in cv["constants"].items()}
    Rr.add("S5.overlaps", C["Omega_S_over_km"] < C["Omega_L_over_km"] and C["zeta_TL"] < C["zeta_SU"] and
           C["zeta_A"] < C["zeta_TU"])
    Rr.add("S5.S_top_above_L", bool(A(C["zeta_SU"]) * OmP / A(H0) > A(C["Omega_L_over_km"]) * km),
           note=f"{(A(C['zeta_SU']) * OmP / A(H0) / km - 8).str(6)} kappa_- above")
    Rr.add("S5.ReplusIm_positive", bool(A(C["Omega_S_over_km"]) * km > e) and bool(A(C["zeta_TL"]) * OmP / A(H0) > e)
           and bool(A(C["zeta_A"]) * OmP / A(H0) > e))
    Lc = json.loads((CERT / "beta" / "part1" / "L_PB.json").read_text())
    Lim = [Fraction(x) for x in Lc["design"]["Im_range_over_km"]]
    Rr.add("S5.L_rows", Lim[0] <= -E and Lim[1] >= -DD and Lc["pass"] and Lc["tiles_ok"] and Lc["grid_boxes"] == 1280
           and [Fraction(x) for x in Lc["design"]["h2_range"]] == [0, H0 * H0])
    T = json.loads((CERT / "beta" / "part1" / "T_PB.json").read_text())
    iv = sorted((Fraction(b["zeta_over_OmegaP"][0]), Fraction(b["zeta_over_OmegaP"][1])) for b in T["boxes"])
    Rr.add("S5.T_tiles", iv[0][0] == C["zeta_TL"] and iv[-1][1] == C["zeta_TU"] and
           all(x[1] == y[0] for x, y in zip(iv, iv[1:])) and all(b["pass"] for b in T["boxes"]) and
           all(nu_ok(b["nu_low"], e) for b in T["boxes"]) and all(b["h0"] == str(H0) for b in T["boxes"]))
    Ac = json.loads((CERT / "beta" / "part1" / "A_PB.json").read_text())
    bnd = ((e * A(H0)) / (A(C["zeta_A"]) * OmP)).atan()
    Rr.add("S5.A", Fraction(Ac["boxes"][0]["zeta_over_OmegaP"][0]) == C["zeta_A"] and Ac["pass"]
           and nu_ok(Ac["boxes"][0]["nu_low"], e) and bool(bnd < A(Fraction(1, 5))))
    for tag in ("plus", "c"):
        S = json.loads((CERT / "beta" / "part1" / f"S_PB_{tag}.json").read_text())
        iv = sorted((Fraction(x["zeta_over_OmegaP"][0]), Fraction(x["zeta_over_OmegaP"][1])) for x in S["leaves"])
        Rr.add(f"S5.S_{tag}_tiles", iv[0][0] == 0 and iv[-1][1] == C["zeta_SU"] and
               all(x[1] == y[0] for x, y in zip(iv, iv[1:])) and S["pass"] and S["depth_over_kappa_minus"] == "23/25")
    l0 = math.ceil(1 / H0 - Fraction(1, 2))
    Rr.add("S5.l0", l0 == 20 and Fraction(2, 2 * 20 + 1) <= H0 < Fraction(2, 2 * 19 + 1))
    # S6 sampled cover
    rng = random.Random(20261027)
    kmf, kpf, OmPf = float(km.mid()), float(kp.mid()), float(OmP.mid())
    af = {l: float(Fraction(afix[l]["a_prime"])) for l in range(20)}
    zaf = {l: float(Fraction(afix[l]["zeta_A_over_OmegaP"])) for l in range(20)}
    bpf, b0f = float(bp_ok), float(b0_ok)
    uncovered = []
    counts = {}
    N = 40000
    for i in range(N):
        if i % 4 == 0:
            l = rng.randint(0, 25)
        else:
            l = int(math.exp(rng.uniform(0, math.log(1e6))))
        y = -rng.uniform(0, float(E)) * kmf                 # Im omega in [-e kappa_-, 0]
        if i % 10 == 0:
            y = -float(E) * kmf                             # the edge itself
        hl = 1 / (l + 0.5)
        scale = max(1.0, (l + 0.5) * OmPf / kmf)
        mode = rng.random()
        if mode < 0.3:
            x = rng.uniform(0, 10) * kmf
        elif mode < 0.9:
            x = rng.uniform(0, 3 * scale) * kmf
        else:                                                # a seam point
            cands = [8 * kmf, 6 * kmf, 0.9 * OmPf / hl, 0.8 * OmPf / hl, 1.3 * OmPf / hl, 1.2 * OmPf / hl]
            if l <= 19:
                cands += [af[l]]
            x = rng.choice(cands)
        if rng.random() < 0.5:
            x = -x
        xr = abs(x)                                          # symmetry
        reg = []
        if y >= -float(DD) * kmf:
            reg.append("S_delta (P2-4)")
        if abs(complex(xr, y) - complex(0, -kpf)) < 1e-15:
            reg.append("omega0 (Lemma E)")
        if l <= 19:
            if xr <= af[l] and -bpf <= y <= -b0f:
                reg.append("band box")
            if xr >= zaf[l] * OmPf / hl and y >= -float(E) * kmf:
                reg.append("A at fixed h")
        else:
            if xr <= 8 * kmf and -0.95 * kmf <= y <= -0.55 * kmf:
                reg.append("L")
            if 6 * kmf <= xr <= 0.9 * OmPf / hl:
                reg.append("S")
            if 0.8 * OmPf / hl <= xr <= 1.3 * OmPf / hl:
                reg.append("T")
            if xr >= 1.2 * OmPf / hl:
                reg.append("A")
        for r_ in reg:
            counts[r_] = counts.get(r_, 0) + 1
        if not reg:
            uncovered.append({"l": l, "re_over_km": x / kmf, "im_over_km": y / kmf})
    Rr.add("S6.sampled_cover", not uncovered, note=f"{N} points; uncovered {len(uncovered)}")
    out["sampled_cover"] = {"points": N, "region_counts": counts, "uncovered": uncovered[:20]}
    out["checks"] = Rr.items
    out["failed"] = [k for k, v in Rr.items.items() if not v["ok"]]
    out["n_checks"] = len(Rr.items)
    out["pass"] = not out["failed"]
    (HERE / "indep_assembly_out.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: out[k] for k in ("n_checks", "failed", "pass", "b_prime", "b0", "b_prime_minus_e_kappa",
                                          "beta", "sampled_cover")}, indent=1))


if __name__ == "__main__":
    main()
