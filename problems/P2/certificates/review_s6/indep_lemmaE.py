"""S6 review, check B7: Lemma E at P_B with the reviewer's own Arb code (indep_rnds.py; nothing from p2).

1. P_B from its definitions; compared with certificates/params.json (read as data).
2. Wt_l(omega0), omega0 = -i kappa_+, for l = 0, 1, 2, 3, 5, 10, 19 and l = 20, 50, 100: the exponent-1
   Frobenius solution Z2 = t + ... at r_+ (B_0 = 0 exactly), validated Taylor continuation, Z_c, and
   Wt_l = [Delta (Z2 Zc' - Z2' Zc) + 2 i omega0 r^2 Z2 Zc]/(r_+ r_c).  Compared with epole_PB.json.
3. The pole: eps W_l(omega0 + eps) against Res = i (2 kappa_+ r_+ + L)/(2 r_+^2) Wt_l for
   eps = 10^-k kappa_- e^{i pi/4}, k = 3, 5, 7 (l = 0, 1, 5, 19), rigorous balls for both sides.
4. omega_*: a certified sign change of the real function W_0(-i y kappa_-) across P2-K6's ball.

Run (from this directory):  python indep_lemmaE.py      -> indep_lemmaE_out.json
"""
from __future__ import annotations

import json
import time
from fractions import Fraction
from pathlib import Path

from flint import acb, arb, ctx

import indep_rnds as R

HERE = Path(__file__).resolve().parent
CERT = HERE.parent


def main():
    ctx.prec = 200
    t0, c0 = time.time(), time.process_time()
    out = {"tool": "review_s6/indep_rnds.py (python-flint Arb, 200 bits); no p2 code"}
    pt = R.Point()
    out["point"] = pt.summary()
    # 1. against params.json
    par = json.loads((CERT / "params.json").read_text())["points"]["PB"]["quantities"]
    cmp = {}
    for mine, key in (("km", "kappa_minus"), ("kp", "kappa_plus"), ("kc", "kappa_c"), ("rP", "r_P"),
                      ("OmegaP", "Omega_P"), ("gamma0", "gamma0")):
        if key in par:
            cmp[key] = bool(getattr(pt, mine).overlaps(arb(par[key]["ball"])))
    out["params_json_overlap"] = cmp
    # 2. Wt_l
    ep = json.loads((CERT / "beta" / "epole_PB.json").read_text())
    stored = {r["l"]: r for r in ep["Wt_l_le_19"] + ep["Wt_samples_l_ge_20"]}
    w0 = acb(0, -pt.kp)
    rows = []
    wt = {}
    for l in (0, 1, 2, 3, 5, 10, 19, 20, 50, 100):
        tl = time.process_time()
        W, info = R.solve_W(pt, w0, l, "Z2")
        wt[l] = W
        st = stored.get(l)
        st_re = arb(st["Wt"].split("] + [")[0] + "]") if st else None
        rows.append({"l": l, "Wt": W.real.str(30, radius=True), "imag_part": W.imag.str(5, radius=True),
                     "nonzero": bool(not W.real.contains(0)), "sign": "+" if bool(W.real > 0) else "-",
                     "author_Wt": st["Wt"].split(" + ")[0] if st else None,
                     "overlaps_author": bool(st_re.overlaps(W.real)) if st else None,
                     "rel_radius": float((W.real.rad() / abs(W.real)).upper()),
                     "info": info, "cpu_s": round(time.process_time() - tl, 2)})
        print("Wt", l, W.real.str(20), rows[-1]["overlaps_author"], flush=True)
    out["Wt"] = rows
    # 3. the pole
    pole = []
    for l in (0, 1, 5, 19):
        L = l * (l + 1)
        res = I_ = acb(0, 1)
        res = I_ * (2 * pt.kp * pt.rp + L) / (2 * pt.rp * pt.rp) * wt[l]
        for k in (3, 5, 7):
            eps = acb(1, 1) / arb(2).sqrt() * arb(10) ** (-k) * pt.km
            W, _ = R.solve_W(pt, w0 + eps, l, "Z+")
            rel = abs(eps * W - res) / abs(res)
            pole.append({"l": l, "eps_over_kappa_minus": f"1e-{k} e^(i pi/4)", "eps_W": (eps * W).str(15),
                         "residue": res.str(15), "rel_diff": rel.str(5), "rel_diff_over_eps":
                         (rel / (arb(10) ** (-k))).str(5)})
            print("pole", l, k, rel.str(5), flush=True)
    out["pole"] = pole
    ratios = {}
    for l in (0, 1, 5, 19):
        r_ = [float(arb(p["rel_diff"]).mid()) for p in pole if p["l"] == l]
        ratios[l] = [r_[0] / r_[1], r_[1] / r_[2]]
    out["pole_ratio_between_eps_steps"] = ratios
    out["pole_simple_O_eps"] = all(99 < x < 101 for v in ratios.values() for x in v)
    # 4. omega_*: W_0(-i y kappa_-) is real on the imaginary axis; sign change across P2-K6's ball
    kz = json.loads((CERT / "part3" / "krawczyk_PB_l0.json").read_text())["enclosure_over_kappa_minus"]
    im = arb(kz["im"])
    y_mid = -im
    sc = {}
    for tag, dy in (("below_ball", -2 * im.rad()), ("above_ball", 2 * im.rad())):
        y = arb(y_mid.mid()) + arb(dy)
        W, _ = R.solve_W(pt, acb(0, -y * pt.km), 0, "Z+")
        sc[tag] = {"y": y.str(35), "W0": W.real.str(10, radius=True), "imag": W.imag.str(3, radius=True),
                   "sign": None if W.real.contains(0) else ("+" if bool(W.real > 0) else "-")}
    sc["sign_change"] = sc["below_ball"]["sign"] is not None and sc["above_ball"]["sign"] is not None and \
        sc["below_ball"]["sign"] != sc["above_ball"]["sign"]
    sc["P2_K6_ball_im"] = kz["im"]
    out["omega_star_sign_change"] = sc
    out["pass"] = bool(all(r["nonzero"] for r in rows) and all(r["overlaps_author"] for r in rows if
                                                                 r["overlaps_author"] is not None)
                       and all(cmp.values()) and out["pole_simple_O_eps"] and sc["sign_change"]
                       and all((r["sign"] == "+") == (r["l"] == 0) for r in rows))
    out["timing_s"] = {"wall": round(time.time() - t0, 1), "cpu": round(time.process_time() - c0, 1)}
    (HERE / "indep_lemmaE_out.json").write_text(json.dumps(out, indent=1) + "\n")
    print("pass", out["pass"], out["timing_s"])


if __name__ == "__main__":
    main()
