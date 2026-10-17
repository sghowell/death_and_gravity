"""S3 review: adversarial checks of the zsol box (omega-disc) enclosure, the part of the pipeline
that no pre-registered control exercises directly (C0 compares point values only).

A. M versus the true sup (all 2351 production discs).  For every stored segment, the reviewer's
   independent float W (indep_w.py) is sampled on a 9 x 9 grid of the box [c +- rho]^2 and on 32
   points of the circle |w - c| = rho; the maximum must not exceed the stored M_upper (M is claimed
   to bound sup |W| on the box).  Also reported: how loose M is (M / max sampled |W|).
B. Containment in tight boxes.  zsol.wronskian_z(pt, l, c, box) on boxes of half-width
   rho in {kappa_-/16, 1e-4 kappa_-, 1e-8 kappa_-} at high-frequency corners and interior points of
   R'_l (l = 1, 5, 10, 19); the box ball must contain
   (i) the independent float W at the 4 corners, 4 edge midpoints and 8 random interior points
       (float error <= 3e-10 relative at the C0 points; 1e-9 relative is added as a radius), and
   (ii) S1's arbw enclosure (rigorous, independent code) at two corners (precision raised until
       the relative radius is below 1e-15).
   Reported: the relative width of the box ball (tightness) and every containment verdict.
Usage: python box_checks.py A | B      (writes box_checks_A.json / box_checks_B.json)
"""
from __future__ import annotations

import json
import random
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from indep_w import Geo, W_l  # noqa: E402

CERT = HERE.parent / "part3"


def part_A():
    g = Geo()
    out = {"per_l": {}}
    worst_ratio = 0.0
    nviol = 0
    for l in range(20):
        st = json.loads((CERT / f"winding_PB_l{l}.json").read_text())
        pts, owners = [], []
        for k, r in enumerate(st["pieces"]):
            c = complex(float(Fr(r["c"][0])), float(Fr(r["c"][1])))
            rho = float(Fr(r["rho"]))
            u = np.linspace(-1, 1, 9)
            grid = (c + rho * (u[:, None] + 1j * u[None, :])).ravel()
            circ = c + rho * np.exp(2j * np.pi * np.arange(32) / 32)
            p = np.concatenate([grid, circ])
            pts.append(p)
            owners += [k] * len(p)
        allp = np.concatenate(pts)
        owners = np.array(owners)
        order = np.argsort(np.abs(allp))
        Wv = np.empty(len(allp), dtype=complex)
        for i in range(0, len(allp), 400):
            idx = order[i:i + 400]
            Wv[idx] = W_l(g, l, allp[idx], rtol=1e-12)
        absW = np.abs(Wv)
        per = []
        for k, r in enumerate(st["pieces"]):
            mx = float(np.max(absW[owners == k]))
            M = r["M_upper"]
            viol = mx > M
            nviol += viol
            per.append({"key": r["key"], "max_sampled_absW": mx, "M_upper": M, "M_over_max": M / mx,
                        "violation": bool(viol)})
            worst_ratio = max(worst_ratio, mx / M)
        out["per_l"][str(l)] = {"segments": len(per), "violations": sum(p["violation"] for p in per),
                                "min_M_over_max": min(p["M_over_max"] for p in per),
                                "median_M_over_max": float(np.median([p["M_over_max"] for p in per])),
                                "max_M_over_max": max(p["M_over_max"] for p in per)}
        print(l, out["per_l"][str(l)], flush=True)
    out["violations_total"] = int(nviol)
    out["max_sampled_over_M"] = worst_ratio
    (HERE / "box_checks_A.json").write_text(json.dumps(out, indent=1) + "\n")


def part_B():
    sys.path.insert(0, str(HERE.parents[1] / "src"))
    from flint import acb, arb, fmpq
    from p2 import zsol
    from p2.arbutil import precision
    from p2.arbw import wronskian as arbw_w
    from p2.params import Point
    import p2.arbw as arbw_mod

    g = Geo()
    pt = Point("PB", 256)
    km = float(pt.km.mid())
    rng = random.Random(20261017)

    def q(x):
        return fmpq(round(x * 2 ** 40), 2 ** 40)

    # arbw's Frobenius start needs more terms at large |w| (the C0 issue); raise min_terms only
    orig = arbw_mod.frobenius

    def frob(ptx, which, w, L, t, min_terms=20):
        kap = ptx.kp if which == "+" else ptx.kc
        ra = float((abs(acb(w)) / (2 * kap)).upper())
        return orig(ptx, which, w, L, t, min_terms=max(min_terms, int(2 * ra) + 8))

    def arbw_ref(l, wq):
        arbw_mod.frobenius = frob
        try:
            for prec in (256, 512, 1024, 2048):
                ptp = Point("PB", prec)
                with precision(prec):
                    Wa = arbw_w(ptp, l, acb(arb(wq[0]), arb(wq[1])))
                    if float(abs(Wa).rad() / abs(Wa).mid()) < 1e-15:
                        return Wa, prec
            return Wa, prec
        finally:
            arbw_mod.frobenius = orig

    cases = []
    for l in (1, 5, 10, 19):
        a = float(Fr(json.loads((CERT / f"winding_PB_l{l}.json").read_text())["region"]["a"])) / km
        for (x, y) in ((a, -0.55), (a, 0.1), (0.5 * a, -0.55), (0.0, 0.1)):
            for rk in (1 / 16, 1e-4, 1e-8):
                cases.append((l, x, y, rk))
    res = []
    for l, x, y, rk in cases:
        t0 = time.time()
        cq = (q(x * km), q(y * km))
        rq = q(rk * km) if rk >= 1e-4 else fmpq(round(rk * km * 2 ** 70), 2 ** 70)
        with precision(128):
            c = acb(arb(cq[0]), arb(cq[1]))
            rho = arb(rq)
            box = acb(c.real + arb(0, 1) * rho, c.imag + arb(0, 1) * rho)
            _, Wb = zsol.wronskian_z(pt, l, c, box, theta=0.25)
            Wc = zsol.wronskian_z(pt, l, c, theta=3.0)[0]
        cf, rf = complex(float(cq[0]), float(cq[1])), float(rq)
        sample = [cf + rf * complex(sx, sy) for sx in (-1, 1) for sy in (-1, 1)]
        sample += [cf + rf * z for z in (1, -1, 1j, -1j)]
        sample += [cf + rf * complex(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(8)]
        Wf = W_l(g, l, np.array(sample), rtol=1e-13)
        float_ok = []
        for wv in Wf:
            tol = 1e-9 * abs(wv)            # float W error is <= 3e-10 relative (C0 points)
            float_ok.append(bool(Wb.real.lower() - tol <= wv.real <= Wb.real.upper() + tol
                                 and Wb.imag.lower() - tol <= wv.imag <= Wb.imag.upper() + tol))
        arbw_ok = []
        for sx, sy in ((-1, -1), (1, 1)):
            wq = (cq[0] + sx * rq, cq[1] + sy * rq)
            Wa, prec = arbw_ref(l, wq)
            with precision(256):
                arbw_ok.append({"contained": bool(Wb.contains(Wa)), "overlap": bool(Wb.overlaps(Wa)),
                                "arbw_prec": prec})
        relw = float(max(Wb.real.rad(), Wb.imag.rad())) / float(abs(Wc).mid())
        rec = {"l": l, "c_over_km": [x, y], "rho_over_km": rk, "box_rel_halfwidth": relw,
               "float_points_contained": sum(float_ok), "float_points": len(float_ok),
               "arbw_corners": arbw_ok, "wall_s": round(time.time() - t0, 1)}
        res.append(rec)
        print(json.dumps(rec), flush=True)
    out = {"cases": res,
           "all_float_contained": all(r["float_points_contained"] == r["float_points"] for r in res),
           "all_arbw_contained": all(x["contained"] for r in res for x in r["arbw_corners"])}
    (HERE / "box_checks_B.json").write_text(json.dumps(out, indent=1) + "\n")
    print({k: out[k] for k in ("all_float_contained", "all_arbw_contained")})


if __name__ == "__main__":
    {"A": part_A, "B": part_B}[sys.argv[1]]()
