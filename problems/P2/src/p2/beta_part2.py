"""S6 part (ii'): regime A at fixed h_l = 1/(l + 1/2) for l <= 19, at depth e kappa_- (pre-registration O3).

S4's regime A (part1_T, design DESIGN_A: phi0 = -1/5, Gamma = 0) on S2's rotated contour, with the
uniform bound F = E0/h0 + E1 + h0 min(0, E2) at h0 := h_l, i.e. valid for every h <= h_l.  On the
zeta-box {zeta_A} x [-e kappa_-, 0] it gives E >= mu/h > 0; above zeta_A it is monotone in Re omega
(S4 Lemma W7, phi <= phi0 < -arctan(e kappa_- h_l/zeta_A)).  So W_l(omega) != 0 for
Re omega >= zeta_A/h_l, -e kappa_- <= Im omega <= 0 (Re omega + Im omega > 0 there: S4 W2).

zeta_A(l) (the rule of section 3): the smallest multiple of 1/20 in [1.1, 12] Omega_P with float margin
>= 0.02 a* and W7, certified in Arb at h0 = h_l; on an Arb failure zeta_A += 1/20, at most 4 times.
a'_l := zeta_A(l) Omega_P / h_l rounded up to 2^-30 (the half-width of the band box B_l).

Run:  python -m p2.beta_part2            (writes certificates/beta/part2/A_fixed_PB.json)
"""
from __future__ import annotations

import math
import sys
from fractions import Fraction

from flint import arb

from .arbutil import Timer, precision
from .beta_common import (CERT, EDGE, L_SMALL, PREC, ZA_FLOAT_MARGIN, ZA_GRID, ZA_RETRIES, FloatTA, commit_tag,
                          deep_point, dy_up, write_json)
from .m3 import fa

OUT = CERT / "part2" / "A_fixed_PB.json"


def w7(pt, h: Fraction, zA: Fraction):
    """Lemma W7 at fixed h: phi0 = -1/5 < -arctan(e kappa_- h / (zeta_A Omega_P)), and phi_inf > -pi."""
    from . import part1_T as P
    bound = -((pt.beta * fa(h)) / (fa(zA) * pt.OmegaP)).atan()
    phi0, phinf = fa(Fraction(P.DESIGN_A["phi0"])), fa(Fraction(P.DESIGN_A["phi_inf"]))
    return {"minus_arctan_bound": bound.str(10), "pass": bool(phinf < phi0 < bound and phinf > -arb.pi())}


def float_zeta_A(f: FloatTA, pt, l: int):
    from . import part1_T as P
    h = Fraction(2, 2 * l + 1)
    for z in ZA_GRID:
        v = f.box(P.DESIGN_A, (z, z), EDGE, h)
        ok7 = -0.2 < -math.atan(float(EDGE) * f.km * float(h) / (float(z) * f.OmP))
        if v["F_min_over_astar"] >= ZA_FLOAT_MARGIN and v["E0_min"] >= 0 and ok7:
            return z, v["F_min_over_astar"]
    return None, None


def certificate():
    from . import part1_T as P
    pt = deep_point()
    f = FloatTA(pt)
    rows = []
    with precision(PREC), Timer() as tm:
        _, rays, tails = P.setup(pt)
        ch = P.TChecker(pt, rays, tails, P.DESIGN_A)
        for l in L_SMALL:
            h = Fraction(2, 2 * l + 1)
            z0, fm = float_zeta_A(f, pt, l)
            if z0 is None:
                rows.append({"l": l, "pass": False, "error": "no float zeta_A in [1.1, 12]"})
                continue
            tries = []
            z = z0
            rec = None
            for k in range(ZA_RETRIES + 1):
                r = P.certify_box(ch, pt, (z, z), P.DESIGN_A, h0=h)
                mono = w7(pt, h, z)
                tries.append({"zeta_A": str(z), "pass": bool(r["pass"] and mono["pass"]), "mu_over_a_star": r["mu_over_a_star"]})
                if r["pass"] and mono["pass"]:
                    rec = (z, r, mono)
                    break
                z = z + Fraction(1, 20)
            if rec is None:
                rows.append({"l": l, "h": str(h), "pass": False, "tries": tries})
                print("A fixed h", l, "FAILED", tries, flush=True)
                continue
            z, r, mono = rec
            aprime = dy_up(fa(z) * pt.OmegaP / fa(h))
            rows.append({"l": l, "h": str(h), "zeta_A_over_OmegaP": str(z), "zeta_A_float": str(z0), "float_margin": fm,
                         "tries": tries, "a_prime": str(aprime), "a_prime_over_kappa_minus": (arb(aprime) / pt.km).str(10),
                         "a_prime_ge_zeta_A_OmegaP_over_h": bool(arb(aprime) >= fa(z) * pt.OmegaP / fa(h)),
                         "W7": mono, "box": {k: r[k] for k in ("zeta_over_OmegaP", "nu_low", "h0", "mu_cert", "mu_over_a_star",
                                                               "E0_low_min", "c_min", "tail_left", "tail_right", "cells",
                                                               "worst_cell", "n_fails", "jost_decay", "pass", "leaves")},
                         "pass": True})
            print("A fixed h", l, z, r["mu_over_a_star"], (arb(aprime) / pt.km).str(6), flush=True)
    ok = len(rows) == len(L_SMALL) and all(x["pass"] for x in rows)
    return {"claim": "S6 part (ii') at P_B: for each l <= 19, W_l(omega) != 0 for Re omega >= zeta_A(l)/h_l and "
                     "-e kappa_- <= Im omega <= 0 (e = 23/25), by S4's regime A at h0 = h_l (uniform in h <= h_l) and "
                     "Lemma W7", "point": "PB", "depth_over_kappa_minus": str(EDGE), "design": P.DESIGN_A,
            "rule": {"zeta_grid": "multiples of 1/20 in [1.1, 12] Omega_P", "float_margin_over_a_star": ZA_FLOAT_MARGIN,
                     "retries": ZA_RETRIES},
            "geometry": "S2's contour (m3.geometry_certificate, A0)", "rows": rows, "pass": ok, "prec_bits": PREC,
            "timing_s": {"wall": round(tm.wall, 1), "cpu": round(tm.cpu, 1)}, "commit": commit_tag()}


def main(argv):
    c = certificate()
    write_json(OUT, c)
    print("part (ii') pass:", c["pass"], c["timing_s"], c["commit"])


if __name__ == "__main__":
    main(sys.argv)
