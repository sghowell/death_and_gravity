"""S6 review, check C8: float argument-principle counts with the reviewer's own W_l (no p2 code).

W_l(w) = [Delta (Z_+ Z_c' - Z_+' Z_c) + 2 i w r^2 Z_+ Z_c]/(r_+ r_c), with Z_+ (a_0 = 1) and Z_c (c_0 = 1) from
Frobenius series at r_+ and r_c and Taylor continuation to r_m (numpy complex128, vectorised over w).
The parameters come from indep_rnds.Point (Arb, from the definitions), rounded to float.

Not rigorous: arg W is tracked along the boundary with adaptive bisection until every step turns by
less than 0.3 rad.  The certified counts are the author's; this is an independent cross-check of
the windings (0 at l = 0, -1 for 1 <= l <= 19), of the pole at omega0, of omega_*, of the symmetry
halving, and of the l = 19/20 seam (l = 20, 21 on boxes reaching into part (i)'s regime A).

Run:  python indep_float.py      -> indep_float_out.json
"""
from __future__ import annotations

import json
import math
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
from flint import acb, ctx

import indep_rnds as R

HERE = Path(__file__).resolve().parent
CERT = HERE.parent


class FloatW:
    def __init__(self, NF=250, NT=30):
        ctx.prec = 128
        pt = R.Point()
        self.ap = pt
        f = lambda x: float(x.mid())
        self.Dc = [f(c) for c in pt.Dc]
        self.rm, self.rp, self.rc = f(pt.rm), f(pt.rp), f(pt.rc)
        self.km, self.kp, self.kc = f(pt.km), f(pt.kp), f(pt.kc)
        self.OmegaP = f(pt.OmegaP)
        self.NF, self.NT = NF, NT
        self.xm = self.rp + 0.5 * (self.rc - self.rp)
        self.t0 = 0.4 * (self.rp - self.rm)            # Frobenius at r_+ (radius r_+ - r_-)
        self.s0 = 0.1 * (self.rc - self.rp)            # Frobenius at r_c (radius r_c - r_+)

    def _dist(self, x):
        return min(abs(x - self.rp), abs(x - self.rc), abs(x - self.rm))

    def _grid(self, x0, x1, sgn, wmax, L):
        """steps with h <= min(d/4, 1/(4 k)), k = (|Delta'| + 2 wmax r^2)/|Delta| + sqrt((2 wmax r + L)/|Delta|)"""
        pts = [x0]
        x = x0
        while (x1 - x) * sgn > 1e-15:
            D = abs(sum(self.Dc[k] * x ** k for k in range(5)))
            Dp = abs(sum(k * self.Dc[k] * x ** (k - 1) for k in range(1, 5)))
            kx = (Dp + 2 * wmax * x * x) / D + math.sqrt((2 * wmax * x + L) / D)
            h = min(0.25 * self._dist(x), 0.25 / kx, 0.05)
            nx = x + sgn * h
            if (nx - x1) * sgn > 0:
                nx = x1
            pts.append(nx)
            x = nx
        return pts

    @staticmethod
    def _shift(p, x0):
        n = len(p)
        return [sum(p[k] * math.comb(k, j) * x0 ** (k - j) for k in range(j, n)) for j in range(n)]

    def _polys(self, w, L, side):
        """Delta, B, C coefficient lists; B, C entries are arrays over w"""
        D = self.Dc
        Dp = [D[k] * k for k in range(1, len(D))]
        B = [Dp[0] + 0 * w, Dp[1] + 0 * w, Dp[2] - 2j * side * w, Dp[3] + 0 * w]
        C = [-L + 0 * w, -2j * side * w]
        return D, B, C

    def _frob(self, w, L, side, rh, t):
        D, B, C = self._polys(w, L, side)
        # shift to rh: Delta (constants), B and C (arrays)
        Ds = self._shift(D, rh)
        Ds[0] = 0.0
        Bs = [sum(B[k] * math.comb(k, j) * rh ** (k - j) for k in range(j, 4)) for j in range(4)]
        Cs = [C[0] + C[1] * rh, C[1]]
        a = [np.ones_like(w)]
        for m in range(1, self.NF + 1):
            acc = 0
            for k in range(2, 5):
                j = m - k + 1
                if 0 <= j < len(a):
                    acc = acc + Ds[k] * ((m - k + 1) * (m - k)) * a[j]
            for k in range(1, 4):
                j = m - k
                if 0 <= j < len(a):
                    acc = acc + Bs[k] * (m - k) * a[j]
            for k in range(2):
                j = m - 1 - k
                if 0 <= j < len(a):
                    acc = acc + Cs[k] * a[j]
            a.append(-acc / (m * (Ds[1] * (m - 1) + Bs[0])))
        Z = np.zeros_like(w)
        dZ = np.zeros_like(w)
        for n in range(self.NF, -1, -1):
            Z = Z * t + a[n]
        for n in range(self.NF, 0, -1):
            dZ = dZ * t + n * a[n]
        return Z, dZ

    def _step(self, w, L, side, x0, Z0, Z1, h):
        D, B, C = self._polys(w, L, side)
        d = self._shift(D, x0)
        b = [sum(B[k] * math.comb(k, j) * x0 ** (k - j) for k in range(j, 4)) for j in range(4)]
        c = [C[0] + C[1] * x0, C[1]]
        z = [Z0, Z1]
        for m in range(2, self.NT + 1):
            acc = 0
            for i in range(1, 5):
                j = m - i
                if j < 0:
                    continue
                t = 0
                if i < len(d):
                    t = t + d[i] * ((m - i) * (m - i - 1))
                if 0 <= i - 1 < len(b):
                    t = t + b[i - 1] * (m - i)
                if 0 <= i - 2 < len(c):
                    t = t + c[i - 2]
                acc = acc + t * z[j]
            z.append(-acc / (d[0] * (m * (m - 1))))
        Z = np.zeros_like(w)
        dZ = np.zeros_like(w)
        for n in range(self.NT, -1, -1):
            Z = Z * h + z[n]
        for n in range(self.NT, 0, -1):
            dZ = dZ * h + n * z[n]
        return Z, dZ

    def W(self, w, l):
        w = np.asarray(w, dtype=complex)
        L = l * (l + 1)
        wmax = float(np.max(np.abs(w)))
        gp = self._grid(self.rp + self.t0, self.xm, +1, wmax, L)
        gc = self._grid(self.rc - self.s0, self.xm, -1, wmax, L)
        Zp, dZp = self._frob(w, L, +1, self.rp, gp[0] - self.rp)
        for x0, x1 in zip(gp[:-1], gp[1:]):
            Zp, dZp = self._step(w, L, +1, x0, Zp, dZp, x1 - x0)
        Zc, dZc = self._frob(w, L, -1, self.rc, gc[0] - self.rc)
        for x0, x1 in zip(gc[:-1], gc[1:]):
            Zc, dZc = self._step(w, L, -1, x0, Zc, dZc, x1 - x0)
        r = self.xm
        Dl = sum(self.Dc[k] * r ** k for k in range(5))
        return (Dl * (Zp * dZc - dZp * Zc) + 2j * w * r * r * Zp * Zc) / (self.rp * self.rc)


def winding(FW, l, corners, n0=400, tol=0.3, maxlev=14):
    """counterclockwise polygon corners (complex, absolute units); returns (winding, min|W|, points, segments)"""
    pts = []
    for z0, z1 in zip(corners, corners[1:] + corners[:1]):
        n = max(8, int(n0 * abs(z1 - z0) / max(abs(b - a) for a, b in zip(corners, corners[1:] + corners[:1]))))
        pts += list(z0 + (z1 - z0) * np.arange(n) / n)
    pts.append(corners[0])
    z = np.array(pts)
    W = FW.W(z, l)
    for lev in range(maxlev):
        d = np.angle(W[1:] / W[:-1])
        bad = np.nonzero(np.abs(d) > tol)[0]
        if len(bad) == 0:
            break
        mids = (z[bad] + z[bad + 1]) / 2
        Wm = FW.W(mids, l)
        z = np.insert(z, bad + 1, mids)
        W = np.insert(W, bad + 1, Wm)
    d = np.angle(W[1:] / W[:-1])
    return float(np.sum(d) / (2 * np.pi)), float(np.min(np.abs(W))), len(z), float(np.max(np.abs(d)))


def main():
    t0, c0 = time.time(), time.process_time()
    FW = FloatW()
    km = FW.km
    out = {"tool": "review_s6/indep_float.py (numpy complex128; parameters from indep_rnds.Point)", "boxes": []}
    # accuracy against the Arb evaluator at a few points
    ctx.prec = 1500
    pt_hi = R.Point()
    acc = []
    for l, w in ((0, complex(0.3, -0.7) * km), (5, complex(20.0, -0.9) * km), (19, complex(100.0, -0.6) * km),
                 (19, complex(105.5, -0.92) * km), (1, complex(0.02, -0.80) * km), (20, complex(131, -0.9) * km),
                 (20, complex(60, -0.92) * km), (0, complex(12.9, -0.92) * km), (19, complex(0.001, -0.5) * km)):
        Wa, _ = R.solve_W(pt_hi, acb(w.real, w.imag), l, "Z+", N_frob=700, N_tay=300, t_frac=0.1)
        Wf = FW.W(np.array([w]), l)[0]
        wa = complex(float(Wa.real.mid()), float(Wa.imag.mid()))
        acc.append({"l": l, "w_over_km": [w.real / km, w.imag / km], "float": [Wf.real, Wf.imag],
                    "arb_1500bit_mid": [wa.real, wa.imag], "arb_rel_radius": float(abs(Wa).rad()) / abs(wa),
                    "rel_err": abs(Wf - wa) / abs(wa)})
    out["float_vs_arb"] = acc
    A_fixed = {r["l"]: r for r in json.loads((CERT / "beta" / "part2" / "A_fixed_PB.json").read_text())["rows"]}
    b0d = json.loads((CERT / "beta" / "part3" / "band_PB_l0.json").read_text())["region"]
    bp, b0 = float(Fraction(b0d["b_prime"])), float(Fraction(b0d["b0"]))
    for l in range(20):
        a = float(Fraction(A_fixed[l]["a_prime"]))
        corners = [complex(-a, -bp), complex(a, -bp), complex(a, -b0), complex(-a, -b0)]
        wn, mn, npts, dmax = winding(FW, l, corners)
        # the right half only (H), doubled: the author's symmetry halving
        H = [complex(0, -bp), complex(a, -bp), complex(a, -b0), complex(0, -b0)]
        z = []
        for z0, z1 in zip(H[:-1], H[1:]):
            n = 300
            z += list(z0 + (z1 - z0) * np.arange(n) / n)
        z.append(H[-1])
        z = np.array(z)
        Wv = FW.W(z, l)
        for _ in range(14):
            d = np.angle(Wv[1:] / Wv[:-1])
            bad = np.nonzero(np.abs(d) > 0.3)[0]
            if len(bad) == 0:
                break
            mids = (z[bad] + z[bad + 1]) / 2
            z = np.insert(z, bad + 1, mids)
            Wv = np.insert(Wv, bad + 1, FW.W(mids, l))
        right = float(np.sum(np.angle(Wv[1:] / Wv[:-1])) / np.pi)
        rec = {"l": l, "box_over_km": [a / km, bp / km, b0 / km], "winding_full": wn, "right_half_over_pi": right,
               "min_abs_W_full": mn, "points": npts, "expected": 0 if l == 0 else -1}
        rec["as_expected"] = abs(wn - rec["expected"]) < 0.05 and abs(right - rec["expected"]) < 0.05
        out["boxes"].append(rec)
        print("band", l, round(wn, 4), round(right, 4), mn, npts, flush=True)
    # KAT-like small boxes and the pole
    kat = []
    for l, bx, expect in ((0, (-0.03, 0.03, -0.86, -0.79), -1), (1, (-0.03, 0.03, -0.86, -0.79), -1),
                          (7, (-0.03, 0.03, -0.86, -0.79), -1), (19, (-0.03, 0.03, -0.86, -0.79), -1),
                          (0, (-0.03, 0.03, -0.93, -0.89), 1), (0, (-0.03, 0.03, -0.93, -0.79), 0)):
        x0, x1, y0, y1 = (v * km for v in bx)
        wn, mn, npts, _ = winding(FW, l, [complex(x0, y0), complex(x1, y0), complex(x1, y1), complex(x0, y1)], n0=200)
        kat.append({"l": l, "box_over_km": bx, "winding": wn, "expected": expect, "ok": abs(wn - expect) < 0.05})
        print("kat", l, bx, round(wn, 4), flush=True)
    out["small_boxes"] = kat
    # omega_*: Newton on W_0 along the imaginary axis (W_0 real there)
    y = 0.9094 * km
    for _ in range(30):
        f = FW.W(np.array([complex(0, -y)]), 0)[0]
        h = 1e-7 * km
        fp = (FW.W(np.array([complex(0, -y - h)]), 0)[0] - FW.W(np.array([complex(0, -y + h)]), 0)[0]) / (2 * h)
        y = y - (f / fp).real
    out["omega_star_float_over_km"] = -y / km
    # the seam: l = 20, 21 on [-1.3 (l + 1/2) Omega_P, 1.3 (l + 1/2) Omega_P] x [-b', -b0] (into regime A)
    seam = []
    for l in (20, 21):
        a = 1.3 * (l + 0.5) * FW.OmegaP
        wn, mn, npts, _ = winding(FW, l, [complex(-a, -bp), complex(a, -bp), complex(a, -b0), complex(-a, -b0)],
                                  n0=1500)
        seam.append({"l": l, "a_over_km": a / km, "winding": wn, "expected": -1, "min_abs_W": mn, "points": npts,
                     "ok": abs(wn + 1) < 0.05})
        print("seam", l, round(wn, 4), mn, npts, flush=True)
    out["seam_l20_21"] = seam
    out["float_accuracy_ok"] = all(a["rel_err"] < 1e-10 for a in acc)
    out["pass"] = bool(out["float_accuracy_ok"] and all(b["as_expected"] for b in out["boxes"]) and all(k["ok"] for k in kat)
                       and all(s["ok"] for s in seam) and abs(out["omega_star_float_over_km"] + 0.909380978994805) < 1e-9)
    out["timing_s"] = {"wall": round(time.time() - t0, 1), "cpu": round(time.process_time() - c0, 1)}
    (HERE / "indep_float_out.json").write_text(json.dumps(out, indent=1) + "\n")
    print("pass", out["pass"], out["timing_s"])


if __name__ == "__main__":
    main()
