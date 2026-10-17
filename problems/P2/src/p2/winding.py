"""KAT-6: validated winding numbers of W_l on small boxes and Krawczyk enclosures of simple zeros.

Why not evaluate W on whole edge segments: interval evaluation of the Frobenius/Taylor
recurrences over an omega-box overestimates the image by a roughly constant factor (~4e7
measured at P_A, l = 1: the radii of the recurrence propagate through the absolute-value
majorant), so segment enclosures would be useless.  Instead (Cauchy interpolation):
  1. M >= sup |W_l| on the disc |w - w_c| <= rho: one box evaluation (overestimation harmless);
  2. S_k = (1/K) sum_j W_l(w_c + rho' e^{2 pi i j/K}) e^{-2 pi i j k/K}, rho' = x rho, from K
     validated point evaluations; since W_l is analytic on the disc,
        W_k rho'^k = S_k - sum_{p>=1} W_{k+pK} rho'^{k+pK},  |W_n| <= M/rho^n,
     so W_k = S_k/rho'^k +- M x^K/((1 - x^K) rho^k)   (aliasing bound);
  3. W_l(w_c + d) in P(d) + disc(M (|d|/rho)^K/(1 - |d|/rho)),  P(d) = sum_{k<K} W_k d^k.
On each boundary segment the enclosure P(segment ball) + remainder must exclude 0 (convex), so
the arg change along the segment is the principal Arg of the ratio of the endpoint values;
the increments add up to 2 pi * winding (a ball containing exactly one integer).

Krawczyk.  For Z = z0 + [-r, r] + i[-r, r] and Y ~ 1/W'(z0):
    K(Z) = z0 - Y W(z0) + (1 - Y W'(Z)) (Z - z0)  subset of int(Z)  =>  unique zero of W_l in Z
(mean-value form; W' over Z enclosed via a central difference at z0 with the analytic remainder
bound M h^2/(R^3 (1 - h^2/R^2)) plus the Cauchy bound 2 M r/(R - r)^2 for W'(Z) - W'(z0), where
M >= |W_l| on the box z0 + [-R, R] + i[-R, R] is itself a validated enclosure).
"""
from __future__ import annotations

import itertools
import json
import time
from pathlib import Path

import mpmath as mp
from flint import acb, arb, fmpq

from .arbutil import git_commit, precision
from .arbw import wronskian
from .params import Point

CERT_DIR = Path(__file__).resolve().parents[2] / "certificates"


def _q(x, bits=80):
    """Exact dyadic rational near the mpmath/float value x."""
    return fmpq(int(mp.nint(mp.mpf(x) * mp.mpf(2) ** bits)), 2 ** bits)


def _seg_ball(z0, z1):
    """acb ball covering the straight segment [z0, z1] (axis-parallel), endpoints exact acb."""
    mid = (z0 + z1) / 2
    hr = abs((z1 - z0).real) / 2
    hi = abs((z1 - z0).imag) / 2
    return acb(mid.real + arb(0, 1) * hr, mid.imag + arb(0, 1) * hi)


def _excludes_zero(z):
    return not z.real.contains(0) or not z.imag.contains(0)


def cauchy_model(pt: Point, l: int, wc, rho, K=64, x="0.5"):
    """Rigorous Taylor model of W_l around wc: (coefficient balls W_k, M, rho)."""
    x = arb(x)
    box = acb(wc.real + arb(0, 1) * rho, wc.imag + arb(0, 1) * rho)       # contains the disc
    Wbox = wronskian(pt, l, box)
    M = arb(abs(Wbox).upper())
    rp_ = x * rho
    vals = []
    for j in range(K):
        th = 2 * arb.pi() * j / K
        vals.append(wronskian(pt, l, wc + rp_ * acb(th.cos(), th.sin())))
    alias = M * x ** K / (1 - x ** K)
    coef = []
    for k in range(K):
        Sk = acb(0)
        for j in range(K):
            th = -2 * arb.pi() * j * k / K
            Sk += vals[j] * acb(th.cos(), th.sin())
        Sk = Sk / K
        err = alias / rho ** k
        coef.append(Sk / rp_ ** k + acb(arb(0, 1) * err, arb(0, 1) * err))
    return coef, M


def _poly_eval(coef, d):
    acc = acb(0)
    for c in reversed(coef):
        acc = acc * d + c
    return acc


def winding(pt: Point, l: int, corners, wc, rho, K=64, nseg=16, max_split=10):
    """Validated winding number of W_l around the rectangle with exact corners (counterclockwise),
    all corners within |w - wc| < rho/2 (checked)."""
    # W_l is analytic on the disc: it avoids E_+ and E_c (poles at -i n kappa, n >= 1)
    for kap in (pt.kp, pt.kc):
        for n in range(1, 200):
            if not (abs(wc - acb(0, -n * kap)) > rho):
                raise ArithmeticError("interpolation disc meets E_+ or E_c")
    coef, M = cauchy_model(pt, l, wc, rho, K)
    total = arb(0)
    nsegs_used = 0
    cyc = corners + corners[:1]
    worst_ratio = 0.0
    for c in corners:
        if not (abs(c - wc) < rho / 2):
            raise ArithmeticError("box not inside the interpolation disc")

    def encl(z):
        d = z - wc
        s_ = arb(abs(d).upper())
        rem = M * (s_ / rho) ** K / (1 - s_ / rho)
        return _poly_eval(coef, d) + acb(arb(0, 1) * rem, arb(0, 1) * rem)

    for a, b in itertools.pairwise(cyc):
        stack = [(a + (b - a) * k / nseg, a + (b - a) * (k + 1) / nseg, 0) for k in range(nseg)][::-1]
        while stack:
            z0, z1, depth = stack.pop()
            E = encl(_seg_ball(z0, z1))
            if not _excludes_zero(E):
                if depth >= max_split:
                    raise ArithmeticError("segment enclosure contains 0 at max depth")
                zm = (z0 + z1) / 2
                stack.append((zm, z1, depth + 1))
                stack.append((z0, zm, depth + 1))
                continue
            inc = (encl(z1) / encl(z0)).arg()
            if not (abs(inc) < arb.pi()):
                raise ArithmeticError("argument increment not certified")
            total += inc
            nsegs_used += 1
            rad = float(max(E.real.rad(), E.imag.rad()).upper())
            worst_ratio = max(worst_ratio, rad / float(abs(E).mid()))
    wn = total / (2 * arb.pi())
    k = round(float(wn.mid()))
    certified = bool(abs(wn - k) < arb("0.25"))
    return {"winding": k if certified else None, "winding_ball": wn.str(10),
            "certified": certified, "segments": nsegs_used, "point_evaluations": K + 1,
            "taylor_model": {"K": K, "rho_over_box_halfdiag": 4, "M_upper": float(M.upper())},
            "max_enclosure_rel_radius": worst_ratio}


def krawczyk(pt: Point, l: int, z0, r, R, h):
    """Unique-zero certificate in Z = z0 + [-r, r] + i[-r, r]; z0 exact acb point; r, R, h arb."""
    Wz0 = wronskian(pt, l, z0)
    box_R = acb(z0.real + arb(0, 1) * R, z0.imag + arb(0, 1) * R)
    WR = wronskian(pt, l, box_R)
    M = arb(abs(WR).upper())
    Wp = wronskian(pt, l, z0 + h)
    Wm = wronskian(pt, l, z0 - h)
    fd = (Wp - Wm) / (2 * h)
    e_fd = M * h * h / (R ** 3 * (1 - (h / R) ** 2))
    # Z is the square z0 + [-r, r] + i[-r, r]: its corners lie at distance sqrt(2) r from z0,
    # so the Cauchy bound on |W'(z) - W'(z0)| must use s = sqrt(2) r (and needs R > s).
    s_corner = arb(2).sqrt() * r
    if not bool(R > s_corner):
        raise ArithmeticError("Krawczyk: need R > sqrt(2) r")
    e_z = 2 * M * s_corner / (R - s_corner) ** 2
    dW = fd + acb(arb(0, 1), arb(0, 1)) * (e_fd + e_z)
    Y = 1 / acb(fd.real.mid(), fd.imag.mid())
    Zc = acb(arb(0, 1) * r, arb(0, 1) * r)
    K = z0 - Y * Wz0 + (1 - Y * dW) * Zc
    inside = bool(abs(K.real - z0.real) < r and abs(K.imag - z0.imag) < r)
    # tight enclosure: K(Z) itself (contains the unique zero)
    return {"certified_unique_zero": inside, "K": K, "M": M, "W_z0": Wz0, "fd": fd,
            "contraction": float(abs(1 - Y * dW).upper())}


def kat6_case(name, l, z_over_km, label, printed=None, hw_over_km="1e-3", prec=256):
    """Winding on the box z +- hw (kappa_- units) and Krawczyk at the float value z."""
    t0, c0 = time.time(), time.process_time()
    pt = Point(name, prec)
    mp.mp.dps = 80
    z_over_km = mp.mpc(*z_over_km)
    with precision(prec):
        km = float(pt.km.mid())
        zc = z_over_km * mp.mpf(pt.km.mid().str(70, radius=False))
        hw = mp.mpf(hw_over_km) * km
        xr = [_q(mp.re(zc) - hw), _q(mp.re(zc) + hw)]
        yi = [_q(mp.im(zc) - hw), _q(mp.im(zc) + hw)]
        corners = [acb(arb(xr[0]), arb(yi[0])), acb(arb(xr[1]), arb(yi[0])),
                   acb(arb(xr[1]), arb(yi[1])), acb(arb(xr[0]), arb(yi[1]))]
        wc = acb(arb(_q(mp.re(zc), 200)), arb(_q(mp.im(zc), 200)))
        rho = arb(_q(4 * mp.sqrt(2) * hw, 200))
        wres = winding(pt, l, corners, wc, rho)
        # Krawczyk centre: the float zero refined by two Newton steps on midpoints of validated
        # values (any centre is admissible; the certificate is the inclusion K(Z) in int(Z))
        absz = abs(zc)
        hfd = arb(_q(absz * mp.mpf("1e-30"), 300))
        z0 = acb(arb(_q(mp.re(zc), 240)), arb(_q(mp.im(zc), 240)))
        for _ in range(2):
            Wz = wronskian(pt, l, z0)
            d = (wronskian(pt, l, z0 + hfd) - wronskian(pt, l, z0 - hfd)) / (2 * hfd)
            step = acb(Wz.real.mid(), Wz.imag.mid()) / acb(d.real.mid(), d.imag.mid())
            z0 = acb(arb(_q(mp.mpf(z0.real.mid().str(80, radius=False)) - mp.mpf(step.real.mid().str(80, radius=False)), 240)),
                     arb(_q(mp.mpf(z0.imag.mid().str(80, radius=False)) - mp.mpf(step.imag.mid().str(80, radius=False)), 240)))
        r = arb(_q(absz * mp.mpf("1e-20"), 300))
        R = arb(_q(absz * mp.mpf("1e-4"), 200))
        kr = krawczyk(pt, l, z0, r, R, hfd)
        K = kr["K"] / acb(pt.km)
        zf = acb(arb(mp.nstr(mp.re(zc), 70)), arb(mp.nstr(mp.im(zc), 70)))
        float_in_Z = bool(abs(zf.real - z0.real) < r and abs(zf.imag - z0.imag) < r)
        in_box = (arb(xr[0]) < z0.real < arb(xr[1])) and (arb(yi[0]) < z0.imag < arb(yi[1]))
        out = {"point": name, "l": l, "mode": label,
               "box_kappa_minus_units": {"centre": [mp.nstr(mp.re(mp.mpc(z_over_km)), 25), mp.nstr(mp.im(mp.mpc(z_over_km)), 25)],
                                         "half_width": hw_over_km},
               "box_corners_exact_absolute": {"re": [str(xr[0]), str(xr[1])], "im": [str(yi[0]), str(yi[1])]},
               "winding": wres,
               "krawczyk": {"certified_unique_zero": kr["certified_unique_zero"],
                            "Z_centre_float": [mp.nstr(mp.re(mp.mpc(z_over_km)), 30), mp.nstr(mp.im(mp.mpc(z_over_km)), 30)],
                            "Z_halfwidth_rel": "1e-20", "contraction": kr["contraction"],
                            "enclosure_over_kappa_minus": {"re": K.real.str(30), "im": K.imag.str(30)},
                            "float_value_in_Z": float_in_Z,
                            "float_value_abs_dist_to_Z_centre_over_kappa_minus": mp.nstr(abs(mp.mpc(mp.mpf(z0.real.mid().str(80, radius=False)), mp.mpf(z0.imag.mid().str(80, radius=False))) - zc) / mp.mpf(km), 5),
                            "Krawczyk_centre_in_winding_box": bool(in_box)},
               "prec_bits": prec}
        if printed is not None:
            pr = acb(arb(printed[0]), arb(printed[1]))
            dre, dim = abs(K.real - pr.real), abs(K.imag - pr.imag)
            out["printed_value"] = printed
            out["printed_within_5e-11"] = bool(dre < arb("5e-11") and dim < arb("5e-11"))
            out["printed_abs_diff_upper"] = [float(dre.upper()), float(dim.upper())]
        out["pass"] = bool(wres["certified"] and wres["winding"] == 1 and kr["certified_unique_zero"] and in_box
                           and out.get("printed_within_5e-11", True))
    out["timing_s"] = {"wall": round(time.time() - t0, 1), "cpu": round(time.process_time() - c0, 1)}
    return out


def run_kat6(write=True):
    cases = [("PA", 1, ("3.242616412611042303686208", "-0.7958326323112821365072551"),
              "Table I l=1 PS", ("3.2426164126", "-0.7958326323")),
             ("PB", 0, ("0", "-0.9093809789948054555249"), "P_B l=0 NE n=0", None)]
    res = {"KAT": "KAT-6", "commit": git_commit(), "cases": []}
    for name, l, z, label, printed in cases:
        mp.mp.dps = 60
        res["cases"].append(kat6_case(name, l, z, label, printed))
        print(json.dumps(res["cases"][-1], indent=1)[:2000], flush=True)
    res["pass"] = all(c["pass"] for c in res["cases"])
    if write:
        (CERT_DIR / "kat6_winding_krawczyk.json").write_text(json.dumps(res, indent=1) + "\n")
    return res


if __name__ == "__main__":
    r = run_kat6()
    print("KAT-6 pass:", r["pass"])
