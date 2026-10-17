"""S3 review: the C0 deviation (3e7fca05) and the l = 0 enclosure (C6 / P2-K6 amendment).

C0a  The ORIGINAL C0 of 89a2336b, replayed: S1 arbw.wronskian called unwrapped at 256 bits at the 24
     pre-registered points; records which points raise, arbw's relative radius, and the overlap
     with zsol.  (Shows what the first run met and whether "overlap" was meaningful.)
C0b  The test point dropped from tests/test_part3.py in 3e7fca05, (l, Re, Im) = (7, 40.0, -0.55),
     run through the current C0 code (part3.zsol_vs_arbw).
K1   part3.krawczyk_z re-run; compared with the stored certificates/part3/krawczyk_PB_l0.json.
K2   Elementary existence check, independent of the Krawczyk code: W_0 is real on the imaginary axis
     (conjugation symmetry), so a certified sign change of Re W_0(i y kappa_-) between the two
     endpoints of the stored enclosure's imaginary interval proves a zero inside it (intermediate
     value theorem).  Done with zsol (S3) and with arbw (S1), at 256 bits.
K3   The reviewer's float W (indep_w.py): the zero located by secant iteration on the imaginary
     axis in double precision.
Writes c0_and_krawczyk.json.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))

from flint import acb, arb, fmpq  # noqa: E402

from p2 import part3, zsol  # noqa: E402
from p2.arbutil import precision  # noqa: E402
from p2.arbw import wronskian as arbw_w  # noqa: E402
from p2.params import Point  # noqa: E402

CERT = HERE.parent / "part3"


def c0a():
    pt = Point("PB", 256)
    out = []
    for l, re_k, im_k in part3.c0_cases("PB"):
        rec = {"l": l, "w_kappa": [re_k, im_k]}
        t0 = time.time()
        with precision(256):
            w = acb(arb(part3._kq(pt, re_k)), arb(part3._kq(pt, im_k)))
            Wz = zsol.wronskian_z(pt, l, w, theta=part3.THETA_POINT)[0]
            try:
                Wa = arbw_w(pt, l, w)
                rec["arbw_raised"] = None
                rec["arbw_relrad"] = float(abs(Wa).rad() / abs(Wa).mid()) if abs(Wa).mid() != 0 else float("inf")
                rec["overlap"] = bool(Wz.overlaps(Wa))
                rec["meaningful_at_1e-10"] = rec["arbw_relrad"] < 1e-10
            except ArithmeticError as ex:
                rec["arbw_raised"] = str(ex)
        rec["wall_s"] = round(time.time() - t0, 1)
        out.append(rec)
        print(json.dumps(rec), flush=True)
    return out


def c0b():
    t0 = time.time()
    r = part3.zsol_vs_arbw("PB", [(7, 40.0, -0.55)])
    r[0]["wall_s"] = round(time.time() - t0, 1)
    print(json.dumps(r), flush=True)
    return r


def k1():
    t0 = time.time()
    kr = part3.krawczyk_z("PB", 0, ("0", "-0.9093809789948054555249"))
    st = json.loads((CERT / "krawczyk_PB_l0.json").read_text())
    same = all(kr[k] == st[k] for k in ("certified_unique_zero", "enclosure_over_kappa_minus", "contraction"))
    out = {"rerun": kr, "identical_to_stored": same, "wall_s": round(time.time() - t0, 1)}
    print(json.dumps(out), flush=True)
    return out


def k2():
    st = json.loads((CERT / "krawczyk_PB_l0.json").read_text())
    ys = st["enclosure_over_kappa_minus"]["im"]
    pt = Point("PB", 256)
    res = {}
    with precision(256):
        yb = arb(ys)
        lo, hi = yb.lower(), yb.upper()        # exact arbs (endpoints of the stored interval)
        out = []
        for name in ("zsol", "arbw"):
            vals = []
            for y in (lo, hi):
                w = acb(arb(0), y * pt.km)     # w = i y kappa_- (y < 0); a tiny ball from kappa_-
                if name == "zsol":
                    W = zsol.wronskian_z(pt, 0, w, theta=part3.THETA_POINT)[0]
                else:
                    W = arbw_w(pt, 0, w)
                vals.append(W)
            s0 = vals[0].real
            s1 = vals[1].real
            change = bool((s0 > 0 and s1 < 0) or (s0 < 0 and s1 > 0))
            out.append({"code": name, "Re_W_at_lower_end": s0.str(5), "Re_W_at_upper_end": s1.str(5),
                        "Im_W_contains_0": [bool(v.imag.contains(0)) for v in vals],
                        "certified_sign_change": change})
        res["interval_im_over_kappa"] = ys
        res["checks"] = out
        res["note"] = ("w = i y kappa_- is evaluated as a ball (kappa_- is a ball); a sign change of the"
                       " real part between the two balls proves a zero of the real-valued W_0(i y kappa_-)"
                       " for y between them, i.e. inside the stored interval up to the kappa_- radius")
    print(json.dumps(res), flush=True)
    return res


def k3():
    from indep_w import Geo, W_l
    g = Geo()
    f = lambda y: W_l(g, 0, [1j * y * g.km], rtol=1e-13)[0]
    y0, y1 = -0.9093, -0.9095
    f0, f1 = f(y0), f(y1)
    for _ in range(30):
        y2 = y1 - f1.real * (y1 - y0) / (f1.real - f0.real)
        y0, f0 = y1, f1
        y1, f1 = y2, f(y2)
        if abs(y1 - y0) < 1e-16:
            break
    res = {"float_zero_im_over_kappa": y1, "W_imag_at_zero": f1.imag, "W_real_at_zero": f1.real,
           "stored_mid": -0.9093809789948054555249, "difference": y1 - (-0.9093809789948054555249)}
    print(json.dumps(res), flush=True)
    return res


def main():
    which = sys.argv[1:] or ["C0a", "C0b", "K1", "K2", "K3"]
    outp = HERE / "c0_and_krawczyk.json"
    out = json.loads(outp.read_text()) if outp.exists() else {}
    fns = {"C0a": c0a, "C0b": c0b, "K1": k1, "K2": k2, "K3": k3}
    for k in which:
        out[k] = fns[k]()
        outp.write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
