"""S3 review: a rigorous recount of W_l around the FULL boundary of R'_l (no conjugation symmetry),
with the reviewer's own counting layer.

What is reused from the author: only the zsol BOX enclosure, for the sup bound M on each disc
(its rigour is reviewed separately; box_checks.py tests it).  Everything else is the reviewer's:
  * the plan: segments of length <= kappa_-/2 (halved near E_+ u E_c), disc radius
    rho = (5/4) * length, so every evaluation point is within 0.4 rho of the centre; the box
    [c +- rho]^2 is certified to avoid E_+ u E_c (corner distance sqrt(2) rho);
  * the point values: S1's arbw.wronskian (an independent code from zsol) at the K nodes
    c + rho/2 e^{2 pi i j/K}, with its Frobenius minimum term count raised as in C0;
  * the coefficients: DFT with the aliasing bound M x^K/((1 - x^K) rho^k), x = 1/2;
  * the enclosure of W on a piece: the Taylor polynomial on the piece (ball Horner) plus
    M (s/rho)^K/(1 - s/rho), s = max |w - c| on the piece;
  * pieces: length <= kappa_-/32 initially, bisected until the enclosure excludes 0; increment =
    principal Arg of the ratio of the end enclosures, certified below pi.
The winding is (sum of all increments)/(2 pi), certified when within 1/4 of an integer.
Usage: nice -n 10 python recount.py 1,2 --workers 3     (writes recount.json; logs/recount_l{l}.jsonl)
"""
from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mpr
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

from flint import acb, arb, fmpq  # noqa: E402

CERT = HERE.parent / "part3"
LOGS = HERE / "logs"
PREC = 128
XN = arb(1) / 2            # node radius / rho
_G = {}


def _init():
    from p2.params import Point
    import p2.arbw as arbw_mod
    _G["pt"] = Point("PB", 256)
    orig = arbw_mod.frobenius

    def frob(ptx, which, w, L, t, min_terms=20):
        kap = ptx.kp if which == "+" else ptx.kc
        ra = float((abs(acb(w)) / (2 * kap)).upper())
        return orig(ptx, which, w, L, t, min_terms=max(min_terms, int(2 * ra) + 8))

    arbw_mod.frobenius = frob
    _G["arbw"] = arbw_mod.wronskian


def poles_clear(pt, c: acb, rho: arb) -> bool:
    s = arb(2).sqrt() * rho
    for kap in (pt.kp, pt.kc):
        for n in range(1, 40):
            if not bool(abs(c - acb(0, -n * kap)) > s):
                return False
    return bool(c.imag > -39 * pt.kp)


def plan(pt, a: Fr, b: Fr, e: Fr):
    km = Fr(float(pt.km.mid())).limit_denominator(2 ** 40)
    corners = [(Fr(0), -b), (a, -b), (a, e), (-a, e), (-a, -b), (Fr(0), -b)]
    segs = []
    for (x0, y0), (x1, y1) in zip(corners[:-1], corners[1:]):
        L = abs(x1 - x0) + abs(y1 - y0)
        hmax = km / 2
        t = Fr(0)
        while t < L:
            h = min(hmax, L - t)
            while True:
                t1 = min(t + h, L)
                p0 = (x0 + (x1 - x0) * t / L, y0 + (y1 - y0) * t / L)
                p1 = (x0 + (x1 - x0) * t1 / L, y0 + (y1 - y0) * t1 / L)
                c = ((p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2)
                rho = Fr(5, 4) * (t1 - t)
                cc = acb(arb(fmpq(c[0].numerator, c[0].denominator)), arb(fmpq(c[1].numerator, c[1].denominator)))
                if poles_clear(pt, cc, arb(fmpq(rho.numerator, rho.denominator)) * Fr(5, 4)):
                    break
                h = h / 2
            segs.append({"p0": [str(p0[0]), str(p0[1])], "p1": [str(p1[0]), str(p1[1])],
                         "c": [str(c[0]), str(c[1])], "rho": str(rho)})
            t = t1
    return segs


def _q(s: str) -> arb:
    f = Fr(s)
    return arb(fmpq(f.numerator, f.denominator))


def segment(args):
    from p2 import zsol
    from p2.arbutil import precision
    l, sg = args
    pt = _G["pt"]
    arbw = _G["arbw"]
    t0 = time.process_time()
    with precision(PREC):
        c = acb(_q(sg["c"][0]), _q(sg["c"][1]))
        rho = _q(sg["rho"])
        p0 = acb(_q(sg["p0"][0]), _q(sg["p0"][1]))
        p1 = acb(_q(sg["p1"][0]), _q(sg["p1"][1]))
        if not poles_clear(pt, c, rho):
            return {"ok": False, "seg": sg, "error": "poles"}
        box = acb(c.real + arb(0, 1) * rho, c.imag + arb(0, 1) * rho)
        _, Wb = zsol.wronskian_z(pt, l, c, box, theta=0.25)
        M = arb(abs(Wb).upper())
        Wc = arbw(pt, l, c)
        if not (abs(Wc).lower() > 0) or not M.is_finite():
            return {"ok": False, "seg": sg, "error": "centre/M"}
        ratio = float((M / abs(Wc).lower()).upper())
        K = 16
        while float(ratio) * 0.5 ** K / (1 - 0.5 ** K) > 1e-9 and K < 240:
            K += 8
        rn = XN * rho
        vals = []
        for j in range(K):
            th = 2 * arb.pi() * j / K
            vals.append(arbw(pt, l, c + rn * acb(th.cos(), th.sin())))
        alias = M * XN ** K / (1 - XN ** K)
        coef = []
        for k in range(K):
            S = acb(0)
            for j in range(K):
                th = -2 * arb.pi() * ((j * k) % K) / K
                S += vals[j] * acb(th.cos(), th.sin())
            err = alias / rho ** k
            coef.append(S / K / rn ** k + acb(arb(0, 1) * err, arb(0, 1) * err))

        def enc(z):
            d = z - c
            s = arb(abs(d).upper())
            if not (s < rho):
                raise ArithmeticError("outside")
            acc = acb(0)
            for cf in reversed(coef):
                acc = acc * d + cf
            rem = M * (s / rho) ** K / (1 - s / rho)
            return acc + acb(arb(0, 1) * rem, arb(0, 1) * rem)

        length = float(abs(p1 - p0).upper())
        n0 = max(1, math.ceil(length / (float(pt.km.mid()) / 32)))
        stack = [(p0 + (p1 - p0) * k / n0, p0 + (p1 - p0) * (k + 1) / n0, 0) for k in range(n0)][::-1]
        total = arb(0)
        npieces = 0
        minabs = None
        while stack:
            z0, z1, dep = stack.pop()
            mid = (z0 + z1) / 2
            hr, hi = abs((z1 - z0).real) / 2, abs((z1 - z0).imag) / 2
            seg_ball = acb(mid.real + arb(0, 1) * hr, mid.imag + arb(0, 1) * hi)
            E = enc(seg_ball)
            if E.real.contains(0) and E.imag.contains(0):
                if dep >= 16:
                    return {"ok": False, "seg": sg, "error": "piece contains 0"}
                stack.append((mid, z1, dep + 1))
                stack.append((z0, mid, dep + 1))
                continue
            inc = (enc(z1) / enc(z0)).arg()
            if not (abs(inc) < arb.pi()):
                return {"ok": False, "seg": sg, "error": "increment"}
            total += inc
            npieces += 1
            lb = max(0.0 if E.real.contains(0) else float(abs(E.real).lower()),
                     0.0 if E.imag.contains(0) else float(abs(E.imag).lower()))
            minabs = lb if minabs is None else min(minabs, lb)
    return {"ok": True, "l": l, "seg": sg, "increment": total.str(25), "K": K, "M_over_absWc": ratio,
            "pieces": npieces, "min_abs_lower": minabs, "cpu": round(time.process_time() - t0, 2)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ls")
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args()
    from p2.params import Point
    from p2.arbutil import precision
    pt = Point("PB", 256)
    LOGS.mkdir(exist_ok=True)
    outp = HERE / "recount.json"
    out = json.loads(outp.read_text()) if outp.exists() else {"per_l": {}}
    for l in [int(x) for x in a.ls.split(",")]:
        reg = json.loads((CERT / f"winding_PB_l{l}.json").read_text())["region"]
        A_, B_, E_ = Fr(reg["a"]), Fr(reg["b"]), Fr(reg["e"])
        lg = json.loads((HERE.parent / "lg" / f"lg_PB_l{l}.json").read_text())
        assert A_ == Fr(lg["W_l_exact"])
        segs = plan(pt, A_, B_, E_)
        f = LOGS / f"recount_l{l}.jsonl"
        done = {}
        if f.exists():
            for line in f.read_text().splitlines():
                r = json.loads(line)
                if r["ok"]:
                    done[json.dumps(r["seg"], sort_keys=True)] = r
        todo = [(l, s) for s in segs if json.dumps(s, sort_keys=True) not in done]
        print(f"l={l}: {len(segs)} segments, {len(todo)} to run", flush=True)
        t0 = time.time()
        ctx = mpr.get_context("spawn")
        with ctx.Pool(a.workers, initializer=_init) as pool:
            for r in pool.imap_unordered(segment, todo):
                with open(f, "a") as fh:
                    fh.write(json.dumps(r) + "\n")
                if not r["ok"]:
                    print("FAILED", r, flush=True)
                else:
                    done[json.dumps(r["seg"], sort_keys=True)] = r
        # chain check and sum
        recs = [done.get(json.dumps(s, sort_keys=True)) for s in segs]
        complete = all(r is not None for r in recs)
        chain = all(Fr(segs[i]["p1"][0]) == Fr(segs[i + 1]["p0"][0]) and Fr(segs[i]["p1"][1]) == Fr(segs[i + 1]["p0"][1])
                    for i in range(len(segs) - 1)) and segs[0]["p0"] == segs[-1]["p1"]
        with precision(PREC):
            total = sum((arb(r["increment"]) for r in recs if r), arb(0))
            wn = total / (2 * arb.pi())
            k = round(float(wn.mid()))
            cert = bool(complete and chain and abs(wn - k) < arb("0.25"))
        out["per_l"][str(l)] = {"segments": len(segs), "complete": complete, "closed_chain": chain,
                                "winding_ball": wn.str(12), "certified": cert, "winding": k if cert else None,
                                "K_max": max(r["K"] for r in recs if r), "min_abs_W_lower": min(r["min_abs_lower"] for r in recs if r),
                                "cpu_s": round(sum(r["cpu"] for r in recs if r), 1), "wall_s": round(time.time() - t0, 1)}
        print(l, out["per_l"][str(l)], flush=True)
        outp.write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
