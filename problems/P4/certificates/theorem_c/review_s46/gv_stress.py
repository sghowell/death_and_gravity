#!/usr/bin/env python3
"""S4-6 review: randomized containment test of tc_layer_centre.GV (graded
values: exact exponents, folding above KMAX, the Neumann inverse with its
remainder, refusal of nonzero negative exponents).

Each case draws GVs with random exponents and complex ball coefficients on a
random g-hull (often touching 0), then picks a point g in the hull and a
point inside every coefficient ball and checks that the GV's value() ball
contains the exact (high-precision) value of sum_a g^a c_a, for x, x*y, x^n,
x.inv() and (x*y).inv().  Negative exponents are given exactly-zero
coefficients (value() must accept them) or a nonzero one (value() must
refuse).

usage: PYTHONPATH=problems/P4/src python gv_stress.py [cases] [out.json]
"""
import json
import os
import random
import sys

from flint import acb, arb

HERE = os.path.dirname(os.path.abspath(__file__))
HQMAX = float(os.environ.get("GV_HQMAX", "-0.3"))   # log10 of the largest inverse g-hull


def main(n, out):
    from p4.validated.arbseries import precision
    from p4.validated.tc_layer_centre import GV, KMAX
    rng = random.Random(4242)
    stats = dict(cases=0, value=0, product=0, power=0, inverse=0, inverse_refused=0,
                 neg_refused=0, neg_zero_ok=0, failures=[])

    def rand_ball():
        m = complex(rng.uniform(-3, 3), rng.uniform(-3, 3))
        r = 10 ** rng.uniform(-12, -2)
        return acb(arb(m.real, r), arb(m.imag, r)), m, r

    def rand_gv(gh, emin=0, emax=7, nterm=4, q0=False):
        c, pts = {}, {}
        es = rng.sample(range(emin, emax + 1), min(nterm, emax - emin + 1))
        if q0 and 0 not in es:
            es.append(0)
        for e in es:
            b, m, r = rand_ball()
            if q0 and e == 0:
                m = complex(rng.choice([-1, 1]) * rng.uniform(2, 4), rng.uniform(-1, 1))
                b = acb(arb(m.real, r), arb(m.imag, r))
            c[e] = b
            # a point inside the ball (componentwise)
            pts[e] = acb(arb(m.real + rng.uniform(-1, 1) * r * 0.999),
                         arb(m.imag + rng.uniform(-1, 1) * r * 0.999))
        return GV(c, gh), pts

    def exact(pts, g):
        acc = acb(0)
        for e, v in pts.items():
            acc += v * g ** e
        return acc

    with precision(256):
        for _ in range(n):
            stats["cases"] += 1
            h = 10 ** rng.uniform(-6, -0.5)
            a = 0.0 if rng.random() < 0.5 else h * rng.uniform(0, 0.9)
            gh = arb(a).union(arb(h))
            g = arb(rng.uniform(a, h))
            x, px = rand_gv(gh)
            y, py = rand_gv(gh)
            vx = exact(px, g)
            if not x.value().contains(vx):
                stats["failures"].append(("value", a, h))
            stats["value"] += 1
            if not (x * y).value().contains(vx * exact(py, g)):
                stats["failures"].append(("product", a, h))
            stats["product"] += 1
            k = rng.randint(2, 4)
            if not (x ** k).value().contains(vx ** k):
                stats["failures"].append(("power", a, h))
            stats["power"] += 1
            # inverses (q0 bounded away from 0, small g-hull so |u| < 1 often)
            hq = 10 ** rng.uniform(-6, HQMAX)
            gq = arb(0).union(arb(hq))
            q, pq = rand_gv(gq, emin=0, emax=9, nterm=5, q0=True)
            gg = arb(rng.uniform(0, hq))
            try:
                qi = q.inv().value()
                if not qi.contains(1 / exact(pq, gg)):
                    stats["failures"].append(("inverse", hq))
                stats["inverse"] += 1
                r2, pr2 = rand_gv(gq, emin=1, emax=6, nterm=3)
                qq = q * (GV.const(1, gq) + r2)
                ex = exact(pq, gg) * (1 + exact(pr2, gg))
                if not qq.inv().value().contains(1 / ex):
                    stats["failures"].append(("inverse2", hq))
            except ArithmeticError:
                stats["inverse_refused"] += 1
            # negative exponents
            z = GV({-2: acb(0), -1: acb(0), 0: acb(1), 3: acb(2)}, gh)
            try:
                z.value()
                stats["neg_zero_ok"] += 1
            except ArithmeticError:
                stats["failures"].append(("neg_zero_refused", a, h))
            nz = GV({-1: acb(arb(0, 1e-40)), 0: acb(1)}, gh)
            try:
                nz.value()
                stats["failures"].append(("neg_nonzero_accepted", a, h))
            except ArithmeticError:
                stats["neg_refused"] += 1
            # folding bookkeeping: no exponent above KMAX + 1 survives
            if any(e > KMAX + 1 for e in (x * y * x).c):
                stats["failures"].append(("fold", a, h))
    stats["ok"] = not stats["failures"]
    stats["failures"] = stats["failures"][:20]
    json.dump(stats, open(out, "w"), indent=1)
    print(json.dumps(stats, indent=1))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20000,
         sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "gv_stress.json"))
