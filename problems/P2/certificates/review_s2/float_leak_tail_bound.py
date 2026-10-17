"""Size of the float leak in m3.tail_bound (m3.py:398-399) (reviewer script, 2026-10-17).

tail_bound accumulates Smin/Smax through float(arb.lower()) / float(arb.upper()), which round to the
NEAREST double, then subtracts/adds 1e-30, far below one ulp of terms of size 5-40.  This script
recomputes the same bound with exact Arb min/max (no float) and reports how far the author's d can
undershoot, then re-evaluates the author's tail check with the disc enlarged by that amount.
Run (repo root):  .venv/bin/python problems/P2/certificates/review_s2/float_leak_tail_bound.py
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from flint import acb, arb  # noqa: E402

from p2 import m3  # noqa: E402
from p2.arbutil import precision  # noqa: E402
from p2.params import Point  # noqa: E402


def tail_bound_exact(pt, side, S, d0):
    """m3.tail_bound with exact Arb min/max in place of the float round trip."""
    xP = m3.rstar(pt, acb(pt.rP))
    rh = pt.rp if side < 0 else pt.rc
    others = [(Ai, ri) for Ai, ri in zip(pt.Acoef, pt.roots) if not ri.overlaps(rh)]
    Ah = next(Ai for Ai, ri in zip(pt.Acoef, pt.roots) if ri.overlaps(rh))
    disc = acb(rh + arb(0, 1) * d0, arb(0, 1) * d0)
    Smin, Smax, worst_ulp = arb(0), arb(0), 0.0
    for Ai, ri in others:
        lo_, hi_ = m3.dist_bounds(disc, ri)
        t1, t2 = Ai * lo_.log(), Ai * hi_.log()
        lo_exact = t1.lower() if t1.lower() < t2.lower() else t2.lower()
        hi_exact = t1.upper() if t1.upper() > t2.upper() else t2.upper()
        # what the author's code adds: arb(float(...)) -/+ 1e-30
        lo_auth = arb(min(float(t1.lower()), float(t2.lower()))) - arb("1e-30")
        hi_auth = arb(max(float(t1.upper()), float(t2.upper()))) + arb("1e-30")
        worst_ulp = max(worst_ulp, float(abs(lo_auth - lo_exact).upper()), float(abs(hi_auth - hi_exact).upper()))
        Smin += arb(lo_exact)
        Smax += arb(hi_exact)
    xS = xP + arb(S) * arb(2).sqrt() / 2
    d = ((xS - arb(Smin.lower())) / Ah).exp() if side < 0 else ((xS - arb(Smax.upper())) / Ah).exp()
    return arb(d.upper()), worst_ulp


pt = Point("PB", 256)
with precision(m3.PREC):
    geo, rays, tails = m3.geometry_certificate(pt)
    for side, S, d0 in ((-1, m3.S_LEFT, "0.06"), (1, m3.S_RIGHT, "0.6")):
        d_auth = m3.tail_bound(pt, side, S, rays[side].at(arb(S)), arb(d0))
        d_ex, ulp = tail_bound_exact(pt, side, S, arb(d0))
        rel = float(((d_ex - d_auth) / d_ex).upper())
        print(f"side {side:+d}: d_author = {d_auth.str(20)}  d_exact = {d_ex.str(20)}  "
              f"max per-term float-rounding error = {ulp:.2e}  (d_exact - d_author)/d = {rel:.2e}")
        # author's tail check with the disc enlarged by 1e-12 (far more than the leak)
        dd = arb((d_ex * (1 + arb("1e-12"))).upper())
        for l in (20, 50, 100):
            wl, a1, a2, b1 = m3.window(pt, l)
            W2 = m3.w2_hull(a1, a2, b1)
            q_low = m3.tail_check(pt, l, side, S, dd, W2)
            print(f"   l={l}: tail lower bound of E with the enlarged disc = {q_low.str(8)}")
