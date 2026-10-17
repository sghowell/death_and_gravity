"""Balls covering float intervals, for the P4 certificates.

The construction  arb((lo + hi) / 2, (hi - lo) / 2)  rounds the midpoint and
the half-width to doubles, so the ball can miss lo or hi by about one ulp
(review docs/review-2026-10-17-p4-layer-engine-adversarial.md, item 2f).
Arb's upward rounding of the radius to a 30-bit mantissa usually hides
this, but not on strips narrower than about 1e-7 relative, nor when the
double radius happens to need no rounding.  The helpers below keep the old
ball whenever it provably contains what it must contain, so every result
that did not depend on a miss stays bit-identical, and otherwise fall back
to an exact hull.

Soundness at any working precision: arb(float) is exact, arb_contains is
exact or conservative (it may answer False for a contained ball, never True
for a non-contained one), and arb_union and ball subtraction round outward.
"""
from __future__ import annotations

from flint import arb


def strip_ball(lo, hi):
    """A ball containing the closed interval [lo, hi] (floats or ints,
    lo <= hi): the float-midpoint ball arb((lo + hi)/2, (hi - lo)/2) if it
    contains both arb(lo) and arb(hi) (a ball is convex, so it then contains
    [lo, hi]), else the hull arb(lo).union(arb(hi))."""
    b = arb((lo + hi) / 2, (hi - lo) / 2)
    a_lo, a_hi = arb(lo), arb(hi)
    if b.contains(a_lo) and b.contains(a_hi):
        return b
    return a_lo.union(a_hi)


def recentred_delta(lo, hi, mid):
    """A ball containing x - mid for every x in [lo, hi], mid an arb: the
    old arb(0, (hi - lo)/2) if it contains arb(lo) - mid and arb(hi) - mid,
    else arb(lo).union(arb(hi)) - mid."""
    d = arb(0, (hi - lo) / 2)
    a_lo, a_hi = arb(lo), arb(hi)
    if d.contains(a_lo - mid) and d.contains(a_hi - mid):
        return d
    return a_lo.union(a_hi) - mid


def half_width(lo, hi, mid):
    """An arb w >= |x - mid| for every x in [lo, hi] (a magnitude for
    Taylor-remainder bounds): the old exact point arb((hi - lo)/2) if both
    |arb(lo) - mid| and |arb(hi) - mid| are provably <= it, else the upper
    bound of max(|arb(lo) - mid|, |arb(hi) - mid|) as an exact arb."""
    w = arb((hi - lo) / 2)
    d_lo, d_hi = abs(arb(lo) - mid), abs(arb(hi) - mid)
    if d_lo <= w and d_hi <= w:
        return w
    return arb(d_lo.max(d_hi).upper())


def radius_cover(r, exact):
    """A radius (for arb(0, .)) at least the exact half-width `exact` (an
    arb): the old double r if arb(0, r) contains the ball `exact`, else
    exact.upper() (an arb, which arb(0, .) rounds up again)."""
    if arb(0, r).contains(exact):
        return r
    return exact.upper()
