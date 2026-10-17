"""Tests of p4.validated.ballutil: balls that contain float intervals
(review docs/review-2026-10-17-p4-layer-engine-adversarial.md, item 2f).

The float-midpoint ball arb((lo + hi)/2, (hi - lo)/2) can miss an endpoint;
strip_ball / recentred_delta / half_width / radius_cover must contain what
they promise at every working precision, and must return the old ball,
bit for bit, whenever the old ball already contained it."""
import random

import pytest
from flint import arb

from p4.validated.arbseries import precision
from p4.validated.ballutil import half_width, radius_cover, recentred_delta, strip_ball
from p4.validated.tc_tm import LC, TM

RELS = (1e-15, 3e-15, 1e-14, 1e-13, 1e-12, 1e-10, 1e-7, 1e-4, 1e-2, 1.0)
PRECS = (30, 53, 128, 320)


def _intervals(n, seed):
    """Random non-dyadic intervals: lo of random magnitude and sign, relative
    widths from 1e-15 to 1, plus [0, hi] and degenerate ones."""
    rnd = random.Random(seed)
    out = []
    for _ in range(n):
        lo = rnd.uniform(0.1, 10.0) * 10.0 ** rnd.randint(-8, 4) * rnd.choice((-1, 1))
        rel = rnd.choice(RELS) * rnd.uniform(0.5, 2.0)
        hi = lo + abs(lo) * rel
        if hi <= lo:
            hi = lo
        out.append((lo, hi))
    out += [(0.0, rnd.uniform(1e-9, 1.0)) for _ in range(20)]
    out += [(x, x) for x in (0.1, 1.0 / 3.0, 224.0, -2.95)]
    return out


def _old(lo, hi):
    return arb((lo + hi) / 2, (hi - lo) / 2)


def _same(a, b):
    return a.mid() == b.mid() and a.rad() == b.rad()


@pytest.mark.parametrize("prec", PRECS)
def test_strip_ball_contains_the_endpoints(prec):
    misses = 0
    with precision(prec):
        for lo, hi in _intervals(4000, prec):
            b = strip_ball(lo, hi)
            assert b.contains(arb(lo)) and b.contains(arb(hi)), (lo, hi)
            old = _old(lo, hi)
            if old.contains(arb(lo)) and old.contains(arb(hi)):
                assert _same(b, old)                     # unchanged when no miss
            else:
                misses += 1
    assert misses > 0                    # the old ball does miss (at every precision)


@pytest.mark.parametrize("prec", PRECS)
def test_recentred_delta_and_half_width(prec):
    with precision(prec):
        for lo, hi in _intervals(3000, 7 + prec):
            mid = arb((lo + hi) / 2)
            d = recentred_delta(lo, hi, mid)
            assert d.contains(arb(lo) - mid) and d.contains(arb(hi) - mid), (lo, hi)
            old = arb(0, (hi - lo) / 2)
            if old.contains(arb(lo) - mid) and old.contains(arb(hi) - mid):
                assert _same(d, old)
            w = half_width(lo, hi, mid)
            assert abs(arb(lo) - mid) <= w and abs(arb(hi) - mid) <= w, (lo, hi)
            wold = arb((hi - lo) / 2)
            if abs(arb(lo) - mid) <= wold and abs(arb(hi) - mid) <= wold:
                assert _same(w, wold)


def test_radius_cover():
    rnd = random.Random(3)
    with precision(320):
        for _ in range(3000):
            h = rnd.uniform(1e-4, 0.1)
            n = rnd.randint(1, 600)
            exact = arb(h) / (2 * n)
            r = radius_cover(h / (2 * n), exact)
            assert arb(0, r).contains(exact) and arb(0, r).contains(-exact)
            if arb(0, h / (2 * n)).contains(exact):
                assert r == h / (2 * n)
        # a double radius that rounds below the exact one is replaced
        exact = arb(1) / 3 + arb(0, 1e-30)
        r = radius_cover(float(arb(1) / 3) * (1 - 2 ** -40), exact)
        assert arb(0, r).contains(exact)


def test_tm_eval_contains_the_endpoint_values():
    """TM.eval and LC.eval on strips of relative width down to 1e-15: the
    enclosure contains the polynomial's values at both endpoints."""
    rnd = random.Random(11)
    with precision(320):
        for k in range(300):
            R = 0.05
            c = [arb(rnd.uniform(-3, 3)) for _ in range(rnd.choice((2, 3, 12, 40)))]
            tm = TM(c, R)
            lo = rnd.uniform(1e-6, 0.04)
            hi = min(R, lo * (1 + rnd.choice(RELS[:6]) * rnd.uniform(0.5, 2.0)))
            v = tm.eval(lo, hi)
            for x in (lo, hi):
                px = arb(0)
                for cc in reversed(c):
                    px = px * arb(x) + cc
                assert v.contains(px), (k, lo, hi)
            lc = LC(tm, -1 if k % 2 else 2)
            lv = lc.eval(lo, hi)
            for x in (lo, hi):
                px = arb(0)
                for cc in reversed(c):
                    px = px * arb(x) + cc
                assert lv.contains(px * arb(x) ** lc.off), (k, lo, hi)
