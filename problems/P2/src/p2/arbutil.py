"""Small helpers around python-flint's Arb types (ball arithmetic, outward rounding).

Comparison convention (Arb): ``x < y`` is True only if it holds for every pair of
points of the two balls ("certainly true"); False means "not certified".
"""
from __future__ import annotations

import os
import subprocess
import time
from contextlib import contextmanager
from fractions import Fraction

from flint import acb, arb, ctx, fmpq

DEFAULT_PREC = 256


@contextmanager
def precision(prec: int = DEFAULT_PREC):
    """Temporarily set the Arb working precision (bits)."""
    old = ctx.prec
    ctx.prec = prec
    try:
        yield
    finally:
        ctx.prec = old


def q(x) -> fmpq:
    """Exact rational from int / str ('3/50') / Fraction / fmpq."""
    if isinstance(x, fmpq):
        return x
    if isinstance(x, Fraction):
        return fmpq(x.numerator, x.denominator)
    if isinstance(x, int):
        return fmpq(x)
    if isinstance(x, str):
        f = Fraction(x)
        return fmpq(f.numerator, f.denominator)
    raise TypeError(f"not an exact rational: {x!r}")


def A(x) -> arb:
    """Exact-as-possible arb from exact rationals (never from a decimal float)."""
    if isinstance(x, arb):
        return x
    if isinstance(x, float):
        raise TypeError("refusing to build an arb from a float; pass an exact rational")
    return arb(q(x))


def lo(x: arb) -> float:
    """Lower endpoint as a float rounded down (for reporting only)."""
    return float(x.lower())


def hi(x: arb) -> float:
    return float(x.upper())


def abs_upper(z) -> arb:
    """Upper bound of |z| as an exact arb (z arb or acb)."""
    return arb(abs(z).upper())


def abs_lower(z) -> arb:
    return arb(abs(z).lower())


def ball_str(x, digits: int = 30) -> str:
    """Arb's own decimal rendering [mid +/- rad]; Arb guarantees the printed ball contains x."""
    return x.str(digits, radius=True)


def mid_str(x, digits: int = 20) -> str:
    """Rounded midpoint, for human reading only (not a bound)."""
    return x.mid().str(digits, radius=False)


def acb_box(re_lo, re_hi, im_lo, im_hi) -> acb:
    """Complex box [re_lo, re_hi] x [im_lo, im_hi] (exact rational/arb endpoints)."""
    rl, rh, il, ih = (A(v) for v in (re_lo, re_hi, im_lo, im_hi))
    return acb((rl + rh) / 2 + arb(0, 1) * ((rh - rl) / 2),
               (il + ih) / 2 + arb(0, 1) * ((ih - il) / 2))


def git_commit() -> str:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                             cwd=os.path.dirname(os.path.abspath(__file__)), check=True)
        sha = out.stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--", "."], capture_output=True,
                               text=True, cwd=os.path.dirname(os.path.abspath(__file__)), check=True)
        return sha + ("+dirty" if dirty.stdout.strip() else "")
    except (OSError, subprocess.CalledProcessError):  # pragma: no cover - outside git
        return "unknown"


class Timer:
    """Wall and process-CPU time of a block."""

    def __enter__(self):
        self.w0, self.c0 = time.time(), time.process_time()
        return self

    def __exit__(self, *exc):
        self.wall = time.time() - self.w0
        self.cpu = time.process_time() - self.c0
        return False
