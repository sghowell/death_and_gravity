"""Opt-in exact SymPy 1.14.0 runner for P8's Gaussian-rational regressions.

Run ``python scripts/p8_exact_regression.py --self-check`` or append pytest
arguments. This changes no source, certificate, dependency or stored formula.
Only univariate Gaussian polynomials with a separate real/imaginary unit
phase per operand descend to the corresponding real coefficient domain.
Both cofactor identities are checked on every descent; all other cases use
the original implementation. The temporary method replacement is restored.
"""

import sys
from collections import Counter
from contextlib import contextmanager

import sympy as sp
from sympy.polys.domains import QQ, QQ_I, ZZ, ZZ_I
from sympy.polys.rings import PolyElement, ring

ORIGINAL_GCD = PolyElement._gcd


def phase(poly):
    if all(not coefficient.y for coefficient in poly.values()):
        return 0
    if all(not coefficient.x for coefficient in poly.values()):
        return 1
    return None


def descended_gcd(left, right, counts):
    parent = left.ring
    domain = parent.domain
    if (domain not in (QQ_I, ZZ_I) or parent.ngens != 1
            or not left or not right or right.ring != parent):
        counts["domain_fallback"] += 1
        return ORIGINAL_GCD(left, right)
    lp, rp = phase(left), phase(right)
    if lp is None or rp is None:
        counts["mixed_fallback"] += 1
        return ORIGINAL_GCD(left, right)
    real_ring = parent.clone(domain=QQ if domain == QQ_I else ZZ)
    a = real_ring.from_dict({m: c.y if lp else c.x for m, c in left.items()})
    b = real_ring.from_dict({m: c.y if rp else c.x for m, c in right.items()})
    gcd, ca, cb = ORIGINAL_GCD(a, b)

    def restore(poly, imaginary=0):
        return parent.from_dict({m: domain.dtype(0, c) if imaginary else domain.dtype(c, 0)
                                 for m, c in poly.items()})

    result = restore(gcd), restore(ca, lp), restore(cb, rp)
    if result[0]*result[1] != left or result[0]*result[2] != right:
        raise ArithmeticError("Exact Gaussian GCD descent failed a cofactor identity")
    counts["exact_descent"] += 1
    return result


def self_check():
    """Compare entire normalized tuples against the untouched Gaussian code."""
    if sp.__version__ != "1.14.0":
        raise RuntimeError("This opt-in private-API adapter is audited only for SymPy 1.14.0")
    counts, compared = Counter(), 0
    for domain in (ZZ_I, QQ_I):
        _, x = ring("x", domain)
        imaginary = domain.dtype(0, 1)
        for left_phase in (1, -1, imaginary, -imaginary):
            for right_phase in (1, -1, imaginary, -imaginary):
                for seed in range(4):
                    denominator = 1 if domain == ZZ_I else seed+1
                    a = left_phase*(x*x+seed+1)*(2*x+3)/denominator
                    b = right_phase*(x*x+seed+1)*(3*x-2)/(denominator+1 if domain == QQ_I else 1)
                    expected = ORIGINAL_GCD(a, b)
                    actual = descended_gcd(a, b, counts)
                    if actual != expected:
                        raise ArithmeticError("Gaussian descent disagrees with the original normalized GCD tuple")
                    compared += 1
    return {"original_tuple_comparisons": compared, **dict(counts)}


@contextmanager
def exact_runner():
    if sp.__version__ != "1.14.0" or PolyElement._gcd is not ORIGINAL_GCD:
        raise RuntimeError("Unexpected SymPy version or existing GCD replacement")
    counts = Counter()

    def replacement(left, right):
        return descended_gcd(left, right, counts)

    PolyElement._gcd = replacement
    try:
        yield counts
    finally:
        PolyElement._gcd = ORIGINAL_GCD


def main(args=None):
    arguments = sys.argv[1:] if args is None else args
    print("Exact GCD adapter self-check:", self_check(), flush=True)
    if arguments == ["--self-check"]:
        return 0
    import pytest

    with exact_runner() as counts:
        try:
            return pytest.main(arguments)
        finally:
            print("Exact GCD runner counters:", dict(counts), flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
