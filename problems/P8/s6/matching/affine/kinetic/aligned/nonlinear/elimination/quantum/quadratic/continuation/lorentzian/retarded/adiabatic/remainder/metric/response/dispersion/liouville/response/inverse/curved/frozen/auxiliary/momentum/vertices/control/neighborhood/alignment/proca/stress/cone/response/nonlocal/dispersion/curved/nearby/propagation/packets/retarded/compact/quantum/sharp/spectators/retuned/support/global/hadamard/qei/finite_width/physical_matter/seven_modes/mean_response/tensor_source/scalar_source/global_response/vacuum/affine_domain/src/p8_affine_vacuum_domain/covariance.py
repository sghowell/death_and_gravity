"""Unrestricted auxiliary Hessian with general gradients, including null gradients."""

from fractions import Fraction
from functools import cache
from itertools import permutations

import sympy as sp
from p8_affine import connection as old


def hessian(v, p, c, f4):
    if not isinstance(v, (list, tuple)) or len(v) != 4:
        raise ValueError("A four-component covector is required")

    def exact(value):
        if isinstance(value, bool) or not isinstance(
            value, (int, Fraction, sp.Rational)
        ):
            raise TypeError(
                "Only finite exact rational Hessian calibration data are supported"
            )
        return sp.Rational(value)

    v = tuple(exact(value) for value in v)
    p, c, f4 = map(exact, (p, c, f4))
    if p <= 0:
        raise ValueError("The calibrated Palatini coefficient must be positive")
    M = sp.zeros(64)

    def add(a, b, weight):
        M[a, b] += weight
        M[b, a] += weight

    ix = old.index
    signs = old.SIGN
    eps = lambda indices: sp.LeviCivita(*indices)
    # Literal p R_Gamma, no scalar rest frame.
    for a in range(4):
        for b in range(4):
            for e in range(4):
                add(ix(a, e, a), ix(e, b, b), p * signs[b])
                add(ix(a, e, b), ix(e, b, a), -p * signs[b])
    # Literal double-dual affine Einstein curvature coefficient.
    for gam, alpha, a, b in permutations(range(4)):
        for beta, cidx, didx in permutations(i for i in range(4) if i != gam):
            weight = (
                sp.Rational(1, 4)
                * signs[a]
                * signs[gam]
                * eps((gam, alpha, a, b))
                * eps((gam, beta, cidx, didx))
                * c
                * v[alpha]
                * v[beta]
            )
            if weight == 0:
                continue
            for e in range(4):
                add(ix(a, e, cidx), ix(e, b, didx), weight)
                add(ix(a, e, didx), ix(e, b, cidx), -weight)
    # Literal epsilon Galileon; nonsymmetric affine Hessian, not a
    # Frobenius contraction. Keep both upper distortion indices.
    for mu, nu, rho, sigma in permutations(range(4)):
        for alpha, beta, gam in permutations(i for i in range(4) if i != sigma):
            weight = (
                f4
                * signs[sigma]
                * eps((mu, nu, rho, sigma))
                * eps((alpha, beta, gam, sigma))
                * v[mu]
                * v[alpha]
            )
            if weight == 0:
                continue
            for a in range(4):
                for b in range(4):
                    add(ix(a, beta, nu), ix(b, gam, rho), weight * v[a] * v[b])
    return sp.ImmutableMatrix(M.applyfunc(sp.factor))


@cache
def checks():
    p = sp.Rational(9, 20)
    e = old.quotient()["embedding"]
    gauge = old.quadratic()["gauge"]
    expected = -(2**52) * p**72 * (8 * p * p - 1) ** 3
    out = {}
    for name, v in (("timelike", (1, 0, 0, 0)), ("spacelike", (0, 1, 0, 0))):
        x = sum(old.SIGN[i] * v[i] ** 2 for i in range(4))
        c = 2 * (p - 2 * p * p) / x
        f4 = (sp.Rational(1, 2) - 2 * p * p) / (x * x)
        M = hessian(v, p, c, f4)
        out[name + "_projective_kernel"] = M * gauge
        Q = e.T * M * e
        out[name + "_complete_quotient_determinant"] = sp.factor(
            Q.det(method="domain-ge") - expected
        )
        if name == "timelike":
            out["general_gradient_timelike_reconstructs_full_literal_matrix"] = (
                M - old.quadratic()["hessian"].subs(old.P, p)
            )
    # A nonzero null gradient is not the constant-field vacuum. All
    # coefficient limits are inserted before evaluating the null matrix.
    p0 = sp.Rational(1, 2)
    M = hessian((1, 1, 0, 0), p0, 0, 512)
    out["nonzero_null_projective_kernel"] = M * gauge
    out["nonzero_null_full_quotient_determinant"] = sp.factor(
        (e.T * M * e).det(method="domain-ge") + 2**52 * p0**72 * (8 * p0 * p0 - 1) ** 3
    )
    M0 = hessian((0, 0, 0, 0), p0, 0, 512)
    out["zero_gradient_full_quotient_determinant"] = sp.factor(
        (e.T * M0 * e).det(method="domain-ge") + 2**52 * p0**72
    )
    return out
