"""Independent whole-tree, recoil, phase-space and analytic-majorant checks."""

import math

import numpy as np
import pytest
import sympy as s
from numpy.polynomial.legendre import leggauss
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import (
    audit,
    rate,
    recoil,
    source,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_identity(name):
    value = ROWS[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).atoms(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_written_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_every_unsupported_scope_input(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "name", tuple(recoil.data()["explicit_nonnegative_majorant_margins"])
)
def test_each_exact_majorant_margin(name):
    value = s.sympify(recoil.data()["explicit_nonnegative_majorant_margins"][name])
    num, den = s.fraction(s.factor(value))
    variables = tuple(sorted(value.free_symbols, key=str)) or (s.Symbol("unused"),)
    assert all(c >= 0 for c in s.Poly(num, *variables).coeffs())
    assert all(c >= 0 for c in s.Poly(den, *variables).coeffs()) and den != 0
    assert not value.atoms(s.Float)


eta = np.diag([1.0, -1.0, -1.0, -1.0])


def dot(a, b):
    return a @ eta @ b


def outer(a, b):
    return np.outer(a, b)


def stress(p, r, m):
    return outer(p, r) + outer(r, p) - eta * (dot(p, r) - m)


def momenta(E, w, nhat, u):
    ep = math.sqrt(E * (E - w))
    r = math.sqrt(E * (E - w) - 1)
    r0 = math.sqrt(E * E - 1)
    c = u @ nhat
    gam = (E - w / 2) / ep
    out = []
    for sign in (1.0, -1.0):
        out.append(
            np.r_[
                E - w / 2 - sign * w * r * c / (2 * ep),
                sign * r * u + (sign * (gam - 1) * r * c - w / 2) * nhat,
            ]
        )
    k = [np.array([-E, 0.0, 0.0, -r0]), np.array([-E, 0.0, 0.0, r0]), *out]
    k0 = [k[0], k[1], np.r_[E, r0 * u], np.r_[E, -r0 * u]]
    return k, k0, np.r_[w, w * nhat], r * E / (ep * r0)


def positive_born(k, n):
    a = [dot(k[i] + k[j], k[i] + k[j]) for i in range(4) for j in range(i + 1, 4)]
    return sum((x - 2) ** 2 / (n - x) for x in a) / (2 * (n - 2) ** 2)


def sew(a, b, nhat):
    pi = np.eye(3) - outer(nhat, nhat)
    return np.trace(a @ pi @ b @ pi) - np.trace(a @ pi) * np.trace(b @ pi) / 2


def amplitudes(E, w, nhat, u, n, check=False):
    k, k0, q, jac = momenta(E, w, nhat, u)
    ji = [outer(p, p) / dot(p, q) for p in k]
    j0 = [outer(p, p) / dot(p, q) for p in k0]
    abar = positive_born(k, n)
    a0 = positive_born(k0, n)
    rem = np.zeros((4, 4))
    literal = (-3 / (n - 2) + 2 / (n - 2) ** 2) * (
        sum([stress(p, p + q, 1) / (2 * dot(p, q)) for p in k]) + eta
    )
    ti = [stress(p, p + q, 1) / (2 * dot(p, q)) for p in k]
    for j in (1, 2, 3):
        L = (0, j)
        R = tuple(i for i in range(4) if i not in L)
        PL = k[0] + k[j]
        PR = sum(k[i] for i in R)
        a = dot(PL, PL)
        b = dot(PR, PR)
        dl = a - n
        dr = b - n
        bb = dot(PL, q) * (sum(ji[i] for i in L) - sum(ji[i] for i in R)) - 2 * outer(
            PL, PL
        )
        # Stable exact denominator subtraction, no subtraction of two close floats.
        rem += bb * (n * (a + b) - a * b) / (n * n * dl * dr)
        literal -= (
            sum(ti[i] for i in L) / dr
            + sum(ti[i] for i in R) / dl
            + stress(PL, PL + q, n) / (dl * dr)
            + eta * (1 / dl + 1 / dr)
        )
    lead = a0 * sum(j0)
    delta = (abar - a0) * sum(ji) + a0 * (sum(ji) - sum(j0)) + rem
    if check:
        assert np.max(np.abs(sum(k) + q)) < 1e-13
        assert max(abs(dot(p, p) - 1) for p in k) < 1e-12
        assert np.max(np.abs(literal @ eta @ q)) < 1e-11
        err = literal - lead - delta
        # The improved remainder differs by pure gauge, tested by complete TT sew.
        assert abs(sew(err[1:, 1:], err[1:, 1:], nhat)) < 1e-18
        assert (
            math.sqrt(max(0.0, sew(delta[1:, 1:], delta[1:, 1:], nhat))) / a0
            < math.sqrt(2) * 330000
        )
        assert 0 <= 1 - jac <= 2 * w
    l = lead[1:, 1:]
    d = delta[1:, 1:]
    density = jac * (2 * sew(l, d, nhat) + sew(d, d, nhat)) + (jac - 1) * sew(
        l, l, nhat
    )
    return density / (a0 * a0)


@pytest.mark.parametrize(
    "E,n,c", ((1.25, 128.0, 0.6), (1.5, 257.0, -0.3), (2.0, 1024.0, 0.2))
)
@pytest.mark.parametrize("w", (0.125, 0.01, 0.001))
def test_independent_literal_tree_recoil_Ward_TT_and_bound(E, n, c, w):
    u = np.array([math.sqrt(1 - c * c), 0.0, c])
    for nhat in (
        np.array([0.0, 0.0, 1.0]),
        np.array([0.6, 0.0, 0.8]),
        np.array([0.0, 1.0, 0.0]),
    ):
        amplitudes(E, w, nhat, u, n, True)


def integrate_real_difference(E, n, c, cut, order):
    u = np.array([math.sqrt(1 - c * c), 0.0, c])
    xx, ww = leggauss(order)
    angular = []
    for z, wz in zip(xx, ww):
        for ph in np.arange(2 * order) * math.pi / order:
            nhat = np.array(
                [
                    math.sqrt(1 - z * z) * math.cos(ph),
                    math.sqrt(1 - z * z) * math.sin(ph),
                    z,
                ]
            )
            angular.append((nhat, wz * math.pi / order))
    total = 0.0
    for x, wx in zip(xx, ww):
        omega = cut * (x + 1) / 2
        avg = sum(weight * amplitudes(E, omega, nhat, u, n) for nhat, weight in angular)
        total += wx * cut / 2 * omega * avg / (2 * (2 * math.pi) ** 3)
    return total


@pytest.mark.parametrize(
    "E,n,c", ((1.25, 128.0, 0.6), (1.5, 257.0, -0.3), (2.0, 1024.0, 0.2))
)
@pytest.mark.parametrize("cut", (0.125, 0.0625))
def test_independent_complete_real_rate_quadrature_calibration(E, n, c, cut):
    coarse = integrate_real_difference(E, n, c, cut, 8)
    fine = integrate_real_difference(E, n, c, cut, 16)
    assert abs(fine - coarse) < 5e-6
    assert abs(fine) < float(rate.finite_real_rate_error(s.Rational(str(cut)), 1))
    # Calibration only: the written analytic majorant proves the domain bound.


@pytest.mark.parametrize(
    "direction,outgoing",
    (
        ([1, 0], [1, 0, 0]),
        ([1, 0, 0], [1, 0]),
        ([2, 0, 0], [1, 0, 0]),
        ([1, 0, 0], [1, 1, 0]),
    ),
)
def test_recoil_rejects_nonunit_or_wrong_dimension(direction, outgoing):
    with pytest.raises(ValueError, match="Require"):
        recoil.momenta(s.Rational(3, 2), s.Rational(1, 16), direction, outgoing)


def test_original_error_bound_exact_without_large_float_conversion():
    numerator = s.Integer(42248192) / 8 + s.Integer(54450000000) / 64
    assert numerator / 18 < 10**8
    assert s.Integer(10) ** 8 / source.KAPPA == s.Rational(1, 10**792)
    assert source.HEAVY_MASS2 >= 128
