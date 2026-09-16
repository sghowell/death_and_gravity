"""Independent whole-D soft, phase-space, virtual and continuity tests."""

import math

import mpmath as mp
import numpy as np
import pytest
import sympy as s
from numpy.polynomial.legendre import leggauss
from p8_vacuum_affine_physical_virtual_soft_pairing import (
    audit,
    continuity,
    inclusive,
    soft,
    source,
)
from scipy.special import gamma

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_each_exact_identity(name):
    value = ROWS[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).atoms(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_each_written_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_each_scope_guard(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "name", tuple(continuity.data()["explicit_nonnegative_margins"])
)
def test_every_majorant_margin(name):
    value = s.sympify(continuity.data()["explicit_nonnegative_margins"][name])
    num, den = s.fraction(s.factor(value))
    variables = tuple(sorted(value.free_symbols, key=str)) or (s.Symbol("unused"),)
    assert all(c >= 0 for c in s.Poly(num, *variables).coeffs())
    assert all(c >= 0 for c in s.Poly(den, *variables).coeffs()) and den != 0


def ratioatanh(b):
    return mp.mpf(1) if b == 0 else mp.atanh(mp.sqrt(b)) / mp.sqrt(b)


def J(e, b, E):
    return mp.hyp2f1(1, mp.mpf("1.5"), mp.mpf("1.5") + e, b) / (E * E)


def kernels(E, c, ep):
    mu = mp.mpf(1)
    ss = 4 * E * E
    tt = -2 * (E * E - 1) * (1 - c)
    uu = -2 * (E * E - 1) * (1 + c)
    ds = (ss - 2) / 2
    dt = 1 - tt / 2
    du = 1 - uu / 2
    pair = lambda d: mp.quad(
        lambda x: J(ep, 1 - (1 + 2 * x * (1 - x) * (d - 1)) / (E * E), E), [0, 1]
    )
    N = lambda d: d * d - 1 / (2 + 2 * ep)
    return 4 * N(mu) * J(ep, 1 - 1 / (E * E), E) + 4 * (
        N(ds) * pair(ds) - N(dt) * pair(dt) - N(du) * pair(du)
    )


def derivative_kernel(E, c):
    ss = 4 * E * E
    tt = -2 * (E * E - 1) * (1 - c)
    uu = -2 * (E * E - 1) * (1 + c)

    def pair(d):
        def f(x):
            mass = 1 + 2 * x * (1 - x) * (d - 1)
            b = 1 - mass / (E * E)
            return (
                mp.mpf(".5") + 2 * (d * d - mp.mpf(".5")) * (1 - ratioatanh(b))
            ) / mass

        return mp.quad(f, [0, 1])

    return (
        6
        - 4 * ratioatanh(1 - 1 / (E * E))
        + 4 * (pair((ss - 2) / 2) - pair(1 - tt / 2) - pair(1 - uu / 2))
    )


@pytest.mark.parametrize("speed", ("0", "0.2", "0.6", "0.8660254037844386"))
def test_independent_hypergeometric_and_log_angular_derivative(speed):
    with mp.workdps(50):
        beta = mp.mpf(speed)
        b = beta * beta
        derivative = mp.diff(lambda ep: J(ep, b, mp.mpf(2)), 0)
        exact = 2 * (1 - ratioatanh(b)) / (4 * (1 - b))
        assert abs(derivative - exact) < mp.mpf("1e-40")
        raw = mp.quad(lambda z: mp.log(1 - z * z) / (1 - beta * z) ** 2, [-1, 0, 1]) / 2
        expected = (2 * mp.log(2) - 2 * ratioatanh(b)) / (1 - b)
        assert abs(raw - expected) < mp.mpf("1e-40")


@pytest.mark.parametrize("hard,cos", (("1.25", "0.6"), ("1.5", "-0.3"), ("2", "0.2")))
@pytest.mark.parametrize("resolution", ("0.125", "0.015625"))
def test_full_D_real_virtual_finite_conversion_and_regulator_bound(
    hard, cos, resolution
):
    with mp.workdps(50):
        E, c, res = mp.mpf(hard), mp.mpf(cos), mp.mpf(resolution)
        k0 = kernels(E, c, 0)
        k1 = derivative_kernel(E, c)
        assert k0 >= 0
        assert abs(mp.diff(lambda ep: kernels(E, c, ep), 0) - k1) < mp.mpf("1e-38")
        delta = (k1 + (mp.euler - 2 - mp.log(mp.pi)) * k0) / (8 * mp.pi**2)
        assert abs(delta) < 112
        errors = []
        for ep in (mp.mpf(".001"), mp.mpf(".0001")):
            phase = (
                (4 * mp.pi) ** (-ep)
                * mp.gamma(mp.mpf("1.5"))
                / mp.gamma(mp.mpf("1.5") + ep)
            )
            paired = (
                res ** (2 * ep) * (phase * kernels(E, c, ep) - k0) / (8 * mp.pi**2 * ep)
            )
            error = abs(paired - delta)
            bound = ep * (16016 * abs(mp.log(res)) + 301824) / (8 * mp.pi**2)
            assert error < bound
            errors.append(error)
        assert errors[1] < errors[0] / 5


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


def sew(a, b, nhat, e=0.0):
    pi = np.eye(3) - outer(nhat, nhat)
    return np.trace(a @ pi @ b @ pi) - np.trace(a @ pi) * np.trace(b @ pi) / (2 + 2 * e)


def amplitudes(E, w, nhat, u, n, check=False, e=0.0):
    k, k0, q, jac = momenta(E, w, nhat, u)
    jac *= ((E * (E - w) - 1) / (E * E - 1)) ** e
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
    density = jac * (2 * sew(l, d, nhat, e) + sew(d, d, nhat, e)) + (jac - 1) * sew(
        l, l, nhat, e
    )
    return density / (a0 * a0)


def integral(E, n, c, cut, ep, order):
    u = np.array([math.sqrt(1 - c * c), 0.0, c])
    xx, ww = leggauss(order)
    angular = []
    for x, weight in zip(xx, ww):
        v = (x + 1) / 2
        r = math.sqrt(1 - v * v)
        for ph in np.arange(2 * order) * math.pi / order:
            nhat = np.array([r * math.cos(ph), v, r * math.sin(ph)])
            angular.append(
                (nhat, weight / 2 * (1 + 2 * ep) * v ** (2 * ep) / (2 * order))
            )
    phase = (4 * math.pi) ** (-ep) * gamma(1.5) / gamma(1.5 + ep)
    total = 0.0
    for x, weight in zip(xx, ww):
        omega = cut * (x + 1) / 2
        avg = sum(w * amplitudes(E, omega, nhat, u, n, e=ep) for nhat, w in angular)
        total += (
            weight * cut / 2 * omega ** (1 + 2 * ep) * avg * phase / (4 * math.pi**2)
        )
    return total


@pytest.mark.parametrize(
    "E,n,c,reference",
    (
        (1.25, 128.0, 0.6, -0.006790638279097058),
        (1.5, 257.0, -0.3, -0.022734235934448472),
        (2.0, 1024.0, 0.2, -0.053423482624840715),
    ),
)
def test_independent_entire_D_recoil_rate_limit(E, n, c, reference):
    baseline = integral(E, n, c, 0.125, 0.0, 16)
    assert abs(baseline - reference) < 3e-6
    errors = []
    for ep in (0.001, 0.0001):
        value = integral(E, n, c, 0.125, ep, 16)
        errors.append(abs(value - baseline))
    assert errors[1] < errors[0] / 5
    # Independent calibration; the positive fixed-domain majorant proves the limit.


@pytest.mark.parametrize("ep", (0.0, 0.001, 0.1))
def test_independent_folded_sphere_normalization(ep):
    x, w = leggauss(96)
    v = (x + 1) / 2
    norm = np.dot(w / 2, (1 + 2 * ep) * v ** (2 * ep))
    assert abs(norm - 1) < 1e-5


@pytest.mark.parametrize("direction", (1, -1))
def test_complete_soft_endpoints_not_a_dropped_phase(direction):
    assert s.factor(soft.kernel(soft.E, direction, source.EP, source.MU)) == 0
    assert s.factor(soft.kernel_first(soft.E, direction, source.MU)) == 0
    assert "Coulomb" in inclusive.data()["whole_physical_sheet_link"]


def test_original_total_reference_error_without_inexact_large_numbers():
    numerator = s.Integer(42248192) / 8 + s.Integer(54450000000) / 64
    assert numerator / 18 + 112 < 10**8
    assert s.Integer(10) ** 8 / source.KAPPA == s.Rational(1, 10**792)
    assert source.HEAVY_MASS2 >= 128
