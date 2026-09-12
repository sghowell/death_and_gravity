"""Literal ADM, finite covariance flow and independent scale-factor Euler tests."""

from functools import cache
from math import comb

import mpmath as mp
import sympy as s


def block(A, B):
    out = mp.zeros(A.rows + B.rows)
    out[: A.rows, : A.cols] = A
    out[A.rows :, A.cols :] = B
    return out


def positive_root(M):
    vals, U = mp.eigsy((M + M.T) / 2)
    return U * mp.diag([mp.sqrt(v) for v in vals]) * U.T


def literal_adm(phi, momentum):
    mass, a = mp.mpf(1000), mp.mpf(25) / 16
    h = a * a * mp.exp(2 * phi / 3) * mp.eye(3)
    volume = mp.sqrt(mp.det(h))
    k = mp.matrix([1, 2, 3]) * mp.mpf(momentum) / mp.sqrt(14)
    C = mp.matrix([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    K = (h + k * k.T / (mass * mass)) / volume
    V = mass * mass * volume * h**-1 + C.T * h * C / volume
    return block(V, K), K, V, k, h


def literal_current_vertex(phi, momentum):
    _M, K, _V, k, h = literal_adm(phi, momentum)
    omega = mp.sqrt(mp.mpf(1000) ** 2 + (k.T * h**-1 * k)[0])
    B = positive_root(K)
    T = block(B / mp.sqrt(omega), mp.sqrt(omega) * B**-1)
    MD = mp.diff(lambda e: literal_adm(e, momentum)[0], phi)
    return T.T * MD * T / omega


def projector_vertex(phi, momentum):
    A = mp.mpf(25) / 16 * mp.exp(phi / 3)
    z = mp.mpf(momentum) ** 2 / (A * A * mp.mpf(1000) ** 2 + mp.mpf(momentum) ** 2)
    n = mp.matrix([1, 2, 3]) / mp.sqrt(14)
    PL = n * n.T
    PT = mp.eye(3) - PL
    return block((1 - 2 * z) * PT / 3 + PL / 3, -PT / 3 - (1 + 2 * z) * PL / 3)


def finite_flow_fixture(momentum):
    _M0, K, _V, k, h = literal_adm(mp.mpf(0), momentum)
    omega = mp.sqrt(mp.mpf(1000) ** 2 + (k.T * h**-1 * k)[0])
    C = block(K / (2 * omega), omega * K**-1 / 2)
    J = block(mp.zeros(3), mp.zeros(3))
    J[:3, 3:] = mp.eye(3)
    J[3:, :3] = -mp.eye(3)
    metric = [
        mp.diff(lambda e: literal_adm(e, momentum)[0], mp.mpf(0), j) for j in range(4)
    ]
    generator = [J * metric[j] for j in range(3)]
    clock = mp.mpf(1) / 10000
    augmented = mp.zeros(18)
    for j in range(3):
        augmented[6 * j : 6 * (j + 1), 6 * j : 6 * (j + 1)] = generator[0]
    augmented[:6, 6:12] = generator[1]
    augmented[6:12, 12:] = generator[1]
    augmented[:6, 12:] = generator[2] / 2
    E = mp.expm(clock * augmented)
    U = [E[:6, :6], E[:6, 6:12], 2 * E[:6, 12:]]
    cov = [
        sum((comb(j, n) * U[n] * C * U[j - n].T for n in range(j + 1)), mp.zeros(6))
        for j in range(3)
    ]

    def trace(A):
        return mp.fsum(A[j, j] for j in range(A.rows))

    rebuilt = [
        -mp.fsum(comb(j, n) * trace(metric[n + 1] * cov[j - n]) for n in range(j + 1))
        / 2
        for j in range(3)
    ]

    def direct(e):
        M = literal_adm(e, momentum)[0]
        Ue = mp.expm(clock * J * M)
        MD = mp.diff(lambda f: literal_adm(f, momentum)[0], e)
        return -trace(MD * Ue * C * Ue.T) / 2

    actual = [mp.diff(direct, mp.mpf(0), j) for j in range(3)]
    contact = -trace(metric[2] * cov[0]) / 2
    return rebuilt, actual, contact


@cache
def literal_euler():
    A = s.symbols("A0:5", positive=True)
    mass = s.Symbol("mass", positive=True)
    rate = A[1] / A[0]
    accel = A[2] / A[0]
    R = 6 * (accel + rate**2)
    Ric = 9 * accel**2 + 3 * (accel + 2 * rate**2) ** 2
    density = s.expand(
        A[0] ** 3
        * (
            s.Rational(5, 2) * mass**4
            + s.Rational(5, 3) * mass**2 * R
            - R**2 / 30
            - Ric / 15
        )
    )
    dt = lambda f: s.expand(sum(s.diff(f, A[i]) * A[i + 1] for i in range(4)))
    Euler = (
        s.diff(density, A[0])
        - dt(s.diff(density, A[1]))
        + dt(dt(s.diff(density, A[2])))
    )
    current = s.factor(A[0] * Euler / 3)
    rates = [rate]
    for _ in range(3):
        rates.append(s.factor(dt(rates[-1])))
    return A, mass, density, current, tuple(rates)


@cache
def compact_hessian_fixture(which):
    A, mass, F, _current, _rates = literal_euler()
    t = s.Symbol("t", real=True)
    a = (1 + t * t) ** 2
    profiles = (
        ((t * t - s.Rational(1, 4)) ** 6, (1 + t) * (t * t - s.Rational(1, 4)) ** 5),
        (
            (1 + 2 * t) * (t * t - s.Rational(1, 4)) ** 6,
            (1 - t + t * t) * (t * t - s.Rational(1, 4)) ** 5,
        ),
    )
    G, D = profiles[which]
    base = {A[j]: s.diff(a, t, j) for j in range(5)}
    DG = [s.diff(a * G / 3, t, j) for j in range(3)]
    DD = [s.diff(a * D / 3, t, j) for j in range(3)]
    cross = [s.diff(a * G * D / 9, t, j) for j in range(3)]
    mixed = sum(s.diff(F, A[i]).subs(base) * cross[i] for i in range(3))
    mixed += sum(
        s.diff(F, A[i], A[j]).subs(base) * DG[i] * DD[j]
        for i in range(3)
        for j in range(3)
    )
    return t, a, G, D, s.factor(mixed.subs(mass, 1000))


def point_modes(e, t, momentum):
    Gamma = (t + mp.mpf(1) / 2) ** 3 / 6
    Gprime = (t + mp.mpf(1) / 2) ** 2 / 2
    a = (1 + t * t) ** 2 * mp.exp(e * Gamma / 3)
    H = 4 * t / (1 + t * t) + e * Gprime / 3
    omega = mp.sqrt(mp.mpf(1000) ** 2 + mp.mpf(momentum) ** 2 / (a * a))
    z = 1 - mp.mpf(1000) ** 2 / (omega * omega)
    return omega, (1 - z) * H / 2, (1 + z) * H / 2


@cache
def literal_dimensional_action():
    dimension, H, H1, q, x = s.symbols("dimension H H1 q x", real=True)
    scale = s.exp(H * x + H1 * x * x / 2)
    omega = s.sqrt(1 + q / scale**2)
    p = s.diff(omega, x) / omega
    rows = []
    for longitudinal in (False, True):
        kinetic = scale ** (2 - dimension) * (omega**2 if longitudinal else 1)
        squeeze = (p - s.diff(kinetic, x) / kinetic) / 2
        next_symbol = s.diff(squeeze, x) - p * squeeze
        rows.append((s.cancel(squeeze.subs(x, 0)), s.cancel(next_symbol.subs(x, 0))))
    (st, tt), (sl, tl) = rows

    def normalize(integrand, alpha):
        numerator, denominator = s.fraction(s.cancel(integrand))
        degree = s.degree(denominator, q)
        constant = denominator.subs(q, 0)
        assert s.expand(denominator - constant * (1 + q) ** degree) == 0
        out = 0
        for powers, co in s.Poly(numerator / constant, q).terms():
            n = powers[0]
            assert n <= degree
            out += (
                co
                * s.rf(dimension / 2, n)
                * s.rf(alpha - dimension / 2, degree - n)
                / s.rf(alpha, degree)
            )
        return s.factor(out)

    second = normalize(((dimension - 1) * st**2 + sl**2) / 4, s.Rational(1, 2))
    fourth = normalize(
        (dimension - 1) * (tt**2 + st**4) + tl**2 + sl**4, s.Rational(3, 2)
    )
    return dimension, H, H1, second, fourth
