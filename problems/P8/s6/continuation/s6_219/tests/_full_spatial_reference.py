"""Independent literal metric-Hamiltonian and exact original-band fixtures."""

from functools import cache

import mpmath as mp
import sympy as s
from p8_vacuum_affine_full_spatial_remainder import shapes


def literal_pair_hamiltonian(Q, k, ell, left, right, a, mass):
    h = a * a * mp.expm(Q)
    inv = h**-1
    volume = mp.sqrt(mp.det(h))
    A, pi = left[:3, :], left[3:, :]
    B, rho = right[:3, :], right[3:, :]
    F = mp.matrix(3, 3)
    E = mp.matrix(3, 3)
    for i in range(3):
        for j in range(3):
            F[i, j] = 1j * (k[i] * A[j] - k[j] * A[i])
            E[i, j] = 1j * (ell[i] * B[j] - ell[j] * B[i])
    magnetic = mp.fsum(
        inv[i, q] * inv[j, r] * F[i, j] * E[q, r]
        for i in range(3)
        for j in range(3)
        for q in range(3)
        for r in range(3)
    )
    return (
        (pi.T * h * rho)[0] / volume
        + volume * mass * mass * (A.T * inv * B)[0]
        + volume * magnetic / 2
        - (k.T * pi)[0] * (ell.T * rho)[0] / (volume * mass * mass)
    )


def symplectic(n):
    J = mp.zeros(2 * n)
    for i in range(n):
        J[i, n + i], J[n + i, i] = 1, -1
    return J


def vacuum_columns(n):
    out = mp.zeros(2 * n, n)
    for i in range(n):
        out[i, i], out[n + i, i] = 1 / mp.sqrt(2), -1j / mp.sqrt(2)
    return out


def pair_covariance(U, V, D, G):
    n = U.rows // 2
    J, v = symplectic(n), vacuum_columns(n)
    source = (J * G + G * J.T) / 2
    covariance = -mp.fsum((D.T * U * source * V.T)[i, i] for i in range(2 * n))
    bd, bg = v.T * U.T * D * V * v, v.T * G * v
    kubo = 2 * mp.im(
        mp.fsum(mp.conj(bd[i, j]) * bg[i, j] for i in range(n) for j in range(n))
    )
    wrong = -2 * mp.im(mp.fsum(bd[i, j] * bg[i, j] for i in range(n) for j in range(n)))
    return covariance, kubo, wrong


@cache
def raw_function(channel):
    rows = shapes.raw_grades()
    keys = [(q, r) for q in range(5) for r in range(q + 1)]
    fn = s.lambdify(
        (shapes.t, shapes.p, shapes.m, shapes.u),
        [rows[channel, q, r] for q, r in keys],
        "mpmath",
        cse=True,
    )
    return keys, fn


def original_band_difference(
    channel, time, transfer, mass, K, source_jets, omit_grazing=False
):
    """Literal ray intersections; no asymptotic shell/shape routine."""
    if not K > transfer + mass:
        raise ValueError("This fixture uses a lost shell above the fixed lower band")
    keys, fn = raw_function(channel)

    def angular(u):
        upper = transfer * u + mp.sqrt(K * K - transfer * transfer * (1 - u * u))
        values = fn(time, transfer, mass, u)
        total = 0
        for (q, r), value in zip(keys, values):
            moment = (
                (K ** (4 - q) - upper ** (4 - q)) / (4 - q)
                if q < 4
                else mp.log(K / upper)
            )
            total -= source_jets[r] * value * moment / (4 * mp.pi**2)
        return total

    endpoint = mp.mpf(0) if omit_grazing else transfer / (2 * K)
    return (
        mp.quad(angular, [-1, 0, endpoint]) if endpoint else mp.quad(angular, [-1, 0])
    )


def expected_nondecaying(channel, time, transfer, mass, K, source_jets):
    cubic, linear = shapes.predicted_shapes(channel)
    substitution = {
        shapes.t: s.Rational(str(time)),
        shapes.p: s.Rational(str(transfer)),
        shapes.m: s.Rational(str(mass)),
    }
    evaluate = lambda value: mp.mpf(str(s.N(value.subs(substitution), mp.mp.dps)))
    return K**3 * evaluate(cubic) * source_jets[0] + K * mp.fsum(
        evaluate(value) * source_jets[r] for r, value in enumerate(linear)
    )
