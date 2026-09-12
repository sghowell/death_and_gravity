"""Independent exact-frequency, time-dependent-kinetic two-mode fixture."""

import mpmath as mp


def dagger(M):
    return M.transpose_conj()


def kron(A, B):
    return mp.matrix(
        [
            [
                A[i // B.rows, j // B.cols] * B[i % B.rows, j % B.cols]
                for j in range(A.cols * B.cols)
            ]
            for i in range(A.rows * B.rows)
        ]
    )


def fundamental(t, w, b):
    c, s = mp.cos(w * t), mp.sin(w * t)
    return mp.matrix(
        [
            [mp.exp(b * t) * c, mp.exp(b * t) * s],
            [mp.exp(-b * t) * (b * c - w * s), mp.exp(-b * t) * (b * s + w * c)],
        ]
    ) / mp.sqrt(w)


def mode(t, w, b):
    f = mp.exp(b * t - 1j * w * t) / mp.sqrt(2 * w)
    return mp.matrix([f, (b - 1j * w) * mp.exp(-2 * b * t) * f])


def covariance(t, w, b):
    u = mode(t, w, b)
    return mp.matrix(
        [[mp.re(u[i] * mp.conj(u[j])) for j in range(2)] for i in range(2)]
    )


def strip_terms(M, w1, b1, w2, b2):
    norm = 1 / (2 * mp.sqrt(w1 * w2))
    return [
        (
            M[i, j]
            * norm
            * ((b1 - 1j * w1) if i else 1)
            * ((b2 - 1j * w2) if j else 1),
            (1 - 2 * i) * b1 + (1 - 2 * j) * b2,
        )
        for i in range(2)
        for j in range(2)
    ]


def stripped(M, t, w1, b1, w2, b2):
    return mp.fsum(c * mp.exp(rate * t) for c, rate in strip_terms(M, w1, b1, w2, b2))


def moment(power, z, T):
    # Exact primitive of s^power exp(z*s) from0 toT; z has nonzero imaginary part.
    return mp.exp(z * T) * mp.fsum(
        (-1) ** j
        * mp.factorial(power)
        / mp.factorial(power - j)
        * T ** (power - j)
        / z ** (j + 1)
        for j in range(power + 1)
    ) - (-1) ** power * mp.factorial(power) / z ** (power + 1)


def source_derivative(M, t, order, w1, b1, w2, b2):
    return mp.fsum(
        c
        * mp.exp(rate * t)
        * mp.fsum(
            mp.binomial(order, h)
            * mp.factorial(6)
            / mp.factorial(6 - h)
            * t ** (6 - h)
            * rate ** (order - h)
            for h in range(min(order, 6) + 1)
        )
        for c, rate in strip_terms(M, w1, b1, w2, b2)
    )


def pair_integral(D, G, T, w1, b1, w2, b2, wrong=False):
    Omega = w1 + w2
    if wrong:
        source = mp.fsum(
            c * moment(6, rate + 1j * Omega, T)
            for c, rate in strip_terms(G, w1, b1, w2, b2)
        )
        return -2 * mp.im(
            mp.conj(stripped(D, T, w1, b1, w2, b2)) * mp.exp(-1j * Omega * T) * source
        )
    source = mp.fsum(
        c * moment(6, rate - 1j * Omega, T)
        for c, rate in strip_terms(G, w1, b1, w2, b2)
    )
    return 2 * mp.im(
        mp.conj(stripped(D, T, w1, b1, w2, b2)) * mp.exp(1j * Omega * T) * source
    )


def endpoint_and_bulk(D, G, T, w1, b1, w2, b2):
    Omega = w1 + w2
    detector = mp.conj(stripped(D, T, w1, b1, w2, b2))
    endpoints = [
        2
        * mp.im(
            detector
            * 1j
            * (-1j) ** j
            * source_derivative(G, T, j, w1, b1, w2, b2)
            / Omega ** (j + 1)
        )
        for j in range(6)
    ]
    source = mp.fsum(
        c
        * mp.fsum(
            mp.binomial(6, h)
            * mp.factorial(6)
            / mp.factorial(6 - h)
            * rate ** (6 - h)
            * moment(6 - h, rate - 1j * Omega, T)
            for h in range(7)
        )
        for c, rate in strip_terms(G, w1, b1, w2, b2)
    )
    bulk = 2 * mp.im(detector * (-1j) ** 6 * mp.exp(1j * Omega * T) * source / Omega**6)
    return endpoints, bulk


def covariance_integrand(D, G, T, t, w1, b1, w2, b2):
    U = fundamental(T, w1, b1) * fundamental(t, w1, b1) ** -1
    V = fundamental(T, w2, b2) * fundamental(t, w2, b2) ** -1
    J = mp.matrix([[0, 1], [-1, 0]])
    source = J * G * covariance(t, w2, b2) + covariance(t, w1, b1) * G * J.T
    matrix = D.T * U * source * V.T
    return -sum(matrix[i, i] for i in range(2)) * t**6


def fock_integrand(D, G, T, t, w1, b1, w2, b2):
    ann = mp.zeros(3)
    ann[0, 1] = 1
    ann[1, 2] = mp.sqrt(2)
    aa = (kron(ann, mp.eye(3)), kron(mp.eye(3), ann))

    def obs(M, t):
        u = (mode(t, w1, b1), mode(t, w2, b2))
        fields = [
            [u[i][j] * aa[i] + mp.conj(u[i][j]) * dagger(aa[i]) for j in range(2)]
            for i in range(2)
        ]
        H = mp.zeros(9)
        for i in range(2):
            for j in range(2):
                H += M[i, j] * fields[0][i] * fields[1][j]
        return H

    hd, hg = obs(D, T), obs(G, t)
    return mp.re(1j * (hd * hg - hg * hd)[0, 0]) * t**6
