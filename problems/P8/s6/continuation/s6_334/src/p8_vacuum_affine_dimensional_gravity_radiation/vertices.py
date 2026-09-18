"""Independent integer-D component expansion of the S304 canonical action."""

from itertools import permutations, product

import sympy as s


def component_engine(dimension):
    if type(dimension) is not int or dimension < 4:
        raise ValueError("Require an integer component dimension>=4")
    ETA = s.diag(1, *([-1] * (dimension - 1)))
    SGN = (1, *([-1] * (dimension - 1)))
    PARTS = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))

    def dot(p, q):
        return (p.T * ETA * q)[0]

    def pair(A, B):
        return sum(A[i, j] * B[i, j] for i, j in product(range(dimension), repeat=2))

    def stress(p, q, mu=1):
        return p * q.T + q * p.T - ETA * (dot(p, q) + mu)

    def field(p, q, mu=1):
        T = stress(p, q, mu)
        return (ETA * T * ETA - ETA * s.trace(ETA * T) / (dimension - 2)) / dot(
            p + q, p + q
        )

    def scalar4(A, B, p, q, mu=1):
        ta, tb = s.trace(ETA * A), s.trace(ETA * B)
        d2 = ta * tb / 4 - s.trace(ETA * A * ETA * B) / 2
        K = (
            ETA * d2
            - (ta * ETA * B * ETA + tb * ETA * A * ETA) / 2
            + ETA * (A * ETA * B + B * ETA * A) * ETA
        )
        return -4 * ((p.T * ETA * K * ETA * q)[0] + mu * d2)

    def connection(H, p):
        pl = ETA * p
        return {
            (r, m, n): s.Rational(SGN[r], 2)
            * (pl[m] * H[r, n] + pl[n] * H[r, m] - pl[r] * H[m, n])
            for r, m, n in product(range(dimension), repeat=3)
        }

    def second(A, G):
        return {
            (r, m, n): -SGN[r] * sum(A[r, z] * G[z, m, n] for z in range(dimension))
            for r, m, n in product(range(dimension), repeat=3)
        }

    def bilinear(M, G, H):
        return sum(
            M[m, n] * (G[r, m, n] * H[z, r, z] - G[r, m, z] * H[z, n, r])
            for m, n, r, z in product(range(dimension), repeat=4)
            if M[m, n]
        )

    def cubic(fields, momenta):
        gs = [connection(A, p) for A, p in zip(fields, momenta)]
        val = 0
        for a, b, c in permutations(range(3)):
            A = fields[a]
            M = s.trace(ETA * A) * ETA / 2 - ETA * A * ETA
            val += (
                bilinear(M, gs[b], gs[c])
                + bilinear(ETA, second(A, gs[b]), gs[c])
                + bilinear(ETA, gs[b], second(A, gs[c]))
            )
        return -4 * val

    def channel(ps, k, eps, left, right, mu=1):
        external = 0
        for i, j, other in (
            (left[0], left[1], right),
            (left[1], left[0], right),
            (right[0], right[1], left),
            (right[1], right[0], left),
        ):
            p, q = ps[i], ps[j]
            external += (
                pair(eps, stress(p, -p - k, mu))
                * pair(field(ps[other[0]], ps[other[1]], mu), stress(p + k, q, mu))
                / (2 * dot(p, k))
            )
        HL, HR = (
            field(ps[left[0]], ps[left[1]], mu),
            field(ps[right[0]], ps[right[1]], mu),
        )
        sg = -scalar4(eps, HR, ps[left[0]], ps[left[1]], mu) - scalar4(
            eps, HL, ps[right[0]], ps[right[1]], mu
        )
        tri = cubic(
            (HL, HR, eps), (ps[left[0]] + ps[left[1]], ps[right[0]] + ps[right[1]], k)
        )
        return tuple(s.factor(v) for v in (external, sg, tri))

    def amplitude(ps, k, eps, mu=1):
        return sum(sum(channel(ps, k, eps, *part, mu)) for part in PARTS)

    def born(ps, mu=1):
        return -sum(
            pair(stress(ps[a], ps[b], mu), field(ps[c], ps[d], mu))
            for (a, b), (c, d) in PARTS
        )

    return {
        "ETA": ETA,
        "dot": dot,
        "stress": stress,
        "field": field,
        "scalar4": scalar4,
        "cubic": cubic,
        "channel": channel,
        "amplitude": amplitude,
        "born": born,
    }
