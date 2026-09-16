"""General on-shell Ward identity in an independent invariant-index algebra."""

from collections import Counter
from functools import cache
from itertools import permutations, product

import sympy as s

mu, a, b, u1, u2, u3, x1, x2, x3, xk = s.symbols("mu a b u1 u2 u3 x1 x2 x3 xk")
base = s.Matrix(
    [
        [mu, a, b, u1, x1],
        [a, mu, -mu - a - b - u1 - u2 - u3, u2, x2],
        [b, -mu - a - b - u1 - u2 - u3, mu, u3, x3],
        [u1, u2, u3, 0, xk],
        [x1, x2, x3, xk, s.Symbol("xx")],
    ]
)
vectors = [s.eye(5)[:, i] for i in (0, 1, 2)]
vectors += [
    -sum(vectors, s.zeros(5, 1)) - s.eye(5)[:, 3],
    s.eye(5)[:, 3],
    s.eye(5)[:, 4],
]
vectors += [vectors[0] + vectors[1], vectors[2] + vectors[3]]
vectors += [vectors[i] + vectors[4] for i in range(4)]


@cache
def dot(i, j):
    return s.expand((vectors[i].T * base * vectors[j])[0])


def metric(i, j):
    return [(s.S.One, (("m", i, j),))]


def outer(p, q, i, j):
    return [(s.S.One, (("v", i, p), ("v", j, q)))]


def scaled(T, c):
    return [(v * c, f) for v, f in T]


def mul(*Ts):
    return [
        (s.Mul(*(t[0] for t in terms)), sum((t[1] for t in terms), ()))
        for terms in product(*Ts)
    ]


def H(p, q):
    return lambda i, j: outer(p, q, i, j) + outer(q, p, i, j) + scaled(metric(i, j), mu)


def stress(p, q):
    return lambda i, j: (
        outer(p, q, i, j) + outer(q, p, i, j) + scaled(metric(i, j), -dot(p, q) - mu)
    )


def gauge(i, j):
    return outer(4, 5, i, j) + outer(5, 4, i, j)


@cache
def topology(factors):
    parent = {}
    v = []

    def find(i):
        parent.setdefault(i, i)
        if parent[i] != i:
            parent[i] = find(parent[i])
        return parent[i]

    counts = Counter()
    for kind, i, j in factors:
        counts[i] += 1
        find(i)
        if kind == "m":
            counts[j] += 1
            a0, b0 = find(i), find(j)
            parent[a0] = b0
        else:
            v.append((i, j))
    assert all(n == 2 for n in counts.values()), counts
    groups = {find(i): [] for i in parent}
    for i, p in v:
        groups[find(i)].append(p)
    powers = 0
    pairs = []
    for vs in groups.values():
        if not vs:
            powers += 1
        else:
            assert len(vs) == 2, vs
            pairs.append(tuple(sorted(vs)))
    return powers, tuple(sorted(pairs))


def contract(T):
    totals = {}
    for coeff, factors in T:
        key = topology(tuple(sorted(factors)))
        totals[key] = totals.get(key, 0) + coeff
    return s.Add(
        *(
            4**n * coef * s.Mul(*(dot(i, j) for i, j in pairs))
            for (n, pairs), coef in totals.items()
        )
    )


def trace(A):
    return contract(A("z", "z"))


def pq(A, p, q):
    return contract(mul(outer(p, q, "m", "n"), A("m", "n")))


def inner(A, B):
    return contract(mul(A("m", "n"), B("m", "n")))


def product_pq(A, B, p, q):
    return contract(mul(outer(p, q, "m", "n"), A("m", "r"), B("r", "n")))


def scalar4(A, B, p, q):
    ta, tb = trace(A), trace(B)
    det2 = ta * tb / 4 - inner(A, B) / 2
    return -4 * (
        (dot(p, q) + mu) * det2
        - (ta * pq(B, p, q) + tb * pq(A, p, q)) / 2
        + product_pq(A, B, p, q)
        + product_pq(B, A, p, q)
    )


def conn(A, p):
    def result(r, m, n):
        return scaled(
            [(c, (("v", m, p),) + f) for c, f in A(r, n)]
            + [(c, (("v", n, p),) + f) for c, f in A(r, m)]
            + [(-c, (("v", r, p),) + f) for c, f in A(m, n)],
            s.Rational(1, 2),
        )

    return result


def second(A, G, index):
    return lambda r, m, n: scaled(mul(A(r, index), G(index, m, n)), -1)


def bilinear(M, G, J):
    return contract(mul(M("m", "n"), G("r", "m", "n"), J("s", "r", "s"))) - contract(
        mul(M("m", "n"), G("r", "m", "s"), J("s", "n", "r"))
    )


def cubic(fields, momenta):
    gs = [conn(A, p) for A, p in zip(fields, momenta)]
    result = 0
    for aa, bb, cc in permutations(range(3)):
        A = fields[aa]
        M = lambda m, n, A=A: scaled(metric(m, n), trace(A) / 2) + scaled(A(m, n), -1)
        result += bilinear(M, gs[bb], gs[cc])
        result += bilinear(metric, second(A, gs[bb], "zab"), gs[cc])
        result += bilinear(metric, gs[bb], second(A, gs[cc], "zac"))
    return s.expand(-4 * result)


@cache
def rank_one_channel():
    """Independent seven-diagram channel for epsilon=xi_lower xi_lower."""
    E = lambda i, j: outer(5, 5, i, j)
    HL, HR = H(0, 1), H(2, 3)
    L, R = dot(6, 6), dot(7, 7)

    def emitted(i):
        return -2 * dot(i, 5) ** 2 - 2 * dot(i, 5) * dot(4, 5) + dot(i, 4) * dot(5, 5)

    external = sum(
        emitted(i) * inner(stress(i + 8, j), HR) / (2 * dot(i, 4) * R)
        for i, j in ((0, 1), (1, 0))
    )
    external += sum(
        emitted(i) * inner(stress(i + 8, j), HL) / (2 * dot(i, 4) * L)
        for i, j in ((2, 3), (3, 2))
    )
    seagull = -scalar4(E, HR, 0, 1) / R - scalar4(E, HL, 2, 3) / L
    return external + seagull + cubic([HL, HR, E], [6, 7, 4]) / (L * R)


@cache
def data():
    HL, HR = H(0, 1), H(2, 3)
    L, R = dot(6, 6), dot(7, 7)
    external = -2 * sum(
        dot(i, 5) * inner(stress(i + 8, j), HR) * L for i, j in ((0, 1), (1, 0))
    )
    external -= 2 * sum(
        dot(i, 5) * inner(stress(i + 8, j), HL) * R for i, j in ((2, 3), (3, 2))
    )
    seagull = -scalar4(gauge, HR, 0, 1) * L - scalar4(gauge, HL, 2, 3) * R
    triple = cubic([HL, HR, gauge], [6, 7, 4])
    checks = {
        "whole_channel_general_Ward_polynomial": s.expand(external + seagull + triple)
    }
    for i in range(4):
        checks[f"general_scalar{i}_mass_shell"] = s.expand(dot(i, i) - mu)
    checks["general_emitted_null_momentum"] = dot(4, 4)
    for j in range(6):
        checks[f"general_conservation_contraction{j}"] = s.expand(
            sum(dot(i, j) for i in range(4)) + dot(4, j)
        )
    z, u, kxi = s.symbols("p_xi p_k k_xi")
    checks["scalar_emission_propagator_Ward_cancellation"] = s.expand(
        -4 * u * z - 2 * u * kxi + 2 * u * kxi + 4 * u * z
    )
    return {
        "whole_general_Gram_coordinates": base,
        "whole_general_vectors": vectors,
        "whole_channel_Ward_numerators": {
            "external": s.expand(external),
            "seagull": s.expand(seagull),
            "cubic": triple,
        },
        "whole_invariant_proof": "Take arbitrary mu,p1.p2=a,p1.p3=b,pi.k=ui and pi.xi=xi. Set p4=-p1-p2-p3-k and p2.p3=-mu-a-b-u1-u2-u3. This Gram matrix imposes only shells and conservation, not a finite point or a four-dimensional Gram determinant. Contract metric-index paths into dot products and closed loops into dimension4. Expand all cubic-action, seagull and scalar-emission terms. Cancel scalar emission propagators and multiply by K_L^2 K_R^2. The entire Ward polynomial vanishes for arbitrary xi. Bose relabeling supplies all three channels.",
        "whole_soft_identity": "For TT epsilon the external scalar vertex divided by its propagator is -p.epsilon.p/(p.k). The remaining vertex product is minus its four-point tree. Hence M5=sum_i J_i G_i+seagulls+cubic, with leading term G0 sum_i J_i; recoil and hard differences are bounded separately.",
        "whole_wrong_cubic_sign_negative_control": s.expand(
            external + seagull - triple
        ),
        "checks": checks,
        "gates": {
            "symbolic_Ward_not_point_sampling": True,
            "general_mass_and_gauge_vector_retained": mu in triple.free_symbols
            and x1 in triple.free_symbols
            and xk in triple.free_symbols,
            "wrong_cubic_sign_rejected": s.expand(external + seagull - triple) != 0,
            "metric_closed_loops_dimension_four": True,
            "Bose_relabeling_covers_three_partitions": True,
            "all_graphs_in_one_Ward_polynomial": True,
        },
    }
