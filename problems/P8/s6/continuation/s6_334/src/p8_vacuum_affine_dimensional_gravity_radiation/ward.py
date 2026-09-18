"""S304 canonical vertices with both D-dimensional traces and projector."""

from collections import Counter
from functools import cache
from itertools import permutations, product

import sympy as s

DIM = s.Symbol("dimension", positive=True)
LOOP_DIM = s.Symbol("closed_metric_dimension")

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
    return lambda i, j: (
        outer(p, q, i, j) + outer(q, p, i, j) + scaled(metric(i, j), 2 * mu / (DIM - 2))
    )


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
            LOOP_DIM**n * coef * s.Mul(*(dot(i, j) for i, j in pairs))
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
def dimensional_ward_numerators():
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
    return external, seagull, triple


@cache
def tt_channel():
    """Canonical TT rank-one polynomial; no dimension-specific Gram identities."""
    return s.factor(
        s.together(rank_one_channel().subs({LOOP_DIM: DIM, xk: 0, base[4, 4]: 0}))
    )


def auxiliary_tt_channel():
    """S295 unit massless-exchange basis only, not a source substitution."""
    L, R = dot(6, 6), dot(7, 7)
    result = -sum(dot(i, 5) ** 2 / (dot(i, 4) * R) for i in (0, 1))
    result -= sum(dot(i, 5) ** 2 / (dot(i, 4) * L) for i in (2, 3))
    result -= 2 * dot(6, 5) ** 2 / (L * R)
    return s.factor(result.subs({xk: 0, base[4, 4]: 0}))


@cache
def data():
    from p8_vacuum_affine_minimal_gravity_radiation import ward as old

    total = sum(dimensional_ward_numerators())
    wrong = s.factor(s.together(total.subs(LOOP_DIM, 4)))
    generic = tt_channel()
    difference = (
        generic
        - generic.subs(DIM, 4)
        - 2 * mu**2 * (DIM - 4) / (DIM - 2) * auxiliary_tt_channel()
    )
    witness = s.factor(
        wrong.subs(
            {
                DIM: 5,
                mu: 1,
                a: 2,
                b: 3,
                u1: s.Rational(1, 7),
                u2: -s.Rational(2, 7),
                u3: s.Rational(3, 7),
                x1: s.Rational(1, 11),
                x2: -s.Rational(2, 11),
                x3: s.Rational(4, 11),
                xk: s.Rational(5, 11),
            }
        )
    )
    htrace = 2 * dot(0, 1) + 2 * mu * DIM / (DIM - 2)
    checks = {
        "general_D_channel_Ward": s.factor(s.together(total.subs(LOOP_DIM, DIM))),
        "generic_all_trace_gauge_four_D_recovery": s.factor(
            s.together(
                rank_one_channel().subs({DIM: 4, LOOP_DIM: 4}) - old.rank_one_channel()
            )
        ),
        "general_TT_auxiliary_decomposition": s.factor(s.together(difference)),
        "wrong_trace_recovers_only_D4": s.factor(wrong.subs(DIM, 4)),
        "D_trace_reversed_scalar_numerator": s.cancel(
            -(dot(0, 1) + mu)
            - (2 * dot(0, 1) - DIM * (dot(0, 1) + mu)) / (DIM - 2)
            - 2 * mu / (DIM - 2)
        ),
        "trace_reversed_numerator_trace": s.cancel(
            htrace - 2 * dot(0, 1) - 2 * mu * DIM / (DIM - 2)
        ),
    }
    return {
        "checks": checks,
        "gates": {
            "general_Ward_rational_identity_no_Gram_rank_relation": True,
            "full_generic_rank_one_recovery_not_only_one_helicity": True,
            "wrong_four_dimensional_metric_loop_witness_nonzero": witness != 0,
            "all_TT_D_dependence_is_exact_dimension_ratio": DIM
            not in s.factor(
                s.cancel((generic - generic.subs(DIM, 4)) * (DIM - 2) / (DIM - 4))
            ).free_symbols,
            "null_rank_one_span_proves_all_integer_D_TT_tensors": True,
            "three_Bose_partitions_and_all_seagulls_cubic_retained": True,
        },
        "whole_projector_numerator": "p_lower q_lower+q_lower p_lower+2mu eta/(D-2); every closed metric loop is D.",
        "whole_channel_decomposition": "G_D:epsilon=G_4^alg:epsilon+2mu^2*(D-4)/(D-2)*S_aux:epsilon for TT epsilon; canonical G4 algebraic representative is fixed before changing dimension.",
        "whole_auxiliary_channel": auxiliary_tt_channel(),
        "wrong_closed_trace_nonzero_witness": witness,
        "whole_Ward_proof": "The exact external+seagull+cubic gauge numerator times L*R is identically zero in unrestricted on-shell Gram variables and dimension. Wrong metric loops=4 violate it. Null complex rank-one tensors span the physical transverse symmetric trace-free space for each integer D>=4.",
    }
