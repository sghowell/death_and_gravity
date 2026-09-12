from functools import cache

import sympy as s

A = s.symbols("a0:5", positive=True)
Dtime = s.symbols("D0:5", real=True)
Gtime = s.symbols("G0:5", real=True)
p = s.Symbol("p", real=True)
ZERO = (s.Integer(0),) * 4


def add(*rows):
    return tuple(s.expand(sum(row[j] for row in rows)) for j in range(4))


def scale(row, value):
    return tuple(s.expand(v * value) for v in row)


def mul(a, b):
    if a == ZERO or b == ZERO:
        return ZERO
    return (
        s.expand(a[0] * b[0]),
        s.expand(a[0] * b[1] + a[1] * b[0]),
        s.expand(a[0] * b[2] + a[2] * b[0]),
        s.expand(a[0] * b[3] + a[1] * b[2] + a[2] * b[1] + a[3] * b[0]),
    )


def scalar(value):
    return (s.sympify(value), 0, 0, 0)


def time_derivative(value):
    return s.expand(
        sum(
            s.diff(value, jets[j]) * jets[j + 1]
            for jets in (A, Dtime, Gtime)
            for j in range(4)
        )
    )


def matrix_product(a, b):
    return [
        [add(*(mul(a[i][k], b[k][j]) for k in range(len(b)))) for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def trace_square(matrix):
    return add(
        *(
            mul(matrix[i][j], matrix[j][i])
            for i in range(len(matrix))
            for j in range(len(matrix))
        )
    )


def literal_curvatures(D, G):
    dim = D.rows
    N = dim + 1
    a, a1, a2 = A[:3]
    S = (D * G + G * D) / 2

    def e(i, j, r=0, z=0):
        if i == 0 or j == 0:
            return scalar(-1) if i == j == 0 and r == z == 0 else ZERO
        i, j = i - 1, j - 1
        if z == 0:
            mixed = sum(
                s.binomial(r, k) * Dtime[k] * Gtime[r - k] for k in range(r + 1)
            )
            return (
                s.Integer(i == j) if r == 0 else 0,
                Dtime[r] * D[i, j],
                Gtime[r] * G[i, j],
                mixed * S[i, j],
            )
        return (
            0,
            (-s.I * p) ** z * Dtime[r] * D[i, j],
            (s.I * p) ** z * Gtime[r] * G[i, j],
            0,
        )

    def metric(i, j, r=0, z=0):
        if r == 0:
            return scale(e(i, j, 0, z), a * a)
        if r == 1:
            return add(scale(e(i, j, 0, z), 2 * a * a1), scale(e(i, j, 1, z), a * a))
        if r == 2:
            return add(
                scale(e(i, j, 0, z), 2 * (a1 * a1 + a * a2)),
                scale(e(i, j, 1, z), 4 * a * a1),
                scale(e(i, j, 2, z), a * a),
            )
        raise ValueError("Only metric jets through two")

    @cache
    def derivative(i, j, indices):
        if any(k not in (0, dim) for k in indices):
            return ZERO
        return metric(i, j, indices.count(0), indices.count(dim))

    inverse = []
    for i in range(N):
        row = []
        for j in range(N):
            if i == 0 or j == 0:
                value = scalar(-1 / a**2) if i == j == 0 else ZERO
            else:
                value = scale(
                    (
                        s.Integer(i == j),
                        -Dtime[0] * D[i - 1, j - 1],
                        -Gtime[0] * G[i - 1, j - 1],
                        Dtime[0] * Gtime[0] * S[i - 1, j - 1],
                    ),
                    a**-2,
                )
            row.append(value)
        inverse.append(row)

    @cache
    def christoffel(i, j, k):
        return scale(
            add(
                derivative(i, k, (j,)),
                derivative(i, j, (k,)),
                scale(derivative(j, k, (i,)), -1),
            ),
            s.Rational(1, 2),
        )

    @cache
    def riemann(i, j, k, l):
        first = scale(
            add(
                derivative(i, l, (j, k)),
                derivative(j, k, (i, l)),
                scale(derivative(i, k, (j, l)), -1),
                scale(derivative(j, l, (i, k)), -1),
            ),
            s.Rational(1, 2),
        )
        terms = []
        for r in range(N):
            for v in range(N):
                if inverse[r][v] != ZERO:
                    terms.append(
                        mul(
                            inverse[r][v],
                            add(
                                mul(christoffel(r, j, k), christoffel(v, i, l)),
                                scale(
                                    mul(christoffel(r, j, l), christoffel(v, i, k)), -1
                                ),
                            ),
                        )
                    )
        return add(first, *terms)

    Ric = [
        [
            add(
                *(
                    mul(inverse[i][k], riemann(i, j, k, l))
                    for i in range(N)
                    for k in range(N)
                )
            )
            for l in range(N)
        ]
        for j in range(N)
    ]
    R = add(*(mul(inverse[i][j], Ric[i][j]) for i in range(N) for j in range(N)))
    Ric2 = trace_square(matrix_product(inverse, Ric))
    pairs = [(i, j) for i in range(N) for j in range(i + 1, N)]
    RB = [[riemann(i, j, k, l) for k, l in pairs] for i, j in pairs]
    GB = [
        [
            add(
                mul(inverse[i][k], inverse[j][l]),
                scale(mul(inverse[i][l], inverse[j][k]), -1),
            )
            for k, l in pairs
        ]
        for i, j in pairs
    ]
    Riem2 = scale(trace_square(matrix_product(GB, RB)), 4)
    volume = (
        s.Integer(1),
        s.trace(D) * Dtime[0] / 2,
        s.trace(G) * Gtime[0] / 2,
        s.trace(D) * s.trace(G) * Dtime[0] * Gtime[0] / 4,
    )
    return {
        name: s.factor(a ** (dim + 1) * mul(volume, row)[3])
        for name, row in (
            ("R_old", R),
            ("R_squared", mul(R, R)),
            ("Ricci_squared", Ric2),
            ("Riemann_squared", Riem2),
        )
    }


def spatial_operator(expression):
    expression = s.expand(expression - expression.subs(p, 0))
    out = 0
    for r in range(3):
        term = s.diff(expression, Dtime[r])
        for _ in range(r):
            term = -time_derivative(term)
        out += term
    return s.factor(out)
