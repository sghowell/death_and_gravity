"""Fraction/Taylor physical Cauchy fixtures and Bernstein conversion.

The frozen independent scalar engine derives its full phase jets by
literal quadratic polarization. This module does not import the new
SymPy/FLINT symbol, principal matrices or degree verdicts. The Bernstein
routine independently converts supplied coefficient records; it does not
claim to independently derive the supplied margin polynomial.
"""

from fractions import Fraction as Q
from itertools import product
from math import comb

from p8_variable_beta.independent import Jet, exact, family
from p8_variable_constraints.independent import (
    inverse,
    matmul,
    scalar_jet_fixture,
    transpose,
)


def add(left, right):
    return [[x+y for x, y in zip(row, other, strict=True)] for row, other in zip(left, right, strict=True)]


def determinant(matrix):
    size = len(matrix)
    rows = [list(row) for row in matrix]
    result = Q(1)
    for j in range(size):
        pivot = next((i for i in range(j, size) if rows[i][j]), None)
        if pivot is None:
            return Q()
        if pivot != j:
            rows[j], rows[pivot] = rows[pivot], rows[j]
            result = -result
        coefficient = rows[j][j]
        result *= coefficient
        for i in range(j+1, size):
            ratio = rows[i][j]/coefficient
            for k in range(j+1, size):
                rows[i][k] -= ratio*rows[j][k]
    return result


def fixture(time, lapse, momentum_squared):
    u, c, K = map(exact, (time, lapse, momentum_squared))
    if not (abs(u) <= Q(1, 10) and 2 < c <= 4 and K > 0):
        raise ValueError("Outside the actual positive-link local microlocal family")
    phase = scalar_jet_fixture(u, c, K)
    f = family(Jet(u, 1), c)
    a, w = f["a"], f["chi_speed"]
    zero = Jet(0)
    pi = [zero, Jet(-1/K), zero, zero, zero, zero]
    xi = [zero, zero, 3*a*a*w/(2*K), zero, -3/(2*a*K), zero]
    chi = [w*x for x in xi]
    chi[2] += 1
    rows = [pi, xi, chi]
    O0 = [[x.value for x in row] for row in rows]
    O1 = [[x.first for x in row] for row in rows]
    O2 = [[2*x.second_coefficient for x in row] for row in rows]
    J = [[Q(int(i < 3 and j == i+3)-int(i >= 3 and j == i-3)) for j in range(6)] for i in range(6)]
    A, Ap = matmul(J, phase["H0"]), matmul(J, phase["H1"])
    C = O0+add(O1, matmul(O0, A))
    inv = inverse(C)
    acceleration = add(add(O2, [[2*x for x in row] for row in matmul(O1, A)]),
                       matmul(O0, add(Ap, matmul(A, A))))
    equation = matmul(acceleration, inv)
    form = matmul(matmul(transpose(inv), J), inv)
    return {"u": u, "c": c, "K": K, "C": C, "inverse_C": inv,
            "Mq": [row[:3] for row in equation], "Mv": [row[3:] for row in equation],
            "symplectic": form, "det_C": determinant(C)}


def fixtures():
    return [fixture(*case) for case in (
        (0, 4, 1), (0, 3, Q(7, 2)), (0, Q(201, 100), 100),
        (Q(1, 100), 4, 7), (Q(-1, 100), 4, 7),
        (Q(1, 10), 4, 500_000),
    )]


def bernstein_from_time_lapse_terms(terms):
    """Independent exact u²=x/100,c=2+2z conversion of supplied terms."""
    converted = {}
    for power, coefficient in terms.items():
        du, dc = map(int, power)
        if du % 2:
            raise ValueError("A margin polynomial must be even in u")
        coefficient = exact(coefficient)*Q(1, 100)**(du//2)
        for j in range(dc+1):
            key = (du//2, j)
            converted[key] = converted.get(key, Q())+coefficient*comb(dc, j)*2**dc
    converted = {key: value for key, value in converted.items() if value}
    nu, nc = max(key[0] for key in converted), max(key[1] for key in converted)
    result = []
    for i, j in product(range(nu+1), range(nc+1)):
        result.append(sum((value*Q(comb(i, p), comb(nu, p))*Q(comb(j, q), comb(nc, q))
                           for (p, q), value in converted.items() if p <= i and q <= j), Q()))
    return (nu, nc), tuple(result)
