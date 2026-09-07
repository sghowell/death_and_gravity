"""Exact time-normalization and principal symmetrizer identities.

The proof supplies finite uniform remainder norms at fixed delta_min; no
numerical value of their compact supremum is fabricated here.
"""

from functools import cache

from . import leading, symbol
from . import rational as ra


def scale(matrix, coefficient):
    return [[coefficient*x for x in row] for row in matrix]


@cache
def derive():
    s, l = symbol.derive(), leading.derive()
    variables = s["variables"]
    def convert(matrix):
        return [[ra.from_sympy(x, variables) for x in row] for row in matrix.tolist()]
    U = ra.from_sympy(l["U"], variables)
    T = [[ra.Rational(int(i == j)) for j in range(3)] for i in range(3)]
    T[0][0] = U
    inverse = ra.inverse(T)
    first = ra.multiply(ra.differentiate(T), inverse)
    second = ra.multiply(ra.differentiate(ra.differentiate(T)), inverse)
    velocity_conjugate = ra.multiply(ra.multiply(T, s["Mv"]), inverse)
    Mv = ra.add(velocity_conjugate, scale(first, 2))
    Mq = ra.add(ra.add(ra.multiply(ra.multiply(T, s["Mq"]), inverse),
                      scale(ra.multiply(velocity_conjugate, first), -1)),
                ra.add(second, scale(ra.multiply(first, first), -2)))
    kinetic, gradient = convert(l["normalized_kinetic"]), convert(l["normalized_gradient"])
    speed = ra.multiply(ra.inverse(kinetic), gradient)
    K = ra.Rational(ra.CONTEXT.gens()[2])
    remainder = ra.add(Mq, scale(speed, K))
    S, principal = ra.zeros(6, 6), ra.zeros(6, 6)
    for i in range(3):
        principal[i][i+3] = ra.Rational(1)
        for j in range(3):
            S[i][j], S[i+3][j+3] = gradient[i][j], kinetic[i][j]
            principal[i+3][j] = -speed[i][j]
    symmetric_residual = ra.add(ra.multiply(S, principal), ra.multiply(ra.transpose(principal), S))
    return {"T": T, "first_connection": first, "second_connection": second,
            "normalized_Mq": Mq, "normalized_Mv": Mv,
            "normalized_Rq": remainder, "kinetic": kinetic, "gradient": gradient,
            "energy_weight": S, "principal_generator": principal,
            "principal_energy_skew_residual": symmetric_residual,
            "energy_rate_formula": "C_delta=4 sup||S'||+512 sup(||R||+||V||), finite at each fixed positive delta_min"}


def checks():
    d = derive()
    if any(x for row in d["principal_energy_skew_residual"] for x in row):
        raise ValueError("The exact order-k energy skew identity failed")
    if any(x.degree() > 0 for matrix in (d["normalized_Rq"], d["normalized_Mv"]) for row in matrix for x in row):
        raise ValueError("A time normalization connection changed the uniform order-one remainder bound")
    return {"energy_skew_identity_count": 36,
            "normalized_Rq_degrees": symbol.degree_table(d["normalized_Rq"]),
            "normalized_Mv_degrees": symbol.degree_table(d["normalized_Mv"]),
            "energy_rate_formula": d["energy_rate_formula"]}
