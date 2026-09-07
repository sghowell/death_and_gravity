"""Generic exact rational Cauchy-symbol construction.

All polynomial coefficients retain u and c as indeterminates. No rational
time fixture, interpolation or numerical rank test is used in this route.
The small native rational engine changes arithmetic, not the pinned scalar
Hamiltonian, canonical boundary or physical observable definitions.
"""

from functools import cache

import sympy as sp
from p8_variable_constraints import action, reduction

from . import leading, physical
from . import rational as ra


@cache
def phase():
    d, r = action.derive(), reduction.derive()
    variables = (d["u"], d["c"], action.K)
    def convert(matrix):
        return [[ra.from_sympy(x, variables) for x in row] for row in matrix.tolist()]
    original = convert(sp.hessian(d["H"], d["state"]))
    transformation = convert(r["old_state"].subs(reduction.Q, 0).jacobian(reduction.EXTENDED))
    boundary = convert(sp.hessian(r["H_time_correction"], reduction.EXTENDED))
    extended = ra.add(ra.multiply(ra.multiply(ra.transpose(transformation), original), transformation), boundary)
    D = extended[6][6]
    expected_D = ra.from_sympy(r["expected_secondary_D"], variables)
    if D != expected_D or D.degree() != 0:
        raise ValueError("Native generic secondary-constraint identity failed")
    H = [[extended[i][j]-extended[i][6]*extended[6][j]/D for j in range(6)] for i in range(6)]
    J = [[ra.Rational(reduction.J[i, j]) for j in range(6)] for i in range(6)]
    O = convert(physical.rows())
    A = ra.multiply(J, H)
    C = O+ra.add(ra.differentiate(O), ra.multiply(O, A))
    return {"variables": variables, "D": D, "H": H, "A": A, "O": O, "C": C}


@cache
def derive():
    d = phase()
    A, O, C = (d[key] for key in ("A", "O", "C"))
    inverse = ra.inverse(C)
    acceleration = ra.add(ra.add(ra.differentiate(ra.differentiate(O)),
        [[2*x for x in row] for row in ra.multiply(ra.differentiate(O), A)]),
        ra.multiply(O, ra.add(ra.differentiate(A), ra.multiply(A, A))))
    equation = ra.multiply(acceleration, inverse)
    J = [[ra.Rational(reduction.J[i, j]) for j in range(6)] for i in range(6)]
    symplectic = ra.multiply(ra.multiply(ra.transpose(inverse), J), inverse)
    return {**d, "inverse_C": inverse, "det_C": ra.determinant(C), "Mq": [row[:3] for row in equation],
            "Mv": [row[3:] for row in equation], "symplectic": symplectic}


def degree_table(matrix):
    return [[x.degree() for x in row] for row in matrix]


@cache
def checks():
    d, l = derive(), leading.derive()
    variables = d["variables"]
    kinetic = [[ra.from_sympy(x, variables) for x in row] for row in l["kinetic"].tolist()]
    gradient = [[ra.from_sympy(x, variables) for x in row] for row in l["gradient"].tolist()]
    speed = ra.multiply(ra.inverse(kinetic), gradient)
    momentum = ra.Rational(ra.CONTEXT.gens()[2])
    remainders = [[d["Mq"][i][j]+momentum*speed[i][j] for j in range(3)] for i in range(3)]
    form_difference = [[d["symplectic"][i][j+3]-kinetic[i][j] for j in range(3)] for i in range(3)]
    determinant_identity = d["det_C"].leading()*ra.determinant(kinetic)-1
    if determinant_identity:
        raise ValueError("The generic physical Cauchy determinant/kinetic identity failed")
    groups = {
        "Mq_plus_K_kinetic_inverse_gradient": (remainders, 0),
        "Mv": (d["Mv"], 0),
        "Sigma_qv_minus_kinetic": (form_difference, -1),
        "Sigma_vv": ([row[3:] for row in d["symplectic"][3:]], -2),
        "Sigma_qq": ([row[:3] for row in d["symplectic"][:3]], 0),
    }
    for name, (matrix, degree) in groups.items():
        if any(x.degree() > degree for row in matrix for x in row):
            raise ValueError(f"The generic uniform symbol degree bound failed: {name}")
    if d["det_C"].degree() != 0:
        raise ValueError("The physical Cauchy determinant lost its generic leading degree")
    return {"groups": {name: {"degrees": degree_table(matrix), "upper_degree": bound}
                       for name, (matrix, bound) in groups.items()},
            "determinant_leading_kinetic_identity": determinant_identity,
            "det_C_degree": d["det_C"].degree()}
