"""Invertible two-trace curl matrices: full longitudinal and sign test.

The fixed rolling Schur matrix is independently interfaced to the literal
connection solve. No open-tube nonlinear 52-component inverse is assumed.
"""
from functools import cache

import sympy as sp

R, S, T = sp.symbols("R S T", real=True)
Q = sp.Symbol("q", positive=True)
D = sp.ImmutableMatrix([[sp.Rational(3, 8), -sp.Rational(21, 8)],
                        [-sp.Rational(21, 8), sp.Rational(3, 8)]])
Z = sp.ImmutableMatrix([[R, T], [T, S]])


def _matrix(value):
    return sp.ImmutableMatrix(value.applyfunc(sp.factor))


@cache
def longitudinal():
    """Eliminate both temporal components before reading the kinetic matrix."""
    velocity = sp.Matrix(sp.symbols("sigma1_dot sigma2_dot", real=True))
    temporal = sp.Matrix(sp.symbols("W01 W02", real=True))
    lagrangian = ((temporal.T*D.inv()*temporal)[0]
                  +Q*((velocity-temporal).T*Z*(velocity-temporal))[0])/2
    hessian = D.inv()+Q*Z
    solution = _matrix(hessian.inv()*Q*Z*velocity)
    reduced = sp.factor(lagrangian.subs(dict(zip(temporal, solution, strict=True)), simultaneous=True))
    kinetic = _matrix(sp.hessian(reduced, tuple(velocity))/2)
    parallel = _matrix((D+Z.inv()/Q).inv()/2)
    minus, plus = sp.Matrix([1, 1]), sp.Matrix([1, -1])
    threshold = sp.factor(2*(minus.T*Z.inv()*minus)[0]/9)
    return {
        "Z": Z, "D": D, "temporal_hessian": hessian,
        "temporal_solution": solution, "kinetic": kinetic,
        "parallel_sum": parallel,
        "threshold": threshold, "negative_test_vector": minus,
        "positive_test_vector": plus,
        "euler_residual": _matrix(hessian*solution-Q*Z*velocity),
        "parallel_sum_residual": _matrix(kinetic-parallel),
        "negative_quadratic": sp.factor((minus.T*(D+Z.inv()/Q)*minus)[0]),
        "positive_quadratic": sp.factor((plus.T*(D+Z.inv()/Q)*plus)[0]),
    }


def exact_matrix(r, s, t):
    values = []
    for value in (r, s, t):
        if isinstance(value, (bool, float, str)):
            raise TypeError("Kinetic matrix entries must be exact finite reals")
        value = sp.sympify(value)
        if (not isinstance(value, sp.Expr) or value.has(sp.Float) or value.free_symbols
                or value.is_real is not True or value.is_finite is not True):
            raise ValueError("Kinetic matrix entries must be exact finite real constants")
        values.append(value)
    return sp.ImmutableMatrix([[values[0], values[2]], [values[2], values[1]]])


def classify(r, s, t):
    """Exact rank/sign chart, not a floating eigenvalue tolerance."""
    matrix = exact_matrix(r, s, t)
    determinant, trace = sp.factor(matrix.det()), sp.trace(matrix)
    if matrix == sp.zeros(2):
        return {"rank": 0, "chart": "unchanged_auxiliary"}
    if determinant.is_zero is True:
        if matrix[0, 0].is_zero is False:
            coupling = matrix[0, 0]
            coefficient = sp.Matrix([1, matrix[0, 1]/coupling])
        elif matrix[0, 0].is_zero is True:
            coupling, coefficient = matrix[1, 1], sp.Matrix([0, 1])
        else:
            raise ValueError("Exact rank-one pivot is unproved")
        if matrix != coupling*coefficient*coefficient.T:
            raise ArithmeticError("The rank-one reconstruction failed")
        return {"rank": 1, "chart": "rank_one_requires_own_Schur_chart",
                "coupling": coupling, "coefficient": coefficient}
    if determinant.is_zero is not False:
        raise ValueError("Exact determinant rank is unproved")
    positive = matrix[0, 0].is_positive is True and determinant.is_positive is True
    if not positive and not (determinant.is_negative is True or trace.is_negative is True):
        raise ValueError("Exact kinetic inertia is unproved")
    return {"rank": 2, "chart": "positive_transverse_test_passed" if positive else "negative_transverse_direction",
            "matrix": matrix, "determinant": determinant}


@cache
def checks():
    data = longitudinal()
    negative, positive = data["negative_test_vector"], data["positive_test_vector"]
    ann, mix1, mix2, theta = sp.symbols("ann mix1 mix2 theta", real=True)
    temporal = data["temporal_hessian"]
    auxiliary = sp.Matrix([[ann, 2*Q*theta, mix1, mix2],
                           [2*Q*theta, 0, 0, 0],
                           [mix1, 0, temporal[0, 0], temporal[0, 1]],
                           [mix2, 0, temporal[1, 0], temporal[1, 1]]])
    return {
        "D_positive_eigenvector": D*positive-3*positive,
        "D_negative_eigenvector": D*negative+sp.Rational(9, 4)*negative,
        "D_nonzero_determinant": sp.factor(D.det()+sp.Rational(27, 4)),
        "temporal_euler": data["euler_residual"],
        "longitudinal_parallel_sum": data["parallel_sum_residual"],
        "negative_ray_threshold": sp.factor(data["negative_quadratic"]
                                              -sp.Rational(9, 2)*(data["threshold"]/Q-1)),
        "positive_ray": sp.factor(data["positive_quadratic"]-6
                                    -(positive.T*Z.inv()*positive)[0]/Q),
        "temporal_chart_equivalence": sp.factor((D.inv()+Q*Z).det()
                                                 -Q**2*Z.det()/D.det()*(D+Z.inv()/Q).det()),
        "full_lapse_shift_temporal_rank": sp.factor(auxiliary.det()+4*Q**2*theta**2*temporal.det()),
    }


@cache
def calibration():
    z = sp.Matrix([[2, 1], [1, 3]])
    point = {R: 2, S: 3, T: 1, Q: 2}
    data = longitudinal()
    return {"Z": z, "threshold": data["threshold"].subs(point),
            "q": 2, "longitudinal_kinetic": _matrix(data["kinetic"].subs(point)),
            "longitudinal_determinant": sp.factor(data["kinetic"].subs(point).det()),
            "rank_two_positive_Z_not_sufficient": True,
            "nonlinear_open_tube_inverse_claim": False}
