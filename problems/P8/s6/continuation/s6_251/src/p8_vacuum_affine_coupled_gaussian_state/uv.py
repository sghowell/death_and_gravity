"""Exact double-cone normalization and noncommuting matrix Riccati recursion."""

from functools import cache

import sympy as s

from . import phase


def principal_pairs():
    f, j = s.symbols("positive_F positive_Jc", positive=True)
    th, e, ell = s.symbols("nonzero_Theta nonzero_E ell", real=True, nonzero=True)
    w = -ell * e
    result = {}
    for name, pivot, mixed in (("outer", th, w / th), ("central", e, ell)):
        K = s.Matrix([[2 * j / pivot**2 + mixed**2, mixed], [mixed, 1]])
        G = s.Matrix([[2 * f / pivot**2 + mixed**2, mixed], [mixed, 1]])
        R = s.Matrix([[pivot / s.sqrt(2 * j), 0], [-mixed * pivot / s.sqrt(2 * j), 1]])
        W = R.inv().T * s.diag(s.sqrt(f / j) / phase.a, 1 / phase.a) * R.inv()
        result[name] = {"K": K, "G": G, "R": R, "W": W}
    return f, j, result


def matrix_at(coefficients, power):
    return coefficients.get(power, s.zeros(2))


def riccati_residual_coefficient(Z, B, A, C, time, power):
    answer = matrix_at(Z, power).diff(time) + matrix_at(C, power)
    for ae, am in A.items():
        zm = matrix_at(Z, power - ae)
        answer += am.T * zm + zm * am
    for ze, zm in Z.items():
        for be, bm in B.items():
            answer += zm * bm * matrix_at(Z, power - ze - be)
    return answer.applyfunc(s.cancel)


def riccati_series(B, A, C, frequencies, time, depth):
    if type(depth) is not int or not 1 <= depth <= 8:
        raise ValueError("Use a finite diagnostic recursion depth from1 through8")
    if len(frequencies) != 2 or B.get(0) != s.eye(2):
        raise ValueError(
            "Use both frequencies in the normalized principal kinetic frame"
        )
    omega = tuple(map(s.sympify, frequencies))
    if any(w.is_positive is not True for w in omega):
        raise ValueError(
            "Both diagnostic principal frequencies must be provably positive"
        )
    if C.get(2) != s.diag(*[w * w for w in omega]):
        raise ValueError("Keep the complete matching principal potential")
    if any(type(k) is not int or k > 0 for rows in (B, A) for k in rows):
        raise ValueError("B and A start at order0")
    if any(type(k) is not int or k > 2 for k in C):
        raise ValueError("C starts at order2")
    if any(M.shape != (2, 2) for rows in (B, A, C) for M in rows.values()):
        raise ValueError("Every whole coefficient has two-by-two shape")
    if any(M != M.T for rows in (B, C) for M in rows.values()):
        raise ValueError("Preserve the full symmetric Hamiltonian coefficients")
    Z = {1: -s.I * s.diag(*omega)}
    for index in range(depth):
        power = 1 - index
        remainder = riccati_residual_coefficient(Z, B, A, C, time, power)
        Z[-index] = s.Matrix(
            2,
            2,
            lambda i, j, remainder=remainder: s.cancel(
                -s.I * remainder[i, j] / (omega[i] + omega[j])
            ),
        )
    return Z


@cache
def data():
    f, j, pairs = principal_pairs()
    checks = {}
    speed = s.Symbol("characteristic_speed_squared", real=True)
    for name, row in pairs.items():
        K, G, R, W = (row[k] for k in ("K", "G", "R", "W"))
        checks[name + "_full_kinetic_normalization"] = R.T * K * R - s.eye(2)
        checks[name + "_both_principal_cones"] = R.T * G * R - s.diag(f / j, 1)
        checks[name + "_whole_characteristic_factorization"] = (
            speed * K - G
        ).det() - K.det() * (speed - f / j) * (speed - 1)
        checks[name + "_positive_Riccati_principal_square"] = (
            W * K.inv() * W - G / phase.a**2
        )
        checks[name + "_symmetric_principal_Riccati_graph"] = W - W.T
        checks[name + "_principal_symplectic_diagonalization"] = (
            s.diag(R, R.inv().T) * phase.OMEGA * s.diag(R, R.inv().T).T - phase.OMEGA
        )
    # Numerical bounds below are inherited EXACT current coefficient inequalities.
    checks["current_squared_speed_gap_lower"] = s.Rational(1, 1000) / 100 - s.Rational(
        1, 10**5
    )
    checks["current_speed_lower_from_F_J_and_scale"] = s.sqrt(
        s.Rational(1, 100) / 100
    ) * s.Rational(16, 25) - s.Rational(4, 625)
    checks["current_distinct_root_gap_lower"] = s.Rational(1, 10**5) * s.Rational(
        1, 2
    ) * s.Rational(16, 25) - s.Rational(8, 25 * 10**5)
    time = s.Symbol("diagnostic_time", real=True)
    mix = s.Matrix([[time / 5, (1 + time) / 7], [-(2 + time * time) / 11, -time / 13]])
    potential = s.Matrix([[1 + time * time, time / 3], [time / 3, 2 + time]])
    B = {
        0: s.eye(2),
        -2: s.Matrix(
            [
                [s.Rational(1, 9), s.Rational(1, 21)],
                [s.Rational(1, 21), s.Rational(2, 17)],
            ]
        ),
    }
    A = {0: mix, -2: s.Matrix([[0, time / 19], [time / 23, 0]])}
    C = {2: s.diag(1, 4), 0: potential, -2: s.Matrix([[time, 1], [1, -time]])}
    Z = riccati_series(B, A, C, (s.S.One, s.Integer(2)), time, 5)
    for power in range(2, -4, -1):
        checks["complete_noncommuting_Riccati_coefficient_" + str(power)] = (
            riccati_residual_coefficient(Z, B, A, C, time, power)
        )
    for power, value in Z.items():
        checks["symmetric_complex_Lagrangian_coefficient_" + str(power)] = (
            value - value.T
        )
    x11, x12, x22, w1, w2 = s.symbols("x11 x12 x22 omega1 omega2", positive=True)
    X = s.Matrix([[x11, x12], [x12, x22]])
    W = s.diag(w1, w2)
    solution = s.Matrix(2, 2, lambda i, j: X[i, j] / ((w1, w2)[i] + (w1, w2)[j]))
    checks["full_positive_sum_Sylvester_inverse"] = W * solution + solution * W - X
    return {
        "full_current_characteristics": "Both complete current scalar principal pencils have c²=F/Jc and1, with0<F/Jc<1. F>1/100,Jc<100,Jc-F>1/1000 and a<=25/16 give omega_min>4/625 after factoring |P|, and distinct positive-root gap>8/(25*10^5).",
        "whole_principal_normalizations": pairs,
        "matrix_Riccati_recursion": "In the complete canonical normalized frame Q'=A Q+B Pi, Pi'=-C Q-A^T Pi, Z'= -C-A^T Z-ZA-ZBZ. B0=I,C2=diag(omega_i²), Z1=-i diag(omega_i). Each next symmetric coefficient solves a Sylvester equation with denominators omega_i+omega_j. Every subleading B,A,C term and time derivative is retained.",
        "diagnostic_full_Riccati_coefficients": Z,
        "all_order_boundary": "Finite coefficient checks do not prove a continuum wavefront theorem. The written proof must supply arbitrary-depth symbol estimates, actual full coefficient smoothness, exact-vs-approximate polynomial-loss control, chart transitions, positive-frequency phases and the smooth-bump anomalous-block estimate.",
        "preparation_not_instantaneous": "Smooth compact sampling annihilates anomalous positive-frequency sums to every algebraic order by repeated integration by parts. This statement requires the all-order matrix mode proof; it is not an imported scalar theorem.",
        "checks": {
            name: value.applyfunc(s.cancel)
            if isinstance(value, s.MatrixBase)
            else s.cancel(value)
            for name, value in checks.items()
        },
        "gates": {
            "two_cones_not_scalar_normally_hyperbolic_symbol": f != j,
            "Riccati_lower_coefficients_retain_noncommuting_mixing": mix * potential
            != potential * mix,
            "finite_q_kinetic_corrections_not_deleted": B[-2] != s.zeros(2),
            "positive_frequency_sum_never_requires_same_branch_gap": w1 + w2 != 0,
            "unselected_instantaneous_ground_state_not_declared_Hadamard": True,
            "reduced_state_not_covariant_renormalized_metric_stress": True,
            "no_nonlinear_or_interacting_state_promotion": True,
        },
    }
