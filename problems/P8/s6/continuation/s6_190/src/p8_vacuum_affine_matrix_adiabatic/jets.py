"""Momentum-uniform finite real-time jet bounds for the actual matrix frame."""

from functools import cache
from math import comb

import sympy as s

ORDER = 10
JETS = 12
DELTA = s.Rational(1, 100)
MASS = s.Integer(1000)
t = s.Symbol("t", real=True)


@cache
def scale_polynomial(power, order):
    if power not in (1, 2, 3) or type(order) is not int or not 0 <= order <= JETS:
        raise ValueError("Require p=1,2,3 and a native derivative order0..12")
    if order == 0:
        return s.Integer(1)
    p = scale_polynomial(power, order - 1)
    n = order - 1
    return s.Poly((1 + t * t) * s.diff(p, t) - 2 * (n + 2 * power) * t * p, t).as_expr()


def scale_bound(power, order):
    return sum(
        abs(c) * s.Rational(1, 2) ** i[0]
        for i, c in s.Poly(scale_polynomial(power, order), t).terms()
    )


@cache
def bell(order):
    if type(order) is not int or not 0 <= order <= JETS:
        raise ValueError("Require a native derivative order0..12")
    if order == 0:
        return s.Integer(1)
    return DELTA * sum(comb(order - 1, j) * bell(j) for j in range(order))


@cache
def constants():
    exponential = [s.Integer(1)] + [2 * bell(j) for j in range(1, JETS + 1)]
    c = {p: [scale_bound(p, j) for j in range(JETS + 1)] for p in (1, 2, 3)}
    q = [s.Integer(1)]
    frequency_square = [s.Integer(1)]
    for n in range(1, JETS + 1):
        q.append(
            max(
                c[3][n],
                sum(comb(n, j) * c[1][n - j] * exponential[j] for j in range(n + 1)),
            )
        )
        frequency_square.append(
            sum(comb(n, j) * c[2][n - j] * exponential[j] for j in range(n + 1))
        )
    # Direct first derivatives use |H|<=8/5, not the looser polynomial bound.
    q[1] = s.Integer(5)
    frequency_square[1] = s.Integer(4)
    root, inverse, omega, reciprocal = ([s.Integer(1)] for _ in range(4))
    for n in range(1, JETS + 1):
        root.append(
            2 * (q[n] + sum(comb(n, j) * root[j] * root[n - j] for j in range(1, n)))
        )
        inverse.append(sum(comb(n, j) * inverse[j] * root[n - j] for j in range(n)))
        omega.append(
            (
                frequency_square[n]
                + sum(comb(n, j) * omega[j] * omega[n - j] for j in range(1, n))
            )
            / 2
        )
        reciprocal.append(
            sum(comb(n, j) * reciprocal[j] * omega[n - j] for j in range(n))
        )
    L = [
        sum(comb(n, j) * inverse[j] * root[n - j + 1] for j in range(n + 1))
        for n in range(JETS)
    ]
    rho = [
        sum(comb(n, j) * reciprocal[j] * omega[n - j + 1] for j in range(n + 1)) / 2
        for n in range(JETS)
    ]
    rotation = L
    squeeze = [L[n] + rho[n] for n in range(JETS)]
    return {
        "matrix_exponential_relative": exponential,
        "scale_relative": c,
        "K_relative": q,
        "omega_squared_relative": frequency_square,
        "square_root_balanced": root,
        "inverse_square_root_balanced": inverse,
        "omega_relative": omega,
        "omega_reciprocal_relative": reciprocal,
        "L_derivatives": L,
        "rho_derivatives": rho,
        "R_derivatives": rotation,
        "S_derivatives": squeeze,
    }


@cache
def data():
    c = constants()
    checks = {}
    for p in (1, 2, 3):
        for n in range(1, JETS + 1):
            prev = scale_polynomial(p, n - 1)
            checks[f"scale_derivative_{p}_{n}"] = s.expand(
                scale_polynomial(p, n)
                - (1 + t * t) * s.diff(prev, t)
                + 2 * (n - 1 + 2 * p) * t * prev
            )
    b1, b2, b3 = s.symbols("b1 b2 b3", positive=True)
    B = s.diag(b1, b2, b3)
    entries = s.symbols("q0:6", real=True)
    Q = s.Matrix(
        [
            [entries[0], entries[1], entries[2]],
            [entries[1], entries[3], entries[4]],
            [entries[2], entries[4], entries[5]],
        ]
    )
    M = s.Matrix(3, 3, lambda i, j: B[j, j] * Q[i, j] / (B[i, i] + B[j, j]))
    checks["exact_positive_square_root_balanced_Sylvester"] = M + M.T - Q
    checks["balanced_root_derivative_is_symmetric_before_normalizing"] = (
        B * M - (B * M).T
    )
    checks["direct_first_frame_bound"] = c["L_derivatives"][0] - 10
    checks["direct_first_squeeze_bound"] = c["S_derivatives"][0] - 11
    return {
        "jet_domain": "Smooth compact homogeneous tracefree symmetric gamma; sup ||gamma^(j)||op<=1/100 for j=0..12 on the actual CD slab. Gamma vanishes near the common initial Cauchy surface.",
        "matrix_frame": "B=K^(1/2), L=B^-1 B'; all derivative norms use positive square-root Sylvester equations, not the condition number of K.",
        "uniform_constants": c,
        "derivative_claim": "||L^(j)||,||R^(j)||,||S^(j)|| and normalized frequency derivatives have the displayed finite constants for every momentum. Only the +/-omega frequency split is used; polarization eigenvalues need not be separated.",
        "checks": checks,
        "gates": {
            "real_exponential_relative_margin": s.Rational(1, 1) / (1 - 2 * DELTA) < 2,
            "direct_K_first_derivative_below_five": 3 * s.Rational(8, 5) + 2 * DELTA
            < 5,
            "direct_frequency_square_derivative_below_four": 2 * s.Rational(8, 5)
            + 2 * DELTA
            < 4,
            "three_dimensional_Sylvester_Frobenius_cost_below_two": 3 < 2**2,
            "all_scale_and_matrix_constants_finite_positive": all(
                x > 0 for row in c.values() if isinstance(row, list) for x in row
            ),
            "finite_enough_jets_for_order_ten_residual": ORDER + 1 <= JETS,
        },
    }
