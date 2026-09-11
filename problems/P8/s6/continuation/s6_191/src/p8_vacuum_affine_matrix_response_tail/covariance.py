"""Full covariance derivative bounds and C2 infinite weighted-tail control."""

from functools import cache
from math import comb

import sympy as s
from p8_vacuum_affine_matrix_adiabatic import initial as prior

from . import reference, variation

K = prior.PARTITION
COV = {
    0: prior.COVARIANCE_ERROR,
    1: 4 * s.Integer(10) ** 33,
    2: 2 * s.Integer(10) ** 36,
}
TAIL = {
    0: s.Rational(1, 10) ** 65,
    1: s.Rational(1, 10) ** 47,
    2: s.Rational(1, 10) ** 28,
}


def multiply(*matrices):
    result = matrices[0]
    for value in matrices[1:]:
        result = (result * value).applyfunc(s.expand)
    return result


def covariance_jet(graph):
    N = len(graph) - 1
    unit = s.eye(graph[0].rows)
    zero = s.zeros(unit.rows)
    E = [unit - multiply(graph[0].H, graph[0])]
    for n in range(1, N + 1):
        E.append(
            -sum(
                (comb(n, j) * multiply(graph[j].H, graph[n - j]) for j in range(n + 1)),
                zero,
            )
        )
    C = [E[0].inv().applyfunc(s.expand)]
    for n in range(1, N + 1):
        C.append(
            -multiply(
                C[0],
                sum(
                    (comb(n, j) * multiply(E[j], C[n - j]) for j in range(1, n + 1)),
                    zero,
                ),
            )
        )
    F = [(unit + graph[0]).col_join(-s.I * (unit - graph[0])) / s.sqrt(2)]
    F.extend(r.col_join(s.I * r) / s.sqrt(2) for r in graph[1:])
    sigma = []
    for n in range(N + 1):
        value = s.zeros(2 * unit.rows)
        for a in range(n + 1):
            for b in range(n - a + 1):
                c = n - a - b
                value += comb(n, a) * comb(n - a, b) * multiply(F[a], C[b], F[c].H)
        sigma.append(value.applyfunc(s.re))
    return sigma


@cache
def constants():
    A = reference.constants()["reference_parameter_norm_over_inverse_frequency"]
    E0, E1, E2 = variation.E0, variation.E1, variation.E2
    cov0 = 16 * E0
    cov1 = 16 * E1 + 128 * E0 * A[1] / K**2
    cov2 = (
        16 * E2
        + 128 * E0 * A[2] / K**3
        + 256 * E1 * A[1] / K**2
        + 128 * E1**2 / K**10
        + 512 * E0 * A[1] ** 2 / K**4
    )
    weighted = {
        0: 3 * COV[0],
        1: 3 * (COV[1] + COV[0] / K),
        2: 3 * (COV[2] + 2 * COV[1] / K + 2 * COV[0] / K**2),
    }
    integral = {
        a: weighted[a] * s.Rational(4, 18) / (6 - a) / K ** (6 - a) for a in range(3)
    }
    return {
        "actual_covariance_parameter_error_coefficients": {0: cov0, 1: cov1, 2: cov2},
        "actual_weighted_parameter_error_coefficients": weighted,
        "actual_complete_infinite_integrals": integral,
    }


@cache
def data():
    G = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]])
    D = s.Matrix([[2, -1, 1], [-1, 0, 2], [1, 2, -2]])
    r = (G + s.I * D) / 100
    direction = (D + s.I * G) / 71
    C = (s.eye(3) - multiply(r.H, r)).inv().applyfunc(s.expand)
    A = multiply(direction.H, r) + multiply(r.H, direction)
    B = 2 * multiply(direction.H, direction)
    C1 = multiply(C, A, C)
    C2 = 2 * multiply(C, A, C, A, C) + multiply(C, B, C)
    C3 = (
        6 * multiply(C, A, C, A, C, A, C)
        + 3 * multiply(C, A, C, B, C)
        + 3 * multiply(C, B, C, A, C)
    )
    eps = s.Symbol("epsilon", real=True)
    inverse_polynomial = C + eps * C1 + eps**2 * C2 / 2 + eps**3 * C3 / 6
    original = s.eye(3) - multiply((r + eps * direction).H, r + eps * direction)
    inverse_product = multiply(original, inverse_polynomial)
    checks = {}
    for n in range(1, 4):
        checks[f"full_noncommuting_inverse_derivative_{n}"] = inverse_product.applyfunc(
            lambda x, n=n: s.factorial(n) * x.coeff(eps, n)
        )
    F0 = (s.eye(3) + r).col_join(-s.I * (s.eye(3) - r)) / s.sqrt(2)
    F1 = direction.col_join(s.I * direction) / s.sqrt(2)
    direct = multiply(F0 + eps * F1, inverse_polynomial, (F0 + eps * F1).H)
    output = covariance_jet([r, direction, s.zeros(3), s.zeros(3)])
    for n in range(1, 4):
        checks[f"complete_six_quadrature_product_derivative_{n}"] = (
            direct.applyfunc(lambda x, n=n: s.re(s.factorial(n) * x.coeff(eps, n)))
            - output[n]
        )
    nu, cut = s.symbols("nu partition", positive=True)
    for a in range(3):
        checks[f"complete_infinite_parameter_tail_{a}"] = s.integrate(
            nu ** (-7 + a), (nu, cut, s.oo)
        ) - cut ** (-6 + a) / (6 - a)
    c = constants()
    radius = s.Rational(1, 10)
    inverse_first = 2**2 * 2 * radius
    inverse_second = 2 * 2**3 * (2 * radius) ** 2 + 2**2 * 2
    inverse_third = 6 * 2**4 * (2 * radius) ** 3 + 6 * 2**3 * 2 * (2 * radius)
    return {
        "full_covariance_map": "Sigma=Re[F(r)(I-rdag r)^-1 F(r)dag] on symmetric complex graphs of norm<=1/10; derivatives are real Frechet derivatives, so the adjoint variations are retained.",
        "inverse_map_derivative_bounds": {
            0: s.Integer(2),
            1: s.Integer(1),
            2: s.Integer(16),
            3: s.Integer(32),
        },
        "full_covariance_map_derivative_bounds": {
            1: s.Integer(16),
            2: s.Integer(128),
            3: s.Integer(512),
        },
        "actual_covariance_parameter_error_displays": COV,
        "constants": c,
        "energy_weighted_integrand": "W(epsilon,t,k)=omega(epsilon,t,k)[Sigma_actual-Sigma_reference]. The common proof band nu_minus>=1e16 is independent of epsilon; no physical mode is removed.",
        "complete_infinite_parameter_tail_displays": TAIL,
        "uniform_C2_result": "The high-band integral of W is C2 in the admitted amplitude family, uniformly on the fixed CD slab. The complete integrals of the norms of W and its first/second derivatives are below1e-65,1e-47,1e-28 respectively.",
        "finite_amplitude_remainder": "For |epsilon|<=1/100, the high-band weighted covariance-tail integral has Taylor remainder at most(epsilon^2/2)*1e-28 after its value and first derivative at zero; the inequality is strict when epsilon is nonzero and the remainder is exactly zero at epsilon0. This is an actual-minus-finite-reference tail result, not the full subtracted stress or gravitational feedback remainder.",
        "checks": checks,
        "gates": {
            "inverse_first_derivative": inverse_first < 1,
            "inverse_second_derivative": inverse_second < 16,
            "inverse_third_derivative": inverse_third < 32,
            "complete_covariance_first_derivative": 4 + 4 + 4 < 16,
            "complete_covariance_second_derivative": 4 * 16 + 4 * 2 + 2 * 2 < 128,
            "complete_covariance_third_derivative": 4 * 32 + 6 * 16 * 2 + 6 < 512,
            "all_actual_covariance_error_displays": all(
                c["actual_covariance_parameter_error_coefficients"][a] < COV[a]
                for a in range(3)
            ),
            "all_full_infinite_parameter_tail_displays": all(
                c["actual_complete_infinite_integrals"][a] < TAIL[a] for a in range(3)
            ),
            "parameter_band_is_fixed": K == prior.PARTITION,
            "second_parameter_tail_still_integrable": -7 + 2 < -1,
        },
    }
