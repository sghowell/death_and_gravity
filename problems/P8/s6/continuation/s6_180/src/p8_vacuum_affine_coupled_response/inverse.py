"""Constructive retained Volterra inverse and quantitative first-row recovery."""

from functools import cache

import sympy as s
from p8_vector_state import wkb

from . import classical


@cache
def first_row():
    tree = classical.adapted()
    A = tree["A"]
    coefficients = tree["complete_derivative_coefficients"]
    checks = {}

    def bound(expr, name):
        b = wkb.box_bound(s.factor(expr.subs(classical.u, wkb.u)))
        checks[name] = b["reconstruction"]
        return b["absolute_upper"]

    Aj = {
        j: bound(s.diff(A, classical.u, j), "A_derivative_" + str(j))
        for j in range(1, 5)
    }
    rows = {}
    for col in (0, 1):
        for order in range(4):
            coeff = coefficients["0" + str(col)][order]
            rows[(col, order)] = {
                j: bound(
                    s.diff(coeff, classical.u, j),
                    f"coefficient_{col}_{order}_derivative_{j}",
                )
                for j in range(order + 1)
            }
    # K_j = sum_k (-1)^k binom(j,k) lag^(3-j+k)c_j^(k)/(3-j+k)!.
    # Separate the A*eta endpoint from I4(A*eta'''') exactly.
    diagonal = 4 * Aj[1] + sum(rows[(col, 3)][0] for col in (0, 1))
    derivative = sum(s.binomial(4, j) * Aj[j] / s.factorial(j - 2) for j in range(2, 5))
    for row in rows.values():
        order = len(row) - 1
        for j, b in row.items():
            degree = 3 - order + j
            if degree >= 1:
                derivative += s.binomial(order, j) * b / s.factorial(degree - 1)
    # gamma times scaled P^(j) = physical P^(j)/kappa < 1e-770,
    # j=0,1 by S6.176. The two actual quantum first-row coefficients
    # are -3H' P/kappa (order0) and -3P/kappa (order1).
    quantum = s.Rational(21, 2) * s.Rational(1, 10**770)
    C1 = Aj[1] + diagonal + derivative + quantum
    return {
        "actual_classical_first_row_coefficient_derivative_envelopes": rows,
        "actual_A_derivative_envelopes": Aj,
        "kernel_diagonal_row_upper": diagonal,
        "kernel_first_output_derivative_integral_upper": derivative + quantum,
        "quantum_first_row_added_bound": quantum,
        "physical_lapse_reconstruction_C1_upper": C1,
        "C1_bound_uses": "The unit slab, |H'|<=4, the exact current coefficients and the fixed S6.176 physical stress/first derivative bounds. It is a numeric bound for the local first-row reconstruction, not for the full inverse.",
        "checks": checks,
        "gates": {
            "C1_positive_finite": C1 > 0,
            "C1_below_15000": C1 < 15000,
            "quantum_first_row_small": quantum < s.Rational(1, 10**768),
        },
    }


@cache
def data():
    t, r = s.symbols("t r", positive=True)
    A = s.Function("A")
    coeff = s.Function("coefficient")
    f = (t - r) ** 3 * A(r) / 6
    remainder = s.diff(f, r, 4)
    explicit = sum(
        (-1) ** j
        * s.binomial(4, j)
        * (t - r) ** (j - 1)
        * s.diff(A(r), r, j)
        / s.factorial(j - 1)
        for j in range(1, 5)
    )
    checks = {
        "variable_fourth_coefficient_primitive_endpoint": s.diff(f, r, 3).subs(r, t)
        + A(t),
        "variable_fourth_coefficient_full_remainder": s.expand(remainder - explicit),
    }
    kernels = {}
    for order in range(4):
        actual = (-1) ** order * s.diff((t - r) ** 3 * coeff(r) / 6, r, order)
        expected = sum(
            (-1) ** j
            * s.binomial(order, j)
            * (t - r) ** (3 - order + j)
            * s.diff(coeff(r), r, j)
            / s.factorial(3 - order + j)
            for j in range(order + 1)
        )
        kernels[order] = expected
        checks["local_primitive_kernel_" + str(order)] = s.expand(actual - expected)
        checks["local_primitive_output_derivative_" + str(order)] = s.expand(
            s.diff(actual - expected, t)
        )
    C, K, gamma = s.symbols("C K_L1 gamma", positive=True)
    invA = s.Rational(15625, 6144)
    beta = C * (invA + K / gamma)
    weight = (4 * beta + 1) ** 2
    # Integral e^-lambda*r(1-log r) <= (2+log lambda)/lambda
    # and log lambda <= 2sqrt(lambda)-2 for lambda>=1.
    contraction = 2 * beta / (4 * beta + 1)
    checks["constructive_weight_positive_gap"] = s.factor(
        s.Rational(1, 2) - contraction - 1 / (2 * (4 * beta + 1))
    )
    checks["simple_parameterized_weight_square"] = weight - (4 * beta + 1) ** 2
    Hmax = s.Rational(8, 5)
    checks["original_force_I4_bound"] = (
        s.Rational(1, 6) + 4 * Hmax / 24 - s.Rational(13, 30)
    )
    checks["original_force_I3_bound"] = (
        s.Rational(1, 2) + 4 * Hmax / 6 - s.Rational(47, 30)
    )
    checks["physical_scale_reconstruction"] = 1 + Hmax - s.Rational(13, 5)
    x = s.symbols("positive_x", positive=True)
    checks["elementary_log_bound_derivative"] = s.factor(
        s.diff(x - 1 - s.log(x), x) - (x - 1) / x
    )
    return {
        "normal_form": "I4(T+gamma Q)=diag(A(t),gamma F_m(partial_t^2))+V, A=-6delta(t)^2, gamma=1/(64pi^2*kappa). All coefficient multipliers remain inside their derivatives.",
        "diagonal_inverse": "B0=diag(A^-1 multiplication,gamma^-1 K_m convolution). The rank-one loop block alone is not inverted as a two-source identity.",
        "weak_log_majorant_parameters": {
            "C": C,
            "K_L1": K,
            "gamma": gamma,
            "A_inverse_upper": invA,
        },
        "integrable_composition_constant": beta,
        "explicit_exponential_weight": weight,
        "contraction_upper": contraction,
        "adapted_inverse_C0_upper_formula": 2 * s.exp(weight) * s.Max(invA, K / gamma),
        "actual_local_primitive_kernels": kernels,
        "variable_fourth_coefficient_remainder": explicit,
        "original_force_I4_upper": s.Rational(13, 30),
        "original_force_I3_upper": s.Rational(47, 30),
        "physical_scale_over_adapted_upper": s.Rational(13, 5),
        "existence": "For each smooth prepared original two-force there is a unique smooth prepared solution of the retained current-parent homogeneous tree-plus-conditional-vector RESPONSE equations on I. The same initial covariance and zero matter-charge perturbation are retained. The original lapse/scale have finite C0 dependence on the original forcing.",
        "uncomputed": "C and the full K_m L1 norm are finite by the matched continuum and scalar-kernel proofs but are not numerically evaluated here. The explicit local C1 bound does not make the full inverse small.",
        "nonstationary_boundary": "This is the derivative of the equations at the uncorrected CD reference, whose conditional vector one-point stress is nonzero. Its residual is not a prepared forcing near the original Cauchy slice. The inverse is not applied to that residual without a compatible preparation/background construction.",
        "checks": checks,
        "gates": {
            "contraction_strictly_below_half": s.factor(
                s.Rational(1, 2) - contraction
            ).is_positive
            is True
        },
    }
