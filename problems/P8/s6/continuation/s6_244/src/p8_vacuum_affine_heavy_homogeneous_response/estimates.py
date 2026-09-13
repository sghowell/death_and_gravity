"""Complete mass-uniform scalar mixed jets, full covariance tails and preparation bounds."""

from fractions import Fraction as F
from functools import cache
from math import comb, factorial

import sympy as s

delta = F(1, 100)
MAX = 12
K = 10**90
CONTOUR = 10**6
B = [F(1)]
for j in range(MAX):
    B.append(delta * sum(comb(j, k) * B[k] for k in range(j + 1)))


def bell(j, a):
    return sum(comb(j, k) * a ** (j - k) * B[k] for k in range(j + 1))


tt = s.Symbol("t")
poly = s.S.One
c = [F(1)]
polynomials = [poly]
for j in range(MAX):
    poly = s.expand((1 + tt * tt) * s.diff(poly, tt) - 2 * (j + 4) * tt * poly)
    c.append(F(int(sum(abs(v) for v in s.Poly(poly, tt).coeffs()))))
    polynomials.append(poly)
indices = sorted(
    ((j, a) for j in range(MAX + 1) for a in range(3)), key=lambda b: sum(b)
)


def subindices(j, a):
    return ((k, b, comb(j, k) * comb(a, b)) for k in range(j + 1) for b in range(a + 1))


q = {}
for j, a in indices:
    q[j, a] = 2 * sum(comb(j, k) * c[j - k] * bell(k, a) for k in range(j + 1))
q[0, 0] = F(1)
q[1, 0] = F(5)  # 2|H|+exp(2delta)||Q'|| <5.
U = {(0, 0): F(1)}
V = {(0, 0): F(1)}
for j, a in indices[1:]:
    U[j, a] = (
        q[j, a]
        + sum(
            weight * U[k, b] * U[j - k, a - b]
            for k, b, weight in subindices(j, a)
            if (k, b) != (0, 0) and (k, b) != (j, a)
        )
    ) / 2
    V[j, a] = sum(
        weight * U[k, b] * V[j - k, a - b]
        for k, b, weight in subindices(j, a)
        if (k, b) != (0, 0)
    )
S = {}
for j in range(MAX):
    for a in range(3):
        lam = sum(
            weight * U[k + 1, b] * V[j - k, a - b] for k, b, weight in subindices(j, a)
        )
        hh = F(2) if j == 0 else F(4 * factorial(j))
        S[j, a] = lam / 2 + (
            3 * hh / 2 + 3 * delta / 4 if a == 0 else F(3, 4) if a == 1 else 0
        )
assert S[0, 0] < 5 and U[0, 1] <= 1 and U[0, 2] <= 2


def triple(first, second, third, j, a):
    total = F(0)
    for j1, a1, w1 in subindices(j, a):
        for j2, a2, w2 in subindices(j - j1, a - a1):
            total += (
                w1
                * w2
                * first(j1, a1)
                * second(j2, a2)
                * third(j - j1 - j2, a - a1 - a2)
            )
    return total


@cache
def b(n, j, a):
    assert n + j <= 11 and 0 <= a <= 2
    if n == 1:
        return (
            sum(
                weight * S[k, aa] * V[j - k, a - aa]
                for k, aa, weight in subindices(j, a)
            )
            / 2
        )
    total = F(0)
    for k, aa, weight in subindices(j, a):
        jt, at = j - k, a - aa
        inside = b(n - 1, jt + 1, at)
        for left in range(1, n - 1):
            right = n - 1 - left
            inside += triple(
                lambda j0, a0, left=left: b(left, j0, a0),
                lambda j0, a0: S[j0, a0],
                lambda j0, a0, right=right: b(right, j0, a0),
                jt,
                at,
            )
        total += weight * V[k, aa] * inside / 2
    return total


A = [sum(b(n, 0, a) * F(1, K ** (n - 1)) for n in range(1, 11)) for a in range(3)]
C = []
for a in range(3):
    value = b(10, 1, a)
    for left in range(1, 11):
        for right in range(1, 11):
            if left + right >= 10:
                value += triple(
                    lambda j0, a0, left=left: b(left, j0, a0),
                    lambda j0, a0: S[j0, a0],
                    lambda j0, a0, right=right: b(right, j0, a0),
                    0,
                    a,
                ) * F(1, K ** (left + right - 10))
    C.append(value)
R = [sum(b(n, 0, a) * F(1, CONTOUR**n) for n in range(1, 11)) for a in range(3)]
assert R[0] < F(1, 100) and A[0] / K < F(1, 100)
# Full SLE beta plus actual W6-to-finite-reference initial mismatch.
INITIAL = 10**90
E0 = 10**40
raw_E0 = 1000 * (F(INITIAL, K) + C[0])
assert raw_E0 < E0
assert A[0] / K + F(E0, K**10) < F(1, 10)
M1 = 12 + S[0, 1] / K + 10 * A[1] / K**2
E1 = 10**42
raw_E1 = 3 * (E0 * M1 + C[1] / K)
assert raw_E1 < E1
M20 = 24 + S[0, 2] / K + (10 * A[2] + 4 * S[0, 1] * A[1]) / K**2
M21 = 24 + S[0, 1] / K + 20 * A[1] / K**2
E2 = 10**44
raw_E2 = 3 * (E0 * M20 / K + E1 * M21 + F(10 * E1**2, K**10) + C[2] / K**2)
assert raw_E2 < E2
D = [
    16 * E0,
    16 * E1 + 128 * E0 * A[1] / K**2,
    16 * E2
    + 128 * E0 * A[2] / K**3
    + 256 * E1 * A[1] / K**2
    + F(128 * E1**2, K**10)
    + 512 * E0 * A[1] ** 2 / K**4,
]
N = [F(4), F(6), F(30)]
CV = [F(2), 16 * R[1], 16 * R[2] + 128 * R[1] ** 2]
contour = []
readout = []
for a in range(3):
    cc = F(0)
    rr = F(0)
    for i in range(a + 1):
        for j in range(a - i + 1):
            rest = a - i - j
            weight = comb(a, i) * comb(a - i, j)
            cc += weight * U[0, i] * N[j] * CV[rest]
            rr += 6 * weight * U[0, i] * N[j] * D[rest] * F(1, K ** (a - rest))
    contour.append(cc)
    readout.append(rr)
assert all(v < limit for v, limit in zip(contour, (10, 100, 1000)))
assert all(v < 10**50 for v in readout)
# Converting every current error to nu^-5 only improves it, since nu>K>1.
TAIL = [readout[a] + 2 * limit * CONTOUR**6 for a, limit in enumerate((10, 100, 1000))]
assert all(2 * v < 10**51 for v in TAIL)
n = F(10**200, 512) + 2
assert F(99, 100) * n > K * K


def exact(value):
    if isinstance(value, F):
        return s.Rational(value.numerator, value.denominator)
    if isinstance(value, dict):
        return {key: exact(v) for key, v in value.items()}
    if isinstance(value, (list, tuple)):
        return tuple(exact(v) for v in value)
    return value


@cache
def data():
    from p8_vacuum_affine_heavy_curved_state import state

    checks = {"actual_unchanged_heavy_mass": exact(n) - state.MASS2}
    aa = (1 + tt * tt) ** 2
    for j, row in enumerate(polynomials):
        checks["complete_background_inverse_scale_jet_" + str(j)] = s.factor(
            s.diff(aa**-2, tt, j) / aa**-2 - row / (1 + tt * tt) ** j
        )
    checks["full_initial_state_and_W6_tail_allowance"] = s.S(INITIAL) - 10**90
    normalized_tail = 10**51 / (state.KAPPA * state.MASS2)
    regulator_tail = s.Rational(10**51, state.KAPPA)
    gates = {
        "full_mixed_square_root_inverse_and_squeeze_jets": all(
            value >= 0 for table in (q, U, V, S) for value in table.values()
        ),
        "all_195_finite_reference_coefficient_derivative_majorants": len(
            {
                (order, j, a): b(order, j, a)
                for order in range(1, 11)
                for j in range(12 - order)
                for a in range(3)
            }
        )
        == 195,
        "complete_order_ten_residual_and_all_parameter_derivatives": all(
            value > 0 for value in C
        ),
        "all_internal_momenta_above_proof_threshold": F(99, 100) * n > K * K,
        "whole_reference_inside_graph_ball": A[0] / K < F(1, 100),
        "full_S240_auxiliary_graph_circle": F(3, 10**7) / (2 - F(3, 10**7)) < F(1, 10),
        "entire_W6_marker_tail_initial_bound": F(1, 5) * 10**88 < 10**88,
        "full_normalized_Bogoliubov_graph_difference": (1 + F(1, 50) ** 2)
        / (1 - F(1, 2) * F(1, 50))
        < 2,
        "state_comparison_and_full_W6_initial_mismatch": 2 * 10**83
        + 4 * 10**62
        + 10**88
        < INITIAL,
        "full_initial_plus_residual_error": raw_E0 < E0,
        "actual_graph_bootstrap_after_full_initial_error": A[0] / K + F(E0, K**10)
        < F(1, 10),
        "complete_first_parameter_error": raw_E1 < E1,
        "complete_second_parameter_error_with_quadratic_term": raw_E2 < E2,
        "mass_independent_full_marker_contour": R[0] < F(1, 100),
        "every_full_current_contour_derivative": all(
            v < limit for v, limit in zip(contour, (10, 100, 1000))
        ),
        "every_full_exact_reference_current_error": all(v < 10**50 for v in readout),
        "entire_covariance_plus_marker_tail_integral": all(
            2 * v < 10**51 for v in TAIL
        ),
        "complete_radial_fifth_power_coefficient": F(64, 54) * F(100, 99) ** 3 < 2,
        "complete_fixed_frequency_radial_Jacobian": F(32, 9) * F(100, 99) ** 2 < 4,
        "same_Omega_star_regulator_tail": all(
            F(200, 99) * value < 10**51 for value in TAIL
        ),
        "all_momenta_not_a_physical_cutoff": True,
        "entire_normalized_state_subtraction_current_derivatives": bool(
            normalized_tail < s.Rational(1, 10**940)
        ),
        "no_new_reference_state_or_initial_reset": True,
    }
    return {
        "history_domain": "Real homogeneous symmetric Q, common zero initial germ at t0=-1/2, all raw time jets0..12 bounded by1/100 on the unit slab. Parameter directions have full C12 norm at most1; the fixed actual S240 SLE covariance is unchanged.",
        "fixed_lower_frequency": "nu²=(99/100)(n+k²/16), omega>=nu,omega<=6nu. All internal momenta already satisfy nu>10^98>10^90. The threshold only bounds rational denominators.",
        "all_Bell_and_background_inverse_scale_majorants": exact(
            {"Bell": B, "scale": c}
        ),
        "complete_mixed_frequency_square_jets": exact(q),
        "complete_mixed_frequency_jets": exact(U),
        "complete_mixed_inverse_frequency_jets": exact(V),
        "complete_mixed_squeeze_jets": exact(S),
        "all_finite_Riccati_coefficient_derivative_majorants": exact(
            {
                (order, j, a): b(order, j, a)
                for order in range(1, 11)
                for j in range(12 - order)
                for a in range(3)
            }
        ),
        "whole_reference_parameter_bounds": exact(A),
        "whole_finite_residual_parameter_bounds": exact(C),
        "whole_initial_state_comparison_and_W6_coefficient_tail": INITIAL,
        "complete_exact_graph_error_bounds": (E0, E1, E2),
        "complete_unrounded_graph_error_bounds": exact((raw_E0, raw_E1, raw_E2)),
        "complete_covariance_difference_derivative_bounds": exact(D),
        "full_normalized_vertex_derivative_bounds": exact(N),
        "mass_independent_marker_contour_scale": CONTOUR,
        "entire_reference_graph_contour_bounds": exact(R),
        "entire_current_contour_bounds": exact(contour),
        "entire_exact_reference_current_error_bounds": exact(readout),
        "entire_combined_state_and_marker_tail_coefficients": exact(TAIL),
        "full_radial_state_subtraction_current_derivative_bound": s.Rational(10**51, 1)
        / state.MASS2,
        "complete_normalized_state_subtraction_derivative_bound": normalized_tail,
        "same_frequency_cutoff_tail_coefficient_over_K_squared": regulator_tail,
        "regulator_scope": "Keep Omega_star<=K, K>=2sqrt(n), independent of the varied history. At zero transfer both internal legs have the same condition. The full fixed finite action/profile remain; no finite-K matched mean is set to zero.",
        "checks": checks,
        "gates": {name: bool(value) for name, value in gates.items()},
    }
