"""Independent exact Legendre, finite Fock, phase and coefficient diagnostics.

Finite occupation diagnostics are not proofs of continuum Wick extension or
physical loop bounds. Quartic vacuum moments below are exactly supported in the
chosen occupation space; evolved-state tests separately compare two cutoffs.
"""

from functools import cache

import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import phase
from p8_vacuum_affine_gaussian_measure_response import (
    audit,
    functional,
    measure,
    response,
)
from p8_vacuum_affine_scalar_tame_propagator import charts
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm
from scipy.sparse import csr_matrix, diags, eye, kron
from scipy.sparse.linalg import expm_multiply


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_every_complete_exact_identity(name, value):
    assert all(
        entry == 0
        for entry in (list(value) if isinstance(value, s.MatrixBase) else [value])
    ), name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_every_scoped_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input_or_scope_promotion(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def independent_symplectic(case):
    r = s.Rational(case + 1, 13)
    linear = s.Matrix([[1 + r, r / 7], [r / 11, 1 + r / 5]])
    lower = s.eye(4)
    upper = s.eye(4)
    lower[2:, :2] = s.Matrix([[r / 17, -r / 19], [-r / 19, r / 23]])
    upper[:2, 2:] = s.Matrix([[r / 29, r / 31], [r / 31, -r / 37]])
    return lower * s.diag(linear, linear.inv().T) * upper


def independent_vertex(case, dimension=4):
    return s.Matrix(
        dimension,
        dimension,
        lambda i, j: (
            s.Rational(
                (-1) ** (i + j) * (case + 1 + i + 2 * j + 2 * i + j),
                11 + i + j + i * j,
            )
            + (s.Rational(case + i + 2, 3) if i == j else 0)
        ),
    )


@pytest.mark.parametrize("case", range(8))
def test_separately_entered_full_lapse_shift_action_and_Dirac_measure(case):
    v, sig, n, b = measure.POSITIONS
    pv, ps, pn, pb = measure.MOMENTA
    vd, sd = measure.vd, measure.sd
    J = s.Rational(case + 3, 5)
    theta = s.Rational(case - 3, 7)
    ell = s.Rational(case + 1, 17)
    E = s.Rational(2 * case - 5, 9)
    A = -s.Rational(case + 1, 101)
    T = s.Rational(case + 2, 103)
    q = s.Rational(case + 4, 3)
    density = s.Integer(10) ** 800 * s.Rational(case + 5, 4) ** 3
    w = -ell * E
    L = density * (
        -3 * vd**2
        + (J + w**2 / 2 - 3 * theta**2) * n**2
        + 6 * theta * n * vd
        + sd**2 / 2
        + w * n * sd
        - 3 * ell * vd * sig
        + 2 * b * (vd - theta * n)
        + ell * b * sig
        + q * v**2
        + 2 * E * q * n * v
        - q * sig**2 / 2
        + 3 * T * n * v
        + s.Rational(9, 2) * A * v**2
    )
    substitution = {
        charts.J: J,
        charts.th: theta,
        charts.l: ell,
        charts.E: E,
        charts.A: A,
        charts.T: T,
        charts.q: q,
        measure.W: density,
    }
    assert s.expand(measure.complete_lagrangian().subs(substitution) - L) == 0
    rates = {
        vd: -(pv / density - 6 * theta * n + 3 * ell * sig - 2 * b) / 6,
        sd: ps / density - w * n,
    }
    H = s.expand((pv * vd + ps * sd - L).subs(rates, simultaneous=True))
    assert s.expand(H - measure.complete_hamiltonian().subs(substitution)) == 0
    constraint = s.Matrix([pn, pb, s.diff(H, n), s.diff(H, b)])
    coordinate = [v, sig, n, b]
    momentum = [pv, ps, pn, pb]
    bracket = lambda f, g: s.expand(
        sum(
            s.diff(f, x) * s.diff(g, p) - s.diff(f, p) * s.diff(g, x)
            for x, p in zip(coordinate, momentum)
        )
    )
    matrix = s.Matrix(4, 4, lambda i, j: bracket(constraint[i], constraint[j]))
    D = s.hessian(H, [n, b])
    assert D.det() == 4 * density**2 * J / 3
    assert matrix.det() == D.det() ** 2
    assert matrix[2, 3] == -density * (2 * E * q + 3 * T) / 3
    inverse = matrix.inv()
    physical = [v, sig, pv, ps]
    for i in range(4):
        for j in range(4):
            corrected = bracket(physical[i], physical[j]) - sum(
                bracket(physical[i], constraint[a])
                * inverse[a, z]
                * bracket(constraint[z], physical[j])
                for a in range(4)
                for z in range(4)
            )
            assert s.cancel(corrected - phase.OMEGA[i, j]) == 0
    configuration_auxiliary = s.hessian(L, [n, b])
    assert configuration_auxiliary.det() == -4 * density**2 * theta**2
    if theta == 0:
        assert configuration_auxiliary.det() == 0 and matrix.det() != 0
    else:
        solution = s.solve([s.diff(L, n), s.diff(L, b)], [n, b])
        reduced = L.subs(solution, simultaneous=True)
        assert s.factor(s.hessian(reduced, [vd, sd]).det()) == (
            2 * density**2 * J / theta**2
        )


@pytest.mark.parametrize("case", range(6))
def test_independent_sourced_finite_Gaussian_Schur_integration(case):
    m = s.Integer(case + 1)
    full = s.Matrix(
        [
            [12 + m, 1, 2, -1],
            [1, 15 + m, 1, 3],
            [2, 1, 19 + m, 2],
            [-1, 3, 2, 23 + m],
        ]
    )
    source = s.Matrix([m / 3, -m / 7, (m + 1) / 11, -1])
    A, C, D = full[:2, :2], full[:2, 2:], full[2:, 2:]
    u, v = source[:2, :], source[2:, :]
    schur = A - C * D.inv() * C.T
    corrected_source = u - C * D.inv() * v
    direct_exponent = (source.T * full.inv() * source)[0]
    sequential_exponent = (v.T * D.inv() * v)[0] + (
        corrected_source.T * schur.inv() * corrected_source
    )[0]
    assert s.cancel(direct_exponent - sequential_exponent) == 0
    assert full.det() == D.det() * schur.det()
    assert full.det() > 0 and D.det() > 0 and schur.det() > 0
    assert s.cancel((u.T * schur.inv() * u)[0] - direct_exponent) != 0


@cache
def fock_phase(cutoff, modes):
    annihilate = diags(
        np.sqrt(np.arange(1, cutoff)), 1, shape=(cutoff, cutoff), format="csr"
    )
    one = eye(cutoff, format="csr")
    all_a = []
    for mode in range(modes):
        product = csr_matrix([[1.0]])
        for slot in range(modes):
            product = kron(product, annihilate if slot == mode else one, format="csr")
        all_a.append(product)
    q = [(a + a.T.conjugate()) / np.sqrt(2) for a in all_a]
    p = [(a - a.T.conjugate()) / (1j * np.sqrt(2)) for a in all_a]
    return tuple(q + p)


def transformed_operators(matrix, operators):
    M = np.asarray(matrix, dtype=complex)
    return tuple(
        sum(
            (M[i, j] * z for j, z in enumerate(operators)),
            csr_matrix(operators[0].shape),
        )
        for i in range(len(operators))
    )


def quadratic_operator(matrix, operators):
    A = np.asarray(matrix, dtype=complex)
    out = csr_matrix(operators[0].shape, dtype=complex)
    for i, left in enumerate(operators):
        for j, right in enumerate(operators):
            out = out + A[i, j] * (left @ right + right @ left) / 4
    return out


@pytest.mark.parametrize("case", range(10))
def test_independent_finite_Fock_entire_coupled_Wick_mean_noise_and_commutator(case):
    S0, L, R = [independent_symplectic(case + offset) for offset in (0, 2, 5)]
    A, B, C = [independent_vertex(case + j) for j in range(3)]
    result = response.kernels(A, B, L, R, S0 * S0.T / 2, C)
    z = fock_phase(5, 2)
    at_left, at_right = (
        transformed_operators(L * S0, z),
        transformed_operators(R * S0, z),
    )
    OA, OB, OC = (
        quadratic_operator(A, at_left),
        quadratic_operator(B, at_right),
        quadratic_operator(C, at_left),
    )
    meanA, meanB = OA[0, 0], OB[0, 0]
    connected = (OA @ OB)[0, 0] - meanA * meanB
    np.testing.assert_allclose(complex(result["mean_A"]), meanA, rtol=2e-12, atol=2e-12)
    np.testing.assert_allclose(complex(result["mean_B"]), meanB, rtol=2e-12, atol=2e-12)
    np.testing.assert_allclose(
        complex(result["ordered_connected_Wick"]), connected, rtol=2e-12, atol=2e-12
    )
    np.testing.assert_allclose(
        float(result["symmetric_noise"]), connected.real, rtol=2e-12, atol=2e-12
    )
    susceptibility = -1j * ((OA @ OB - OB @ OA)[0, 0])
    np.testing.assert_allclose(
        float(result["observable_retarded_susceptibility_before_step"]),
        susceptibility,
        rtol=2e-12,
        atol=2e-12,
    )
    np.testing.assert_allclose(
        float(result["effective_action_seagull"]), -OC[0, 0], rtol=2e-12, atol=2e-12
    )
    centered = [OA - meanA * eye(25), OB - meanB * eye(25), OC - OC[0, 0] * eye(25)]
    gram = np.array([[(left @ right)[0, 0] for right in centered] for left in centered])
    assert np.linalg.eigvalsh(gram).min() > -2e-10
    weights = np.array([1.0, -0.7, 0.23])
    combined = sum(
        (weight * op for weight, op in zip(weights, centered)), csr_matrix((25, 25))
    )
    np.testing.assert_allclose(
        weights @ gram.real @ weights,
        (combined @ combined)[0, 0].real,
        rtol=2e-12,
        atol=2e-12,
    )
    # The vacuum quartic contractions only visit occupations <=2 per mode.
    # N=5 therefore introduces no occupation truncation in these moments.


@pytest.mark.parametrize("case", range(6))
def test_independent_TT_Fock_response_and_two_polarizations(case):
    r = s.Rational(case + 1, 11)
    S0 = s.Matrix([[1, r], [r / 7, 1 + r**2 / 7]])
    L = s.Matrix([[1 + r, 0], [r / 13, 1 / (1 + r)]])
    R = s.Matrix([[1, r / 3], [0, 1]])
    A, B = independent_vertex(case, 2), independent_vertex(case + 2, 2)
    result = response.kernels(A, B, L, R, S0 * S0.T / 2)
    operators = fock_phase(5, 1)
    OA = quadratic_operator(A, transformed_operators(L * S0, operators))
    OB = quadratic_operator(B, transformed_operators(R * S0, operators))
    connected = (OA @ OB)[0, 0] - OA[0, 0] * OB[0, 0]
    np.testing.assert_allclose(
        complex(result["ordered_connected_Wick"]), connected, rtol=2e-12, atol=2e-12
    )
    OA2 = kron(OA, eye(5)) + kron(eye(5), OA)
    OB2 = kron(OB, eye(5)) + kron(eye(5), OB)
    doubled = (OA2 @ OB2)[0, 0] - OA2[0, 0] * OB2[0, 0]
    np.testing.assert_allclose(doubled, 2 * connected, rtol=2e-12, atol=2e-12)


def finite_CTP_overlap(case, cutoff):
    operators = fock_phase(cutoff, 2)
    O = np.array(phase.OMEGA, dtype=float)
    H0 = np.array(independent_vertex(case), dtype=float)
    H1 = np.array(independent_vertex(case + 1), dtype=float)
    H2 = np.array(independent_vertex(case + 2), dtype=float)
    H3 = np.array(independent_vertex(case + 3), dtype=float)
    # Different full mixed vertices, not commuting occupation Hamiltonians.
    H2[0, 2] += 0.37
    H2[2, 0] += 0.37
    H3[1, 3] -= 0.29
    H3[3, 1] -= 0.29
    assert np.linalg.norm((O @ H1) @ (O @ H2) - (O @ H2) @ (O @ H1)) > 0.1
    duration = np.array([0.031, 0.043, 0.027, 0.038, 0.052]) / (case + 1)
    seed = np.zeros(cutoff**2, dtype=complex)
    seed[0] = 1
    prepared = expm_multiply(
        -1j * duration[0] * quadratic_operator(H0, operators), seed
    )
    plus = expm_multiply(
        -1j * duration[1] * quadratic_operator(H1, operators), prepared
    )
    plus = expm_multiply(-1j * duration[2] * quadratic_operator(H2, operators), plus)
    minus = expm_multiply(
        -1j * duration[3] * quadratic_operator(H2, operators), prepared
    )
    minus = expm_multiply(-1j * duration[4] * quadratic_operator(H3, operators), minus)
    S0 = expm(duration[0] * O @ H0)
    Splus = expm(duration[2] * O @ H2) @ expm(duration[1] * O @ H1)
    Sminus = expm(duration[4] * O @ H3) @ expm(duration[3] * O @ H2)
    relative = np.linalg.solve(S0, np.linalg.solve(Sminus, Splus @ S0))
    alpha = (
        relative[:2, :2] + relative[2:, 2:] + 1j * (relative[2:, :2] - relative[:2, 2:])
    ) / 2
    beta = (
        relative[:2, :2] - relative[2:, 2:] + 1j * (relative[2:, :2] + relative[:2, 2:])
    ) / 2
    np.testing.assert_allclose(
        alpha @ alpha.conjugate().T - beta @ beta.conjugate().T, np.eye(2), atol=2e-12
    )
    # These short paths remain in the identity branch. The full-turn test below
    # separately prevents treating this endpoint branch as a global formula.
    formula = np.exp(-np.log(np.linalg.det(alpha.conjugate())) / 2)
    return (
        np.vdot(minus, plus),
        formula,
        np.linalg.det(alpha @ alpha.conjugate().T).real ** -0.25,
    )


@pytest.mark.parametrize("case", range(5))
def test_independent_two_cutoff_Fock_evolution_full_CTP_overlap_and_modulus(case):
    small, formula, modulus = finite_CTP_overlap(case, 10)
    large, second_formula, second_modulus = finite_CTP_overlap(case, 14)
    np.testing.assert_allclose(small, large, atol=2e-11, rtol=2e-11)
    np.testing.assert_allclose(large, formula, atol=2e-11, rtol=2e-11)
    np.testing.assert_allclose(formula, second_formula, atol=2e-13, rtol=2e-13)
    np.testing.assert_allclose(abs(formula), modulus, atol=2e-13, rtol=2e-13)
    assert second_modulus <= 1 + 1e-13
    assert abs(formula.imag) > 1e-3
    assert abs(large - formula.conjugate()) > 1e-3


@pytest.mark.parametrize("turns", (1, 2, 3, 4, 5))
def test_closed_classical_path_requires_continued_metaplectic_phase(turns):
    angle = 2 * s.pi * turns
    rotation = s.Matrix([[s.cos(angle), s.sin(angle)], [-s.sin(angle), s.cos(angle)]])
    alpha, beta = functional.bogoliubov(rotation)
    assert alpha == s.eye(1) and beta == s.zeros(1)
    # Direct oscillator spectrum H|n>=(n+1/2)|n>, not endpoint det(alpha).
    evolved_vacuum_phase = s.exp(-s.I * angle / 2)
    assert evolved_vacuum_phase == (-1) ** turns
    if turns % 2:
        assert evolved_vacuum_phase != 1 / s.sqrt(alpha.conjugate().det())


@pytest.mark.parametrize("case", range(4))
def test_entire_noncommuting_time_dependent_Duhamel_response_against_finite_variation(
    case,
):
    O = np.array(phase.OMEGA, dtype=float)
    H0 = np.array(independent_vertex(case), dtype=float)
    H1 = np.array(independent_vertex(case + 3), dtype=float)
    A = np.array(independent_vertex(case + 1), dtype=float)
    B = np.array(independent_vertex(case + 2), dtype=float)
    prepare = np.array(independent_symplectic(case), dtype=float)
    V0 = prepare @ prepare.T / 2
    finish = 0.19

    def source(t):
        return (1 + t) * np.exp(-3 * t * t)

    def integrate(epsilon):
        def rhs(t, flat):
            S = flat.reshape(4, 4)
            return (O @ (H0 + t * H1 + epsilon * source(t) * B) @ S).ravel()

        solution = solve_ivp(
            rhs,
            (0, finish),
            np.eye(4).ravel(),
            method="DOP853",
            rtol=2e-12,
            atol=2e-13,
            dense_output=True,
        )
        assert solution.success
        return solution

    baseline = integrate(0)
    end = baseline.y[:, -1].reshape(4, 4)
    Apulled = end.T @ A @ end

    def integrated_commutator(t):
        right = baseline.sol(t).reshape(4, 4)
        Bpulled = right.T @ B @ right
        return (
            source(t)
            * np.trace((Apulled @ O @ Bpulled - Bpulled @ O @ Apulled) @ V0)
            / 2
        )

    retarded, error = quad(integrated_commutator, 0, finish, epsabs=2e-12, epsrel=2e-12)
    assert error < 1e-9

    def mean(epsilon):
        S = integrate(epsilon).y[:, -1].reshape(4, 4)
        return np.trace(A @ S @ V0 @ S.T) / 2

    step = 2e-5
    finite = (mean(step) - mean(-step)) / (2 * step)
    finer = (mean(step / 2) - mean(-step / 2)) / step
    np.testing.assert_allclose(finer, retarded, rtol=2e-7, atol=2e-8)
    np.testing.assert_allclose(finite, finer, rtol=2e-7, atol=2e-8)


@cache
def independent_coefficient_hessian():
    Q, sigma, P, momentum = s.symbols(
        "independent_Q independent_sigma independent_P independent_momentum", real=True
    )
    radial = functional.complete_coefficient_hessian()[0]
    a = phase.a
    J, theta, E, ell, A, T = charts.J, charts.th, charts.E, charts.l, charts.A, charts.T
    w = -ell * E
    q = radial**2 / a**2
    pv, ps = P / a**3, momentum / a**3
    numerator = theta * pv - w * ps + 3 * theta * ell * sigma - (2 * E * q + 3 * T) * Q
    h = a**3 * (
        ps**2 / 2
        - ell * pv * sigma / 2
        + (q / 2 - 3 * ell**2 / 4) * sigma**2
        - q * Q**2
        - s.Rational(9, 2) * A * Q**2
        + numerator**2 / (4 * J)
    )
    return s.hessian(h, [Q, sigma, P, momentum]).applyfunc(s.cancel)


PARAMETERS = (charts.J, charts.th, charts.E, charts.l, charts.A, charts.T, phase.a)


def test_separately_entered_complete_canonical_current_Hessian():
    whole = functional.complete_coefficient_hessian()[1]
    assert (whole - independent_coefficient_hessian()).applyfunc(s.cancel) == s.zeros(4)


@pytest.mark.parametrize("parameter", PARAMETERS)
def test_all_seven_complete_first_coefficient_vertices_from_independent_action(
    parameter,
):
    expected = independent_coefficient_hessian().diff(parameter)
    actual = functional.data()["whole_first_coefficient_vertices"][str(parameter)]
    assert (expected - actual).applyfunc(s.cancel) == s.zeros(4)


@pytest.mark.parametrize("i,j", [(i, j) for i in range(7) for j in range(i, 7)])
def test_all_twenty_eight_complete_second_vertices_from_independent_action(i, j):
    expected = independent_coefficient_hessian().diff(PARAMETERS[i], PARAMETERS[j])
    actual = functional.data()["whole_second_coefficient_vertices"][
        str(PARAMETERS[i]) + "__" + str(PARAMETERS[j])
    ]
    assert (expected - actual).applyfunc(s.cancel) == s.zeros(4)


@pytest.mark.parametrize("case", range(6))
def test_both_TT_full_scale_vertices_from_independent_canonical_action(case):
    a, P = phase.a, functional.complete_coefficient_hessian()[0]
    h1, h2, p1, p2 = s.symbols("TT1 TT2 TTmomentum1 TTmomentum2", real=True)
    h = (p1**2 + p2**2) / (2 * a**3) + a * P**2 * (h1**2 + h2**2) / 2
    actual = functional.data()["each_complete_TT_Hessian_first_and_second_scale_vertex"]
    for order in range(3):
        full = s.hessian(s.diff(h, a, order), [h1, p1, h2, p2])
        assert full == s.diag(actual[order], actual[order])
        values = {a: s.Rational(case + 4, 3), P: s.Rational(case + 2, 5)}
        assert full.subs(values) == s.diag(
            actual[order].subs(values), actual[order].subs(values)
        )


@pytest.mark.parametrize("case", range(5))
def test_noncommuting_formal_logdet_against_independent_scalar_determinant_series(case):
    epsilon = s.Symbol("independent_formal_epsilon")
    m = s.Integer(case + 1)
    D0 = s.Matrix([[5 + m, 1, -1], [1, 7 + m, 2], [-1, 2, 11 + m]])
    K = s.Matrix([[m, 2, 3], [2, -m, 1], [3, 1, m + 1]])
    F = s.Matrix([[1, m, 2], [m, -2, -1], [2, -1, 3]])
    Q = s.Matrix([[2, -1, m], [-1, 1, 2], [m, 2, -2]]) / 3
    matrix = D0 + epsilon * K + epsilon**2 * Q
    determinant = s.expand(matrix.det() / D0.det())
    scalar_series = s.series(s.log(determinant), epsilon, 0, 6).removeO().expand()
    candidate = functional.logdet_coefficients([D0, K, Q], 5)
    for order in range(1, 6):
        assert s.cancel(candidate[order - 1] - scalar_series.coeff(epsilon, order)) == 0
    assert D0.inv() * K * D0.inv() * F != D0.inv() * F * D0.inv() * K
    wrong = functional.logdet_coefficients([D0, K + F, Q], 5)
    assert wrong != candidate


@pytest.mark.parametrize("case", range(5))
def test_time_dependent_full_symplectic_chart_keeps_connection_and_boundary(case):
    t = s.Symbol("chart_time", real=True)
    r = s.Rational(case + 1, 7)
    scale = s.diag(1 + t * t, 1 + r * t * t)
    lower = s.eye(4)
    lower[2:, :2] = s.Matrix([[r * t, t * t / 5], [t * t / 5, -t / 11]])
    transform = lower * s.diag(scale, scale.inv().T)
    O = phase.OMEGA
    H = independent_vertex(case) + t * independent_vertex(case + 1)
    new_H = (
        transform.inv().T * H * transform.inv()
        - O * transform.diff(t) * transform.inv()
    )
    direct_generator = (transform.diff(t) + transform * O * H) * transform.inv()
    assert (new_H - new_H.T).applyfunc(s.cancel) == s.zeros(4)
    assert (direct_generator - O * new_H).applyfunc(s.cancel) == s.zeros(4)
    assert (O * transform.diff(t) * transform.inv()).subs(
        t, s.Rational(1, 9)
    ) != s.zeros(4)
    # The symmetric canonical one-form differs from p dq by d(p.q)/2.
    z = s.Matrix(s.symbols("new_q1 new_q2 new_p1 new_p2", real=True))
    zdot = s.Matrix(s.symbols("new_dq1 new_dq2 new_dp1 new_dp2", real=True))
    old = transform.inv() * z
    olddot = transform.inv() * zdot + transform.inv().diff(t) * z
    boundary = (old[:2, :].dot(old[2:, :]) - z[:2, :].dot(z[2:, :])) / 2
    boundary_dot = s.diff(boundary, t) + sum(
        s.diff(boundary, x) * dx for x, dx in zip(z, zdot)
    )
    old_action = old[2:, :].dot(olddot[:2, :]) - (old.T * H * old)[0] / 2
    new_action = z[2:, :].dot(zdot[:2, :]) - (z.T * new_H * z)[0] / 2
    assert s.cancel(old_action - new_action - boundary_dot) == 0


@pytest.mark.parametrize("order", range(6))
def test_finite_normal_Taylor_subtraction_removes_all_required_jets(order):
    normal = s.symbols("normal0:4")
    polynomial = sum(
        (j + 1) * (normal[0] + 2 * normal[1] - normal[2] + 3 * normal[3]) ** j
        for j in range(order + 3)
    )
    full = s.Poly(polynomial, *normal)
    removed = sum(
        coefficient * s.prod(x**power for x, power in zip(normal, powers))
        for powers, coefficient in full.terms()
        if sum(powers) <= order
    )
    remainder = s.Poly(polynomial - removed, *normal)
    assert all(
        sum(powers) > order
        for powers, coefficient in remainder.terms()
        if coefficient != 0
    )
    lam = s.Symbol("positive_scale", positive=True)
    rescaled = s.expand(
        remainder.as_expr().subs({x: lam * x for x in normal}, simultaneous=True)
    )
    assert all(rescaled.coeff(lam, j) == 0 for j in range(order + 1))
    assert rescaled.coeff(lam, order + 1) != 0


def test_written_continuum_bound_not_physical_counterterm_selection_or_P8_closure():
    assert functional.bilinear_scaling_bound(4, 4) == 20
    assert 20 - 4 == 16
    orders = functional.data()["complete_chart_symbol_orders_in_q"]
    assert orders["outer"]["B"] == [1, 0, "zero_polynomial", 0]
    assert orders["central"]["B"] == [1, 0, 0, -1]
    assert all(
        orders[chart][field] == [0, 0, 0, 0] for chart in orders for field in ("K", "G")
    )
    assert audit.frontier() == audit.previous.frontier()
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 108
    assert audit.matching()[:-1] == audit.previous.matching()
    assert (
        "NOT_COVARIANT_COUNTERFUNCTIONAL_PHYSICAL_LOOP_MEAN_BOUND"
        in audit.ITEM["status"]
    )
    assert audit.require_state(audit.previous.STATE) == audit.previous.STATE


@pytest.mark.parametrize("case", range(4))
def test_full_squeezed_annihilator_norm_and_Weyl_Schrodinger_phase(case):
    transport = independent_symplectic(case + 4)
    alpha, beta = functional.bogoliubov(transport)
    K = (alpha.conjugate().T.inv() * beta.T).applyfunc(s.cancel)
    assert (K - K.T).applyfunc(s.cancel) == s.zeros(2)
    assert (K - beta * alpha.conjugate().inv()).applyfunc(s.cancel) == s.zeros(2)
    norm_matrix = (s.eye(2) - K.conjugate().T * K).applyfunc(s.cancel)
    assert norm_matrix[0, 0].is_positive is True
    assert s.cancel(norm_matrix.det()).is_positive is True
    assert s.cancel(norm_matrix.det() * (alpha * alpha.conjugate().T).det() - 1) == 0
    H = independent_vertex(case + 1)
    qq, qp, pp = H[:2, :2], H[:2, 2:], H[2:, 2:]
    h = (qq + pp + s.I * (qp.T - qp)) / 2
    pairing = (qq - pp + s.I * (qp + qp.T)) / 2
    assert h == h.conjugate().T and pairing == pairing.T
    derivative = phase.OMEGA * H * transport
    alpha_derivative = (
        derivative[:2, :2]
        + derivative[2:, 2:]
        + s.I * (derivative[2:, :2] - derivative[:2, 2:])
    ) / 2
    assert (
        alpha_derivative + s.I * (h * alpha + pairing * beta.conjugate())
    ).applyfunc(s.cancel) == s.zeros(2)
    determinant_phase = (
        -s.trace(alpha.conjugate().inv() * alpha_derivative.conjugate()) / 2
    )
    direct_vacuum_Schrodinger = (
        -s.I * (s.trace(h) + s.trace(pairing.conjugate() * K)) / 2
    )
    assert s.cancel(determinant_phase - direct_vacuum_Schrodinger) == 0
