"""Independent interval, canonical, sampled-state and physical-observable diagnostics."""

import math
from functools import cache

import numpy as np
import pytest
import sympy as s
from p8_affine import verify as certificate
from p8_vacuum_affine_coupled_gaussian_state import gaussian, phase
from p8_vacuum_affine_quantitative_gaussian_window import (
    audit,
    physical,
    reference,
    source,
    state,
)
from p8_vacuum_affine_scalar_tame_propagator import charts
from scipy import integrate, linalg, special


@pytest.mark.parametrize("name", tuple(audit.packets()))
def test_all_complete_exact_packets(name):
    data = audit.packets()[name]
    assert certificate.certify_residuals(data["checks"])
    assert all(bool(value) for value in data["gates"].values())


@pytest.mark.parametrize(
    "name,call,args",
    audit.bad_cases(),
    ids=lambda value: value if isinstance(value, str) else None,
)
def test_every_unsupported_input_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_unchanged_frontier_and_explicit_quantum_probability_boundary():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 120
    assert audit.matching()[:-1] == audit.previous.matching()
    assert all(audit.gates().values())
    assert "OPEN" in audit.observable()["original_problem"]
    assert "not" in audit.observable()["probability_scope"]
    assert "S263" in audit.observable()["original_problem"]


@pytest.mark.parametrize("P", [source.LOW, source.HIGH, 10**120])
def test_complete_preparation_estimate_high_momentum_domain(P):
    assert audit.require_momentum(P) == P


@pytest.mark.parametrize("threshold", [source.EPSILON, 2 * source.EPSILON, 1])
def test_stated_marginal_tail_domain(threshold):
    assert audit.require_threshold(threshold) == threshold


@pytest.mark.parametrize("name,pivot,mixed,lo,hi", source.CHARTS)
def test_principal_LMIs_against_independent_generalized_eigenvalues(
    name, pivot, mixed, lo, hi
):
    data = reference.principal()["whole_actual_profile_dependent_principal_matrices"][
        name
    ]
    # Numerical eigenvalues are diagnostics only; all nonzero actual
    # profile terms are retained by the separate exact interval proof.
    functions = s.lambdify(
        (source.u, *source.PROFILE_VARIABLES),
        [data[key] for key in ("K", "G", "Kdot_plus_6HK", "Gdot_minus_2HG")],
        "numpy",
        cse=True,
    )
    for t in np.linspace(float(lo), float(hi), 37):
        K, G, VK, VG = functions(t, 0, 0, 0, 0)
        for Q, V in [(K, VK), (G, VG)]:
            eig = linalg.eigvalsh(Q, check_finite=True)
            assert eig[0] > 0.5 and eig[-1] < 32
            assert np.max(np.abs(linalg.eigvalsh(V, Q, check_finite=True))) < 20
    assert reference.principal()["whole_exact_smallest_positive_minor_enclosures"][
        name
    ] > s.Rational(1, 100)


@pytest.mark.parametrize("name,pivot,mixed,lo,hi", source.CHARTS[:2])
def test_complete_profile_chain_against_explicit_nonzero_polynomial_fixture(
    name, pivot, mixed, lo, hi
):
    # Deliberately finite nonzero profiles test algebra, not the actual
    # tiny profile hypothesis or a different certified preparation.
    K, _ = source.principal(pivot, mixed)
    u = source.u
    profiles = {source.rho: (1 + u * u) / 31, source.pressure: (1 + u**3) / 37}
    jets = {**profiles, source.rho1: 2 * u / 31, source.pressure1: 3 * u * u / 37}
    for entry in K:
        actual = source.dt(entry).subs(jets, simultaneous=True)
        independent = s.diff(entry.subs(profiles, simultaneous=True), u)
        assert s.factor(actual - independent) == 0


def independent_covariance(M):
    eigen, Q = linalg.eigh(M, check_finite=True)
    assert eigen[0] > 0
    half = (Q * np.sqrt(eigen)) @ Q.T
    inverse_half = (Q * (1 / np.sqrt(eigen))) @ Q.T
    omega = np.array(phase.OMEGA, float)
    frequencies, U = linalg.eigh(1j * half @ omega @ half, check_finite=True)
    modulus = (U * np.abs(frequencies)) @ U.conj().T
    covariance = inverse_half @ modulus @ inverse_half / 2
    assert np.max(np.abs(covariance.imag)) < 1e-8
    return covariance.real


def two_mode_formula(M):
    omega = np.array(phase.OMEGA, float)
    product = math.sqrt(linalg.det(M, check_finite=True))
    total = math.sqrt(-np.trace((omega @ M) @ (omega @ M)) / 2 + 2 * product)
    return (product * linalg.inv(M, check_finite=True) - omega @ M @ omega) / (
        2 * total
    )


@pytest.mark.parametrize("case", range(12))
def test_full_positive_covariance_matrix_equation_and_condition_bound_independently(
    case,
):
    _, _, exact = gaussian.fixture_preparation(case)
    M = np.array(exact, float)
    expected = np.array(gaussian.ground_covariance(exact), float)
    V = independent_covariance(M)
    np.testing.assert_allclose(V, expected, rtol=3e-9, atol=3e-9)
    omega = np.array(phase.OMEGA, float)
    np.testing.assert_allclose(4 * V @ M @ V, -omega @ M @ omega, rtol=2e-8, atol=2e-8)
    spectrum = linalg.eigvalsh(M, check_finite=True)
    upper = math.sqrt(spectrum[-1] / spectrum[0]) / 2
    assert linalg.eigvalsh(V, check_finite=True)[-1] <= upper * (1 + 1e-10)
    assert linalg.eigvalsh(V + 0.5j * omega, check_finite=True)[0] > -1e-9
    assert linalg.eigvalsh(0.4 * np.eye(4) + 0.5j * omega, check_finite=True)[0] < 0


@cache
def bare_gradient_function():
    return s.lambdify(source.u, source.F, "numpy", cse=True)


def fixture_coefficients(t, P, seed):
    a = (1 + t * t) ** 2
    H = 4 * t / (1 + t * t)
    delta = 1 / (2 * (1 + t * t) ** 3)
    E = 1 - 3 * delta
    theta = H - t / (1 + t * t) ** 4
    ell = 1 / (10 * (1 + t * t) ** 6)
    rho = seed * 1e-5 * (1 + t * t)
    pressure = seed * 1e-5 * (0.5 + t**3)
    F = float(bare_gradient_function()(t))
    J = (
        F
        + 1 / (50 * (1 + t * t) ** 6)
        - (21 * delta * delta - 3 * delta) * pressure / 2
        - (1 - 6 * delta) * (rho + pressure) / 2
    )
    return (
        a,
        P * P / (a * a),
        theta,
        E,
        ell,
        F,
        J,
        -pressure,
        rho - 3 * delta * pressure,
    )


@cache
def original_source_functions():
    variables = (
        phase.a,
        charts.q,
        charts.th,
        charts.E,
        charts.l,
        charts.J,
        charts.A,
        charts.T,
        charts.jets[0],
        charts.jets[1],
        charts.H,
    )
    matrices = [phase.canonical_generator(), phase.selection_matrix()]
    assert all(matrix.free_symbols <= set(variables) for matrix in matrices)
    return s.lambdify(variables, matrices, "numpy", cse=True)


def source_fixture_arguments(t, P, seed):
    a, q, theta, E, ell, _F, J, A, T = fixture_coefficients(t, P, seed)
    theta1 = 4 * (1 - t * t) / (1 + t * t) ** 2 - (1 - 7 * t * t) / (1 + t * t) ** 5
    E1 = 9 * t / (1 + t * t) ** 4
    return a, q, theta, E, ell, J, A, T, theta1, E1, 4 * t / (1 + t * t)


def independent_original_hamiltonian(x, t, P, seed):
    a, q, theta, E, ell, _F, J, A, T = fixture_coefficients(t, P, seed)
    v, sigma, Pv, Ps = x
    pv, ps = Pv / a**3, Ps / a**3
    lc = theta * pv + ell * E * ps + 3 * theta * ell * sigma - (2 * E * q + 3 * T) * v
    return a**3 * (
        ps * ps / 2
        - ell * pv * sigma / 2
        + (q / 2 - 3 * ell * ell / 4) * sigma * sigma
        - q * v * v
        - 4.5 * A * v * v
        + lc * lc / (4 * J)
    )


def independent_original_generator(t, P, seed):
    a, q, theta, E, ell, _F, J, A, T = fixture_coefficients(t, P, seed)
    Hess = np.zeros((4, 4))
    Hess[3, 3] = a**-3
    Hess[1, 2] = Hess[2, 1] = -ell / 2
    Hess[1, 1] = a**3 * (q - 1.5 * ell * ell)
    Hess[0, 0] = -(a**3) * (2 * q + 9 * A)
    lc = np.array([-(2 * E * q + 3 * T), 3 * theta * ell, theta / a**3, ell * E / a**3])
    Hess += a**3 * np.outer(lc, lc) / (2 * J)
    return np.array(phase.OMEGA, float) @ Hess


def independent_selection(t, P, seed):
    if not -15 / 32 < t < -13 / 32:
        raise ValueError(
            "Evaluate the selection norm only on its original outer support"
        )
    a, q, theta, E, ell, F, J, _A, _T = fixture_coefficients(t, P, seed)
    mixed = -ell * E / theta
    K = np.array([[2 * J / theta**2 + mixed * mixed, mixed], [mixed, 1]])
    G = np.array([[2 * F / theta**2 + mixed * mixed, mixed], [mixed, 1]])
    velocity = independent_original_generator(t, P, seed)[:2, :]
    return a**3 * (
        velocity.T @ K @ velocity + linalg.block_diag(q * G + K, np.zeros((2, 2)))
    )


@pytest.mark.parametrize("P", [3.0, 7.0, 11.0])
@pytest.mark.parametrize("seed", [0, 1])
def test_independent_full_hamiltonian_polarization_and_original_sampling_flow(P, seed):
    # Moderate-P finite algebra/ODE fixtures, not a computation of actual
    # giant-P covariances or a change of the fixed P8 state.
    basis = np.eye(4)
    t = -0.43
    singles = np.array([independent_original_hamiltonian(x, t, P, seed) for x in basis])
    Hess = np.empty((4, 4))
    for i in range(4):
        Hess[i, i] = 2 * singles[i]
        for j in range(i):
            Hess[i, j] = Hess[j, i] = (
                independent_original_hamiltonian(basis[i] + basis[j], t, P, seed)
                - singles[i]
                - singles[j]
            )
    omega = np.array(phase.OMEGA, float)
    np.testing.assert_allclose(
        omega @ Hess, independent_original_generator(t, P, seed), rtol=2e-12, atol=2e-12
    )
    original_generator, original_selection = original_source_functions()(
        *source_fixture_arguments(t, P, seed)
    )
    np.testing.assert_allclose(
        original_generator,
        independent_original_generator(t, P, seed),
        rtol=2e-12,
        atol=2e-12,
    )
    np.testing.assert_allclose(
        original_selection, independent_selection(t, P, seed), rtol=2e-12, atol=2e-12
    )

    def rhs(t, z):
        S = z[:16].reshape(4, 4)
        A = independent_original_generator(t, P, seed)
        x = 32 * (t + 7 / 16)
        if abs(x) < 1:
            f = math.exp(-1 / (1 - x * x))
            Md = f * f * S.T @ independent_selection(t, P, seed) @ S
        else:
            Md = np.zeros((4, 4))
        return np.concatenate(((A @ S).ravel(), Md.ravel()))

    initial = np.concatenate((np.eye(4).ravel(), np.zeros(16)))
    sol = integrate.solve_ivp(
        rhs,
        (-0.5, 0.0),
        initial,
        method="DOP853",
        max_step=1 / 256,
        rtol=2e-11,
        atol=2e-13,
    )
    assert sol.success
    S = sol.y[:16, -1].reshape(4, 4)
    M = sol.y[16:, -1].reshape(4, 4)
    M = (M + M.T) / 2
    assert linalg.eigvalsh(M, check_finite=True)[0] > 0
    np.testing.assert_allclose(S @ omega @ S.T, omega, rtol=1e-8, atol=1e-8)
    V = independent_covariance(M)
    np.testing.assert_allclose(V, two_mode_formula(M), rtol=3e-8, atol=3e-8)
    inverse = linalg.inv(S, check_finite=True)
    Mt = inverse.T @ M @ inverse
    center = np.array(
        [[0, 0, -1 / (2 * P * P), 0], [0, 1, 0, 0], [2 * P * P, 0, 0, 0], [0, 0, 0, 1]]
    )
    balance = (
        np.diag([math.sqrt(P), math.sqrt(P), 1 / math.sqrt(P), 1 / math.sqrt(P)])
        @ center
    )
    inverse_balance = linalg.inv(balance, check_finite=True)
    Mr = inverse_balance.T @ Mt @ inverse_balance
    Vr = balance @ S @ V @ S.T @ balance.T
    np.testing.assert_allclose(
        independent_covariance((Mr + Mr.T) / 2), Vr, rtol=3e-8, atol=3e-8
    )
    np.testing.assert_allclose(independent_covariance(7 * M), V, rtol=2e-9, atol=2e-9)


@pytest.mark.parametrize("P", [1e64, 1.5e64, 2e64])
def test_complete_physical_rows_against_direct_constraint_with_both_momenta(P):
    # The deliberately enlarged Tcorr fixture tests the complete algebra,
    # not membership in the actual tiny-profile interval box.
    data = physical.rows()
    symbols = {
        str(x): x
        for key in (
            "whole_hat_log_scale_row_times_sqrt_kappa",
            "whole_physical_lapse_row_times_sqrt_kappa",
        )
        for x in data[key].free_symbols
    }
    a = 1.2
    J = 1.52
    E = -0.49
    theta = 3e-60
    ell = 0.099
    T = 1e-200
    values = {
        "whole_comoving_P": P,
        "scale": a,
        "Jc": J,
        "E": E,
        "Theta": theta,
        "ell": ell,
        "Tcorr": T,
        "whole_reference_Sbar00": 0.7e-49,
        "whole_reference_Sbar01": -0.3e-49,
        "whole_reference_Sbar11": 0.2e-49,
    }
    bind = {symbols[name]: value for name, value in values.items()}
    x = np.array([0.8, -0.4, 0.3, 0.9])
    yb, ym = x[:2] / P
    Pb = x[2] + P * P * (
        values["whole_reference_Sbar00"] * yb + values["whole_reference_Sbar01"] * ym
    )
    Ps = x[3] + P * P * (
        values["whole_reference_Sbar01"] * yb + values["whole_reference_Sbar11"] * ym
    )
    v = Pb / (2 * a * P * P)
    pv = -2 * P * P * yb / (a * a)
    n = (
        theta * pv
        + 3 * theta * ell * ym
        - (2 * E * P * P / (a * a) + 3 * T) * v
        + ell * E * Ps / a**3
    ) / (2 * J)
    vr = np.array(
        data["whole_hat_log_scale_row_times_sqrt_kappa"].subs(bind), float
    ).ravel()
    nr = np.array(
        data["whole_physical_lapse_row_times_sqrt_kappa"].subs(bind), float
    ).ravel()
    assert vr @ x == pytest.approx(v, rel=2e-13, abs=1e-140)
    assert nr @ x == pytest.approx(n, rel=2e-13, abs=1e-12)
    assert abs(n) > 1e10


@pytest.mark.parametrize("threshold", [0.2, 0.8, 1.7, 2.5])
def test_single_Gaussian_spectral_tail_by_independent_integral(threshold):
    sigma = 0.4
    density = lambda x: (
        math.exp(-x * x / (2 * sigma * sigma)) / (sigma * math.sqrt(2 * math.pi))
    )
    exact = (
        2 * integrate.quad(density, threshold, np.inf, epsabs=1e-13, epsrel=1e-12)[0]
    )
    wanted = special.erfc(threshold / (math.sqrt(2) * sigma))
    assert exact == pytest.approx(wanted, rel=2e-9, abs=1e-13)
    assert exact <= 2 * math.exp(-threshold * threshold / (2 * sigma * sigma))
    assert exact > 0


def test_small_single_point_tail_does_not_make_uniform_support():
    # Independent Gaussian-sequence diagnostic, not an asserted iid P8 field.
    single_tail = special.erfc(2 / math.sqrt(2))
    assert 0 < single_tail < 0.05
    assert (1 - single_tail) ** 1000 < 1e-19
    assert "not" in physical.moments()["whole_tail_boundary"].lower()
    exponent = physical.moments()["whole_marginal_Gaussian_tail_exponent_lower"]
    assert exponent > 10**29
    # Never numerically underflow exp(-1e29) and call its probability zero.


@pytest.mark.parametrize("scale", [0.2, 0.6, 1.0])
def test_complete_Weyl_volume_contact_with_nonzero_cross_covariance(scale):
    # Gaussian cubature of the FINITE polynomial only, not of a nonlinear
    # volume function over an unbounded Gaussian.
    covariance = np.array([[0.06, 0.017], [0.017, 0.04]]) * scale
    L = linalg.cholesky(covariance, lower=True, check_finite=True)
    nodes, weights = np.polynomial.hermite.hermgauss(7)
    value = 0.0
    for i, x in enumerate(nodes):
        for j, y in enumerate(nodes):
            v, n = math.sqrt(2) * L @ np.array([x, y])
            value += (
                weights[i]
                * weights[j]
                * (4.5 * v * v + 4.5 * v * n + 0.375 * n * n)
                / math.pi
            )
    expected = (
        4.5 * covariance[0, 0] + 4.5 * covariance[0, 1] + 0.375 * covariance[1, 1]
    )
    assert value == pytest.approx(expected, rel=2e-13)
    deleted = 4.5 * covariance[0, 0] + 0.375 * covariance[1, 1]
    assert abs(value - deleted) > 0.01 * scale


@pytest.mark.parametrize("case", [0, 3, 7])
def test_full_covariance_bound_also_at_degenerate_preparation_frequencies(case):
    S = gaussian.symplectic_fixture(case)
    M = S.inv().T * S.inv()
    actual = independent_covariance(np.array(M, float))
    expected = np.array(S * S.T / 2, float)
    np.testing.assert_allclose(actual, expected, rtol=3e-8, atol=3e-8)
    np.testing.assert_allclose(
        np.array(gaussian.ground_covariance(M), float), expected, rtol=3e-12, atol=3e-12
    )


@pytest.mark.parametrize("case", range(6))
def test_complete_noncommuting_energy_with_both_clean_momenta(case):
    K = np.array([[1.2 + case / 10, 0.3], [0.3, 2.0]])
    G = np.array([[1.7, 0.25], [0.25, 0.8 + case / 30]])
    assert np.linalg.norm(K @ G - G @ K) > 0.01
    P = 10.0 + case
    a = 1.1 + case / 20
    alpha = 0.4
    A = np.array([[0, alpha], [-alpha, 0]])
    x = np.array([0.7, -0.3, 0.6, 0.9])
    velocity = linalg.solve(
        K, x[2:] / a**3 - A @ x[:2] / P, assume_a="pos", check_finite=True
    )
    expected = (velocity @ K @ velocity + x[:2] @ G @ x[:2] / a**2) / 2
    expression = state.energy()["whole_scalar_fixed_reference_weighted_energy_metric"]
    values = {
        "positive_P": P,
        "positive_reference_scale": a,
        "full_K00": K[0, 0],
        "full_K01": K[0, 1],
        "full_K11": K[1, 1],
        "full_G00": G[0, 0],
        "full_G01": G[0, 1],
        "full_G11": G[1, 1],
        "full_skew_boundary": alpha,
    }
    Q = np.array(
        expression.subs(
            {symbol: values[str(symbol)] for symbol in expression.free_symbols}
        ),
        float,
    )
    assert x @ Q @ x / 2 == pytest.approx(expected, rel=2e-13)
    assert linalg.eigvalsh(Q, check_finite=True)[0] > 1e-4
    r = x / math.sqrt(P)
    assert P * r @ Q @ r / 2 == pytest.approx(expected, rel=2e-13)


@pytest.mark.parametrize("P", [64.0, 128.0, 512.0])
def test_both_complete_boundary_maps_from_original_generator_and_legendre_momenta(P):
    t = -3 / 16
    seed = 1
    arguments = source_fixture_arguments(t, P, seed)
    a, q = arguments[:2]
    generator, _ = original_source_functions()(*arguments)
    H = 4 * t / (1 + t * t)
    swap = np.array(
        [
            [0, 0, -1 / (2 * a**3 * q), 0],
            [0, 1, 0, 0],
            [2 * a**3 * q, 0, 0, 0],
            [0, 0, 0, 1],
        ]
    )
    swapdot = np.zeros((4, 4))
    swapdot[0, 2] = -H * swap[0, 2]
    swapdot[2, 0] = H * swap[2, 0]
    center = (swapdot + swap @ generator) @ linalg.inv(swap, check_finite=True)
    shears = []
    for flow in (generator, center):
        K = linalg.inv(a**3 * flow[:2, 2:], check_finite=True)
        B = -K @ flow[:2, :2]
        shear = np.eye(4)
        shear[2:, :2] = -(a**3) * (B + B.T) / 2
        shears.append(shear)
    direct = shears[1] @ swap @ linalg.inv(shears[0], check_finite=True)
    expression = phase.clean_transition()
    bindings = dict(
        zip(
            (
                phase.a,
                charts.q,
                charts.th,
                charts.E,
                charts.l,
                charts.J,
                charts.A,
                charts.T,
                charts.jets[0],
                charts.jets[1],
                charts.H,
            ),
            arguments,
        )
    )
    bindings[charts.H] = H
    expected = np.array(expression.subs(bindings), float)
    np.testing.assert_allclose(direct, expected, rtol=3e-10, atol=3e-9)
    balance = np.diag([math.sqrt(P), math.sqrt(P), 1 / math.sqrt(P), 1 / math.sqrt(P)])
    balanced = balance @ direct @ linalg.inv(balance, check_finite=True)
    assert linalg.svdvals(balanced, check_finite=True)[0] < 16
    assert (
        linalg.svdvals(linalg.inv(balanced, check_finite=True), check_finite=True)[0]
        < 16
    )
