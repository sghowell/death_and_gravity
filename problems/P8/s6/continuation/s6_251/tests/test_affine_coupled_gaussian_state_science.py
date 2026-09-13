"""Independent complete canonical, covariance, mode and preparation diagnostics."""

from functools import cache

import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import audit, gaussian, phase, uv
from scipy.integrate import solve_ivp


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_every_entire_exact_identity(name, value):
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    ), name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_every_scoped_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input_or_physical_promotion(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def independent_symplectic(case):
    u = s.Rational(case + 1, 9)
    L = s.Matrix([[1 + u, u / 5], [-u / 7, (1 - u * u / 35) / (1 + u)]])
    C = s.Matrix([[u / 3, -u / 11], [-u / 11, (1 + u) / 13]])
    D = s.Matrix([[u / 17, u / 19], [u / 19, -u / 23]])
    low = s.eye(4)
    low[2:, :2] = C
    up = s.eye(4)
    up[:2, 2:] = D
    return low * s.diag(L, L.inv().T) * up


@pytest.mark.parametrize("case", range(12))
def test_independent_exact_noncommuting_Gaussian_minimum_and_uncertainty(case):
    T = independent_symplectic(case)
    n1 = s.Rational(case + 2, 5)
    n2 = n1 if case % 3 == 0 else s.Rational(case + 3, 7)
    D = s.diag(n1, n2, n1, n2)
    M = T.inv().T * D * T.inv()
    V = gaussian.ground_covariance(M)
    J = s.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]])
    assert T * J * T.T == J
    assert V == T * T.T / 2
    assert V * J * V == J / 4
    assert s.trace(M * V) / 2 == (n1 + n2) / 2
    complex_test = s.Matrix(
        [
            1 + s.I * (case + 1),
            s.Rational(2, 3) - s.I,
            -2 + s.I / 7,
            s.Rational(case + 2, 5) + 2 * s.I,
        ]
    )
    value = s.simplify(
        (complex_test.conjugate().T * (V + s.I * J / 2) * complex_test)[0]
    )
    assert value.is_real is True and value >= 0
    squeeze = s.diag(
        s.Rational(6, 5), s.Rational(7, 4), s.Rational(5, 6), s.Rational(4, 7)
    )
    excited = T * squeeze * squeeze.T * T.T / 2
    assert s.trace(M * excited) / 2 > (n1 + n2) / 2
    assert gaussian.ground_covariance((case + 3) * M) == V
    if case < 6:
        U = independent_symplectic(case + 13)
        assert gaussian.ground_covariance(U.inv().T * M * U.inv()) == U * V * U.T


def spectral_covariance(M):
    eigenvalues, eigenvectors = np.linalg.eigh(M)
    assert np.min(eigenvalues) > 0
    root = (eigenvectors * np.sqrt(eigenvalues)) @ eigenvectors.T
    invroot = (eigenvectors / np.sqrt(eigenvalues)) @ eigenvectors.T
    J = np.array(phase.OMEGA, dtype=float)
    B = root @ J @ root
    d, q = np.linalg.eigh(-(B @ B))
    assert np.min(d) > 0
    absolute = (q * np.sqrt(d)) @ q.T
    return invroot @ absolute @ invroot / 2


@pytest.mark.parametrize("case", range(10))
def test_independent_full_matrix_spectral_route(case):
    T = independent_symplectic(case + 2)
    D = s.diag(
        s.Rational(case + 2, 7),
        s.Rational(case + 4, 3),
        s.Rational(case + 2, 7),
        s.Rational(case + 4, 3),
    )
    M = T.inv().T * D * T.inv()
    candidate = np.array(gaussian.ground_covariance(M), dtype=float)
    independent = spectral_covariance(np.array(M, dtype=float))
    np.testing.assert_allclose(candidate, independent, rtol=3e-10, atol=3e-10)


@pytest.mark.parametrize("case", range(8))
def test_independent_tensor_Gaussian_reference(case):
    T = s.Matrix(
        [
            [1, s.Rational(case + 1, 5)],
            [s.Rational(case + 2, 9), 1 + s.Rational((case + 1) * (case + 2), 45)],
        ]
    )
    nu = s.Rational(case + 3, 11)
    M = T.inv().T * (nu * s.eye(2)) * T.inv()
    V = gaussian.single_mode_covariance(M)
    J = s.Matrix([[0, 1], [-1, 0]])
    assert V == T * T.T / 2 and V * J * V == J / 4
    assert s.trace(M * V) / 2 == nu / 2
    assert gaussian.single_mode_covariance(7 * M) == V


@pytest.mark.parametrize("case", range(8))
def test_real_covariance_positivity_does_not_replace_complex_uncertainty(case):
    epsilon = s.Rational(case + 1, 30)
    V = epsilon * s.eye(4)
    J = phase.OMEGA
    z = s.Matrix([1, 0, s.I, 0])
    assert V.is_positive_definite is True
    assert s.simplify((z.conjugate().T * (V + s.I * J / 2) * z)[0]) < 0


@cache
def independent_dynamics():
    """Separately entered full current-shape action; finite diagnostic stress, not an approximation to the actual tiny profiles."""
    t, k = s.symbols("independent_time independent_momentum", real=True)
    a = (1 + t * t) ** 2
    H = s.diff(a, t) / a
    delta = 1 / (2 * (1 + t * t) ** 3)
    ell = 1 / (10 * (1 + t * t) ** 6)
    theta = H - t / (1 + t * t) ** 4
    E = 1 - 3 * delta
    w = -ell * E
    F = (
        theta * s.diff(E, t)
        - E * s.diff(theta, t)
        + H * E * theta
        - theta**2
        - w * w / 2
    )
    rho = s.Rational(1, 10**8) * (1 + t + 2 * t * t)
    pressure = s.Rational(1, 10**8) * (2 - t + t * t)
    A = -pressure
    Bp = -(rho + pressure) / 2
    Jc = (
        F
        + 1 / (50 * (1 + t * t) ** 6)
        + (21 * delta**2 - 3 * delta) * A / 2
        + (1 - 6 * delta) * Bp
    )
    Tc = rho - 3 * delta * pressure
    q = k * k / a**2
    v, sigma, pv, ps = s.symbols("metric matter metric_p matter_p", real=True)
    variables = s.Matrix([v, sigma, pv, ps])
    lc = theta * pv - w * ps + 3 * theta * ell * sigma - (2 * E * q + 3 * Tc) * v
    h = (
        ps**2 / 2
        - ell * pv * sigma / 2
        + (q / 2 - 3 * ell**2 / 4) * sigma**2
        - q * v * v
        - s.Rational(9, 2) * A * v * v
        + lc * lc / (4 * Jc)
    )
    omega = s.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]])
    weighted = omega * s.hessian(h, variables) - 3 * H * s.diag(0, 0, 1, 1)
    R = s.diag(1, 1, a**3, a**3)
    canonical = R * weighted * R.inv() + 3 * H * s.diag(0, 0, 1, 1)
    K = s.Matrix([[(2 * Jc + w * w) / theta**2, w / theta], [w / theta, 1]])
    G = s.Matrix([[(2 * F + w * w) / theta**2, w / theta], [w / theta, 1]])
    return {
        "t": t,
        "k": k,
        "a": a,
        "H": H,
        "delta": delta,
        "ell": ell,
        "theta": theta,
        "E": E,
        "F": F,
        "Jc": Jc,
        "A": A,
        "Tc": Tc,
        "q": q,
        "weighted": weighted,
        "canonical": canonical,
        "K": K,
        "G": G,
    }


@cache
def numerical_functions():
    d = independent_dynamics()
    t, k = d["t"], d["k"]
    return {
        name: s.lambdify((t, k), d[name], "numpy", cse=True)
        for name in ("canonical", "weighted", "K", "G", "a", "q")
    }


@cache
def exact_full_trajectory(momentum):
    funcs = numerical_functions()

    def rhs(t, z):
        return (
            np.asarray(funcs["canonical"](t, momentum), dtype=float) @ z.reshape(4, 4)
        ).ravel()

    sol = solve_ivp(
        rhs,
        (-0.5, 0.5),
        np.eye(4).ravel(),
        method="DOP853",
        rtol=2e-12,
        atol=2e-13,
        dense_output=True,
    )
    assert sol.success
    return sol


def integrated_preparation(momentum, order):
    funcs = numerical_functions()
    sol = exact_full_trajectory(momentum)
    points, weights = np.polynomial.legendre.leggauss(order)
    center, width = -7 / 16, 1 / 32
    Q = np.eye(4)[:2]
    M = np.zeros((4, 4))
    for point, weight in zip(points, weights, strict=True):
        t = center + width * point
        f = np.exp(-1 / (1 - point * point))
        a = float(funcs["a"](t, momentum))
        q = float(funcs["q"](t, momentum))
        K = np.asarray(funcs["K"](t, momentum), dtype=float)
        G = np.asarray(funcs["G"](t, momentum), dtype=float)
        A = np.asarray(funcs["canonical"](t, momentum), dtype=float)
        V = Q @ A
        E = a**3 * (V.T @ K @ V + Q.T @ (q * G + K) @ Q)
        S = sol.sol(t).reshape(4, 4)
        M += width * weight * f * f * (S.T @ E @ S)
    return (M + M.T) / 2


@pytest.mark.parametrize("momentum", (0.0, 0.5, 1.5, 3.0))
def test_full_bounce_trajectory_and_integrated_matrix_state(momentum):
    sol = exact_full_trajectory(momentum)
    J = np.array(phase.OMEGA, dtype=float)
    for t in np.linspace(-0.5, 0.5, 11):
        S = sol.sol(t).reshape(4, 4)
        np.testing.assert_allclose(S @ J @ S.T, J, rtol=2e-8, atol=2e-8)
    M = integrated_preparation(momentum, 96)
    independently_refined = integrated_preparation(momentum, 160)
    np.testing.assert_allclose(M, independently_refined, rtol=5e-9, atol=5e-10)
    V = spectral_covariance(M)
    p = np.sqrt(np.linalg.det(M))
    total = np.sqrt(-np.trace((J @ M) @ (J @ M)) / 2 + 2 * p)
    direct = (p * np.linalg.inv(M) - J @ M @ J) / (2 * total)
    np.testing.assert_allclose(V, direct, rtol=2e-8, atol=2e-8)
    np.testing.assert_allclose(V @ J @ V, J / 4, rtol=3e-8, atol=3e-8)
    assert np.linalg.eigvalsh(V + 1j * J / 2).min() > -3e-9
    S = sol.sol(0.4).reshape(4, 4)
    T = sol.sol(-0.2).reshape(4, 4)
    W = S @ (V + 1j * J / 2) @ T.T
    np.testing.assert_allclose(
        W - np.conjugate(W), 1j * S @ J @ T.T, rtol=1e-8, atol=1e-8
    )
    assert np.linalg.norm(V - np.eye(4) / 2) > 1e-3


@pytest.mark.parametrize(
    "time",
    (-s.Rational(2, 5), -s.Rational(3, 16), 0, s.Rational(3, 16), s.Rational(2, 5)),
)
@pytest.mark.parametrize("momentum", (0, 2, 100))
def test_independent_entire_weighted_physical_phase_at_and_away_from_crossing(
    time, momentum
):
    from p8_vacuum_affine_scalar_tame_propagator import charts

    d = independent_dynamics()
    at = {d["t"]: time, d["k"]: momentum}
    mapping = {
        charts.th: d["theta"],
        charts.E: d["E"],
        charts.l: d["ell"],
        charts.J: d["Jc"],
        charts.A: d["A"],
        charts.T: d["Tc"],
        charts.H: d["H"],
        charts.q: d["q"],
        phase.a: d["a"],
    }
    values = {key: s.factor(value.subs(at)) for key, value in mapping.items()}
    values[phase.kappa] = s.Integer(10**800)
    actual = phase.canonical_generator().subs(values)
    independent = d["canonical"].subs(at)
    assert (actual - independent).applyfunc(s.cancel) == s.zeros(4)
    weighted = d["weighted"].subs(at)
    density = values[phase.a] ** 3
    C = 10**400 * s.diag(1, 1, density, density)
    connection = 3 * values[charts.H] * s.diag(0, 0, 1, 1)
    assert (C * weighted * C.inv() + connection - independent).applyfunc(
        s.cancel
    ) == s.zeros(4)
    assert C * phase.OMEGA * C.T == 10**800 * density * phase.OMEGA


@pytest.mark.parametrize("case", range(6))
def test_independent_variable_frequency_full_Riccati_residual(case):
    t = s.Symbol("variable_frequency_time", real=True)
    k = s.Symbol("large_momentum", positive=True)
    omega = (s.Rational(2 + case, 5), 2 + t * t + s.Rational(case, 7))
    B = {
        0: s.eye(2),
        -2: s.Matrix([[s.Rational(1, 17), t / 23], [t / 23, s.Rational(2, 19)]]),
    }
    A = {0: s.Matrix([[t / 7, (t + case + 1) / 11], [-(t * t + 1) / 13, -t / 17]])}
    C = {
        2: s.diag(*[w * w for w in omega]),
        0: s.Matrix([[1 + t * t, t / 19], [t / 19, 2 + case + t]]),
    }
    coefficients = uv.riccati_series(B, A, C, omega, t, 3)
    # Extract coefficients of the WHOLE independent matrix equation before
    # inserting rational t-jets; this avoids a needless multivariate GCD.
    symbolic = {}
    jets = {}
    for exponent, matrix in coefficients.items():
        block = s.Matrix(
            2,
            2,
            lambda i, j, exponent=exponent: s.Function(f"full_z_{exponent}_{i}_{j}")(t),
        )
        symbolic[exponent] = block
        for i in range(2):
            for j in range(2):
                jets[block[i, j]] = matrix[i, j]
                jets[s.diff(block[i, j], t)] = s.diff(matrix[i, j], t)
    Z = sum((value * k**power for power, value in symbolic.items()), s.zeros(2))
    BM = sum((value * k**power for power, value in B.items()), s.zeros(2))
    AM = A[0]
    CM = C[2] * k * k + C[0]
    full = Z.diff(t) + CM + AM.T * Z + Z * AM + Z * BM * Z
    for entry in full:
        polynomial = s.Poly(s.expand(entry * k**6), k)
        for power in (2, 1, 0, -1):
            assert (
                s.cancel(polynomial.coeff_monomial(k ** (power + 6)).xreplace(jets))
                == 0
            )
    omitted = coefficients[1].diff(t) + AM.T * coefficients[1] + coefficients[1] * AM
    assert omitted != s.zeros(2)
    assert coefficients[0] != s.zeros(2)


@pytest.mark.parametrize("time", (-s.Rational(3, 16), s.Rational(3, 16)))
def test_complete_clean_chart_transition_has_bounded_inverse_and_expected_limit(time):
    from p8_vacuum_affine_scalar_tame_propagator import charts

    d = independent_dynamics()
    fields = (charts.th, charts.E, charts.l, charts.J, charts.A, charts.T, charts.H)
    values = (d["theta"], d["E"], d["ell"], d["Jc"], d["A"], d["Tc"], d["H"])
    mapping = {
        key: s.factor(value.subs(d["t"], time))
        for key, value in zip(fields, values, strict=True)
    }
    mapping[phase.a] = s.factor(d["a"].subs(d["t"], time))
    th, E = mapping[charts.th], mapping[charts.E]
    target = np.diag([float(-E / th), 1, float(-th / E), 1])
    norms = []
    for k in (10000, 20000):
        now = dict(mapping)
        now[charts.q] = s.Integer(k * k) / mapping[phase.a] ** 2
        T = np.array(phase.clean_transition().subs(now), dtype=float)
        scale = np.diag([np.sqrt(k), np.sqrt(k), 1 / np.sqrt(k), 1 / np.sqrt(k)])
        clean = scale @ T @ np.linalg.inv(scale)
        np.testing.assert_allclose(
            clean @ np.array(phase.OMEGA, dtype=float) @ clean.T,
            np.array(phase.OMEGA, dtype=float),
            atol=1e-10,
        )
        norms.append(np.linalg.norm(clean - target))
        assert (
            np.linalg.norm(clean) < 100 and np.linalg.norm(np.linalg.inv(clean)) < 100
        )
    assert norms[1] < 0.51 * norms[0]


@pytest.mark.parametrize("derivative", range(5))
@pytest.mark.parametrize("endpoint", (-1, 1))
def test_entire_smooth_sampling_has_no_boundary_terms(derivative, endpoint):
    x = s.Symbol("bump_coordinate", real=True)
    amplitude = s.exp(-2 / (1 - x * x))
    differentiated = s.diff(amplitude, x, derivative)
    assert s.limit(differentiated, x, endpoint, dir="+" if endpoint < 0 else "-") == 0


@pytest.mark.parametrize("frequency", (5, 15, 40))
def test_full_nonstationary_sampling_integration_by_parts(frequency):
    x = s.Symbol("integration_coordinate", real=True)
    omega = s.Integer(frequency)
    theta = omega * (x + x**3 / 10)
    amplitude = (1 + x / 7) * s.exp(-2 / (1 - x * x))
    first = s.diff(amplitude / (s.I * s.diff(theta, x)), x)
    second = s.diff(first / (s.I * s.diff(theta, x)), x)
    functions = [
        s.lambdify(x, expr, "numpy", cse=True) for expr in (amplitude, first, second)
    ]
    pts, wts = np.polynomial.legendre.leggauss(320)
    oscillation = np.exp(-1j * frequency * (pts + pts**3 / 10))
    values = [np.dot(wts, fun(pts) * oscillation) for fun in functions]
    np.testing.assert_allclose(values, values[0], rtol=2e-8, atol=3e-12)


def test_explicit_state_boundary_and_preserved_frontier():
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.matching()) == 107 and len(audit.frontier()) == 9
    assert audit.require_state(audit.STATE) == audit.STATE
    assert audit.require_sampling(-s.Rational(15, 32), -s.Rational(13, 32))
    assert audit.require_physical_normalization(10**800) == 10**800
    assert all(audit.require_observable(label) == label for label in audit.OBSERVABLES)
