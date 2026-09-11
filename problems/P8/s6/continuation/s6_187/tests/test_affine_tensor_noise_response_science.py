"""Independent tensor normalization, TT, actual-wave and boundary fixtures."""

from functools import cache

import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_tensor_noise_response import (
    audit,
    detectors,
    energy,
    projector,
    tree,
)
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm

t = s.Symbol("time", real=True)
bump_expression = s.exp(-1 / (1 - (5 * t / 2) ** 2))
bump_derivatives = [
    s.lambdify(t, s.diff(bump_expression, t, j), "numpy") for j in range(4)
]


def bump(time, order=0):
    return float(bump_derivatives[order](time)) if abs(time) < 0.4 else 0.0


def geometry(time):
    a = (1 + time * time) ** 2
    H = 4 * time / (1 + time * time)
    Hp = 4 * (1 - time * time) / (1 + time * time) ** 2
    return a, H, Hp


def forcing(time, k, compatible):
    if not compatible:
        return bump(time), bump(time, 1)
    a, H, Hp = geometry(time)
    omega2 = k * k / a**2
    q = bump(time, 2) + 3 * H * bump(time, 1) + omega2 * bump(time)
    qp = (
        bump(time, 3)
        + 3 * Hp * bump(time, 1)
        + 3 * H * bump(time, 2)
        - 2 * H * omega2 * bump(time)
        + omega2 * bump(time, 1)
    )
    return q, qp


@cache
def advanced(k, compatible):
    def equation(time, y):
        a, H, _ = geometry(time)
        return [
            y[1],
            forcing(time, k, compatible)[0] - 3 * H * y[1] - k * k / a**2 * y[0],
        ]

    result = solve_ivp(
        equation,
        (0.5, -0.5),
        [0.0, 0.0],
        method="DOP853",
        rtol=2e-11,
        atol=2e-12,
        dense_output=True,
        max_step=min(0.01, 0.25 / max(k, 1)),
    )
    assert result.success
    return result.sol


@cache
def fundamental(k):
    def equation(time, y):
        a, H, _ = geometry(time)
        return [
            y[1],
            -3 * H * y[1] - k * k / a**2 * y[0],
            y[3],
            -3 * H * y[3] - k * k / a**2 * y[2],
        ]

    result = solve_ivp(
        equation,
        (-0.5, 0.5),
        [1.0, 0.0, 0.0, 1.0],
        method="DOP853",
        rtol=2e-11,
        atol=2e-12,
        dense_output=True,
        max_step=min(0.01, 0.25 / max(k, 1)),
    )
    assert result.success
    return result.sol


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_exact_residual(name):
    value = audit.residuals()[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    ), name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_unsupported_scope_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("name", list(audit.packets()))
def test_all_continuous_proof_packets(name):
    assert all(audit.packets()[name].get("gates", {}).values())


def test_scope_counts_and_gates():
    assert all(v is True for v in audit.gates().values())
    assert len(audit.residuals()) == 40 and audit.scalar_entry_count() == 79
    assert len(audit.gates()) == 32 and audit.rejected_inputs() == 109
    assert len(audit.controls()) == 9
    assert audit.validate_scope(audit.frontier(), audit.matching())


@pytest.mark.parametrize("fixture", range(6))
def test_independent_unimodular_tensor_kinetic(fixture):
    rng = np.random.default_rng(9100 + fixture)
    n = rng.normal(size=3)
    n /= np.linalg.norm(n)
    P = np.eye(3) - np.outer(n, n)
    A = rng.normal(size=(3, 3))
    A = (A + A.T) / 2
    gamma = P @ A @ P - P * np.trace(P @ A) / 2
    # One full matrix direction is sufficient for an independent exact exponential path.
    H = 0.7 - fixture / 5
    for epsilon in (0.001, 0.01, 0.1):
        spatial = expm(epsilon * gamma)
        inverse = expm(-epsilon * gamma)
        derivative = 2 * H * spatial + epsilon * gamma @ spatial
        K = inverse @ derivative / 2
        assert abs(np.linalg.det(spatial) - 1) < 1e-12
        assert abs(np.trace(K) - 3 * H) < 1e-12
        kinetic = -(np.trace(K) ** 2 - np.trace(K @ K)) / 2 + 3 * H**2
        assert np.isclose(kinetic, epsilon**2 * np.trace(gamma @ gamma) / 8, atol=1e-12)
    assert tree.data()["gates"]["actual_tensor_speed_one"]


@pytest.mark.parametrize("fixture", range(10))
def test_rotated_TT_projection_is_orthogonal_contraction(fixture):
    rng = np.random.default_rng(9300 + fixture)
    n = rng.normal(size=3)
    n /= np.linalg.norm(n)
    P = np.eye(3) - np.outer(n, n)
    A = rng.normal(size=(3, 3))
    A = (A + A.T) / 2
    B = rng.normal(size=(3, 3))
    B = (B + B.T) / 2
    TT = lambda M: P @ M @ P - P * np.trace(P @ M) / 2
    T = TT(A)
    assert np.linalg.norm(T - T.T) < 1e-12
    assert np.linalg.norm(n @ T) < 1e-12 and abs(np.trace(T)) < 1e-12
    assert np.linalg.norm(TT(T) - T) < 1e-12
    assert abs(np.sum(T * B) - np.sum(A * TT(B))) < 1e-12
    assert np.linalg.norm(T) <= np.linalg.norm(A) + 1e-12
    k = (fixture + 1) * n
    Q = np.dot(k, k) * np.eye(3) - np.outer(k, k)
    local = Q @ A @ Q - Q * np.trace(Q @ A) / 2
    assert np.allclose(local, np.dot(k, k) ** 2 * T)


def test_local_TT_differential_operator_independently():
    x, y, z = s.symbols("x y z", real=True)
    coordinates = (x, y, z)
    potential = s.Matrix(
        [
            [x * x * y * y * z**4, x * y**3 * z**2, 0],
            [x * y**3 * z**2, x**4 * y * z * z, 0],
            [0, 0, x * x * y**4 * z * z],
        ]
    )
    laplace = lambda f: sum(s.diff(f, variable, 2) for variable in coordinates)
    Q = lambda i, j, f: (
        s.diff(f, coordinates[i], coordinates[j]) - (laplace(f) if i == j else 0)
    )
    # This polynomial is an algebra fixture. The proof applies the same local operator to C-infinity compact potentials.
    local = s.Matrix(
        3,
        3,
        lambda i, j: sum(
            Q(i, k, Q(j, l, potential[k, l])) - Q(i, j, Q(k, l, potential[k, l])) / 2
            for k in range(3)
            for l in range(3)
        ),
    )
    assert any(value != 0 for value in local)
    assert s.expand(s.trace(local)) == 0
    assert all(s.expand(value) == 0 for value in local - local.T)
    for j in range(3):
        assert s.expand(sum(s.diff(local[i, j], coordinates[i]) for i in range(3))) == 0
    assert "compact" in projector.data()["local_compact_construction"]


@pytest.mark.parametrize("k", (0.0, 1e-6, 0.5, 4.0, 40.0, 400.0))
def test_actual_compact_compatible_wave_detector(k):
    v = advanced(k, True)
    grid = np.linspace(-0.5, 0.5, 301)
    expected = np.array([[bump(time), bump(time, 1)] for time in grid]).T
    assert np.allclose(v(grid), expected, atol=3e-8, rtol=2e-7)
    assert np.max(np.abs(v(-0.5))) < 3e-8
    assert np.max(np.abs(v(0.5))) == 0
    # The forcing is not the zero test even though both homogeneous moments vanish.
    assert abs(forcing(0.0, k, True)[0]) > 1e-8


@pytest.mark.parametrize("k", (0.0, 0.1, 2.0, 20.0, 200.0))
def test_actual_weighted_Wronskian_and_moment_initial_data(k):
    f = fundamental(k)
    v = advanced(k, False)
    grid = np.linspace(-0.5, 0.5, 101)
    values = f(grid)
    a0 = (1 + 0.5**2) ** 2
    for i, time in enumerate(grid):
        a, _, _ = geometry(time)
        wronskian = values[0, i] * values[3, i] - values[1, i] * values[2, i]
        assert np.isclose(a**3 * wronskian, a0**3, rtol=1e-8, atol=1e-8)
    moments = [
        quad(
            lambda time, j=j: geometry(time)[0] ** 3 * f(time)[2 * j] * bump(time),
            -0.4,
            0.4,
            epsabs=1e-10,
            epsrel=1e-10,
            limit=1000,
        )[0]
        for j in range(2)
    ]
    assert np.allclose(
        v(-0.5), [moments[1] / a0**3, -moments[0] / a0**3], atol=2e-9, rtol=2e-7
    )


@pytest.mark.parametrize("k", (0.0, 1e-6, 0.5, 4.0, 40.0, 400.0))
def test_full_actual_adjoint_energy_and_smearing_norm(k):
    v = advanced(k, False)
    qnorm2 = quad(lambda time: bump(time) ** 2, -0.4, 0.4, epsabs=1e-12)[0]
    qp_norm2 = quad(lambda time: bump(time, 1) ** 2, -0.4, 0.4, epsabs=1e-12)[0]
    D2 = (1 + k * k + k**4) * qnorm2 + qp_norm2
    grid = np.linspace(-0.5, 0.5, 3001)
    samples = v(grid)
    norm_integrand = []
    for i, time in enumerate(grid):
        a, H, Hp = geometry(time)
        w, wt = samples[:, i]
        q, qp = forcing(time, k, False)
        wtt = q - 3 * H * wt - k * k / a**2 * w
        wttt = (
            qp
            - 3 * H * wtt
            - 3 * Hp * wt
            - k * k / a**2 * wt
            + 2 * H * k * k / a**2 * w
        )
        E2 = a**3 * wt * wt + a * k * k * w * w + a**3 * w * w
        assert E2 < 54**2 * qnorm2
        norm_integrand.append(w * w + wt * wt + wtt * wtt + wttt * wttt + k * k * w * w)
    # Numeric modal fixtures supplement the continuous full-momentum proof, not a momentum cutoff.
    integral = np.trapezoid(norm_integrand, grid)
    assert integral < 5000**2 * D2
    assert energy.data()["actual_second_time_majorant"] < 400
    assert energy.data()["actual_third_time_majorant"] < 4000


def test_generic_detector_does_not_satisfy_initial_compatibility():
    v = advanced(0.0, False)
    assert v(-0.5)[0] > 0.001
    assert v(-0.5)[1] < -0.001
    # The zero-extension of this advanced function is not even continuous at t0.
    # This rejects applicability of the proven H3-time stress-test bound, not the entire Proca/CD row.
    assert abs(v(-0.5)[0]) > 0.001
    assert "restricted" in detectors.data()["initial_state_boundary"]


@pytest.mark.parametrize("k", (0.0, 0.1, 2.0, 20.0))
def test_compact_detector_annihilates_both_actual_homogeneous_solutions(k):
    f = fundamental(k)
    a0 = geometry(-0.5)[0]
    v = advanced(k, True)
    for j in range(2):
        moment = quad(
            lambda time, j=j: (
                geometry(time)[0] ** 3 * f(time)[2 * j] * forcing(time, k, True)[0]
            ),
            -0.4,
            0.4,
            epsabs=1e-9,
            epsrel=1e-9,
            limit=1000,
        )[0]
        assert abs(moment) < 1e-7
    assert max(abs(x) for x in v(-0.5)) * a0**3 < 1e-7


@pytest.mark.parametrize("scaling", (1, 2, 7, 100))
def test_same_stress_noise_and_canonical_response_normalization(scaling):
    packet = detectors.data()
    canonical = packet["leading_h_variance_upper_per_D_squared"] * scaling**2
    metric = packet["leading_gamma_variance_upper_per_D_squared"] * scaling**2
    assert canonical < (5 * s.Rational(scaling, 10**372)) ** 2
    assert metric < s.Rational(scaling**2, 10**1542)
    assert metric == 4 * canonical / s.Integer(10) ** 800
    assert "finite-coupling remainder" in packet["finite_order_boundary"]
