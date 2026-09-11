"""Independent actual-state, covariant full-source and normalization fixtures."""

from functools import cache

import mpmath as mp
import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_cd_source_noise import (
    audit,
    covariance,
    derivatives,
    force,
    modes,
)
from p8_vector_hadamard import preparation
from p8_vector_state import wkb
from scipy.integrate import solve_ivp


def mpq(value):
    r = s.Rational(value)
    return mp.mpf(int(r.p)) / int(r.q)


def switch_parts(x, n=1024):
    q = (1 - x) / x
    power = q**n
    complement = power / (1 + power)
    w = n * x * x * mp.exp(-n * x * x)
    wx = 2 * n * x * mp.exp(-n * x * x) * (1 - n * x * x)
    Tx = n * q ** (n - 1) / (x * x * (1 + power) ** 2)
    return 1 + complement * (w - 1), Tx * (1 - w) + complement * wx, complement, w


def covariant_jet(time, jet):
    a = (1 + time * time) ** 2
    H = 4 * time / (1 + time * time)
    p = [1 + jet[1], *jet[2:5]]
    h = [[mp.mpf(0) for _ in range(4)] for _ in range(4)]
    at = 5
    for i in range(4):
        for j in range(i, 4):
            h[i][j] = h[j][i] = jet[at]
            at += 1
    raw = [row[:] for row in h]
    for i in range(1, 4):
        h[0][i] = h[i][0] = h[0][i] - H * p[i]
        h[i][i] -= a * a * H * p[0]
    inv = [1, -1 / (a * a), -1 / (a * a), -1 / (a * a)]
    X = sum(inv[i] * p[i] ** 2 for i in range(4))
    box = sum(inv[i] * h[i][i] for i in range(4))
    Z = sum(inv[i] * inv[j] * p[i] * p[j] * h[i][j] for i in range(4) for j in range(4))
    return a, H, p, raw, X, box, Z


def literal_source(time, jet):
    _, _, p, _, X, box, Z = covariant_jet(time, jet)
    u = time + jet[0]
    B, BX, _, _ = switch_parts(X)
    Rminus = B * (X - 1) / (1 + u * u) ** 3
    R = 1 + Rminus
    RX = (BX * (X - 1) + B) / (1 + u * u) ** 3
    Ru = -6 * u * Rminus / (1 + u * u)
    H = 4 * u / (1 + u * u)
    Q = (
        (box - 3 * H * X) / X
        + 3 * Ru / (4 * R)
        + (-1 / X**2 + 3 * RX / (2 * R * X)) * Z
    )
    return [pi * Rminus * Q for pi in p]


def directional(time, jet, eta, component):
    return mp.diff(
        lambda e: literal_source(time, [x + e * y for x, y in zip(jet, eta)])[
            component
        ],
        mp.mpf(0),
    )


@cache
def actual_evolution(kind, nu):
    initial = preparation.initial_data(kind, 1000, nu)
    W = float(initial["all_order_frequency"].evalf(30))
    rate = float(initial["all_order_half_log_rate"].evalf(30))
    k2 = float(initial["comoving_momentum_squared"])
    f = 1 / np.sqrt(2 * W)
    initial_y = np.array([f, (-1j * W - rate) * f], dtype=complex)

    def coefficients(t):
        a = (1 + t * t) ** 2
        H = 4 * t / (1 + t * t)
        Hd = 4 * (1 - t * t) / (1 + t * t) ** 2
        omega2 = 1000**2 + k2 / a**2
        z = (k2 / a**2) / omega2
        d = H / 2 if kind == "transverse" else (0.5 + z) * H
        U = (
            Hd / 2 + H * H / 4
            if kind == "transverse"
            else (0.5 + z) * Hd + (0.25 - z + 3 * z * z) * H * H
        )
        return a, omega2, d, U

    def rhs(t, y):
        _, omega2, _, U = coefficients(t)
        return np.array([y[1], -(omega2 - U) * y[0]])

    grid = np.linspace(-0.5, 0.5, 151)
    result = solve_ivp(
        rhs,
        (-0.5, 0.5),
        initial_y,
        method="DOP853",
        rtol=2e-10,
        atol=1e-12,
        t_eval=grid,
    )
    assert result.success
    return initial, result, coefficients


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
    assert all(audit.packets()[name]["gates"].values())


def test_complete_scope_gates():
    assert all(v is True for v in audit.gates().values())
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert len(audit.frontier()) == 9


@pytest.mark.parametrize("kind", ("transverse", "longitudinal"))
@pytest.mark.parametrize("nu", (1000, 3000))
def test_actual_all_order_Cauchy_and_evolved_modes(kind, nu):
    initial, result, coefficients = actual_evolution(kind, nu)
    assert initial["active_higher_derivative_orders"] == ([] if nu == 1000 else [6])
    if nu == 3000:
        assert initial["all_order_frequency"] != initial["frozen_frequency"]
        assert initial["included_correction_terms"][6]["frequency_correction"] != 0
    for index, time in enumerate(result.t):
        a, omega2, d, _ = coefficients(time)
        f, fd = result.y[:, index]
        p = fd - d * f
        omega = np.sqrt(omega2)
        assert abs(f) ** 2 < 9 / omega
        assert abs(p) ** 2 < 27 * omega
        assert abs((f * np.conjugate(fd) - np.conjugate(f) * fd) - 1j) < 2e-6
        assert 1 <= a <= 25 / 16
    u, z = wkb.u, wkb.z
    U = wkb.frequency(kind)["U"]
    for time in (-0.5, 0, 0.5):
        a, omega2, _, literal = coefficients(time)
        q = omega2 - 1000**2
        actual = float(U.subs({u: time, z: q / omega2}))
        assert abs(actual - literal) < 1e-12


@pytest.mark.parametrize("nu", (1000, 3000))
@pytest.mark.parametrize("index", (0, 25, 75, 125, 150))
def test_full_actual_covariance_retains_cross_components(nu, index):
    initial, tmode, coeff = actual_evolution("transverse", nu)
    _, lmode, lcoeff = actual_evolution("longitudinal", nu)
    time = tmode.t[index]
    a, omega2, dT, _ = coeff(time)
    _, _, dL, _ = lcoeff(time)
    fT, fdT = tmode.y[:, index]
    fL, fdL = lmode.y[:, index]
    pL = fdL - dL * fL
    k = np.sqrt(float(initial["comoving_momentum_squared"]))
    gL = np.sqrt(a * 1000**2 / omega2)
    U = np.zeros((4, 3), dtype=complex)
    U[2, 0] = U[3, 1] = fT / (a * np.sqrt(a))
    U[1, 2] = fL / (a * gL)
    U[0, 2] = -1j * k * gL * pL / (a**3 * 1000**2)
    C = U @ U.conjugate().T
    eig = np.linalg.eigvalsh(C)
    assert eig[0] > -1e-14
    bound = (
        27 / np.sqrt(omega2) + 36 * k * k / (a * a * 1000**2 * np.sqrt(omega2))
    ) / a**3
    assert np.trace(C).real < bound
    if k:
        assert abs(C[0, 1]) > 0 and C[0, 0].real > 0
    else:
        assert C[0, 0] == 0
    f = np.array([1 + 1j, 2 - 1j, -1 + 2j, 3 + 0.25j])
    actual = np.vdot(f, C @ f).real
    assert -1e-12 < actual < bound * np.vdot(f, f).real
    assert a**6 * bound <= 4 * (27 / 1000 + 36 * k * k / 1000**3)
    assert abs(fdT - dT * fT) ** 2 < 27 * np.sqrt(omega2)


@pytest.mark.parametrize("fixture", range(8))
def test_complete_complex_switch_without_complement_cancellation(fixture):
    with mp.workdps(100):
        x = 1 + mp.mpf("0.1") * mp.exp(2j * mp.pi * (fixture + mp.mpf("0.125")) / 8)
        B, BX, C, w = switch_parts(x)
        assert 0 < abs(C) < 2 * (mp.mpf(1) / 9) ** 1024
        assert abs(B) < 3 and abs(w) < 1
        # Differentiate the tiny correction itself, not a rounded B-1.
        actual = mp.diff(lambda y: switch_parts(y)[2] * (switch_parts(y)[3] - 1), x)
        assert abs((actual - BX) / BX) < mp.mpf("1e-75")
        inner = 1 + (x - 1) / 2
        assert abs(switch_parts(inner)[1]) < 60


@pytest.mark.parametrize("time", (-s.Rational(1, 2), 0, s.Rational(1, 2)))
@pytest.mark.parametrize("fixture", range(4))
def test_full_covariant_source_on_complex_jet_ball(time, fixture):
    with mp.workdps(90):
        t = mpq(time)
        jet = [
            mp.mpf("0.02") * mp.exp(2j * mp.pi * (i + fixture + mp.mpf("0.25")) / 19)
            for i in range(15)
        ]
        a, H, p, raw, X, box, Z = covariant_jet(t, jet)
        expanded_Z = (
            p[0] ** 2 * raw[0][0]
            - 2 * p[0] * sum(p[i] * raw[0][i] for i in range(1, 4)) / a**2
        )
        expanded_Z += (
            sum(p[i] * p[j] * raw[i][j] for i in range(1, 4) for j in range(1, 4))
            / a**4
        )
        expanded_Z += H * p[0] * sum(p[i] ** 2 for i in range(1, 4)) / a**2
        assert abs(Z - expanded_Z) < mp.mpf("1e-80")
        assert abs(
            box
            - (raw[0][0] + 3 * H * p[0] - sum(raw[i][i] for i in range(1, 4)) / a**2)
        ) < mp.mpf("1e-80")
        assert abs(X - 1) < mp.mpf("0.05") and abs(Z) < mp.mpf("0.03")
        assert abs(box) < 5
        assert max(abs(v) for v in literal_source(t, jet)) < 10


@pytest.mark.parametrize("time", (-s.Rational(1, 2), 0, s.Rational(1, 2)))
def test_exact_reference_and_complete_nearby_directional_source(time):
    with mp.workdps(100):
        t = mpq(time)
        zero = [mp.mpf(0)] * 15
        eta = [mp.cos(i + 1) for i in range(15)]
        assert literal_source(t, zero) == [0] * 4
        assert max(abs(directional(t, zero, eta, c)) for c in range(4)) < mp.mpf(
            "1e-70"
        )
        for delta in (mp.mpf("0.01"), mp.mpf("1e-14")):
            jet = [delta * mp.sin(i + 1) for i in range(15)]
            value = [directional(t, jet, eta, c) for c in range(4)]
            testnorm = mp.sqrt(sum(x * x for x in eta))
            assert mp.sqrt(sum(abs(x) ** 2 for x in value)) < 2048 * delta * testnorm
            # One actual spatial-jet direction including derivative of eta.
            dx = [delta * mp.cos(i + 2) for i in range(15)]
            eta_dx = [mp.sin(i + 2) for i in range(15)]
            fullnorm = mp.sqrt(sum(x * x for x in eta + eta_dx))

            def moved(e, component, jet=jet, dx=dx, eta_dx=eta_dx):
                return directional(
                    t,
                    [x + e * y for x, y in zip(jet, dx)],
                    [x + e * y for x, y in zip(eta, eta_dx)],
                    component,
                )

            mixed = [mp.diff(lambda e, c=c: moved(e, c), mp.mpf(0)) for c in range(4)]
            assert (
                mp.sqrt(sum(abs(x) ** 2 for x in mixed))
                < mpq(derivatives.GRAD) * delta * fullnorm
            )


@pytest.mark.parametrize("delta", (s.Rational(1, 100), s.Rational(1, 10**14)))
@pytest.mark.parametrize("test_norm", (1, 2, s.Rational(3, 7)))
def test_both_force_normalizations_with_exact_large_hierarchy(delta, test_norm):
    m, K = modes.MASS, modes.KAPPA
    D = 108 * m * 2048**2 + 144 * derivatives.GRAD**2 / m
    raw_norm = test_norm / s.sqrt(K)
    normalized_variance = D * delta**2 * raw_norm**2 / K
    canonical_variance = s.cancel(K * K * normalized_variance)
    assert canonical_variance == D * delta**2 * test_norm**2
    assert canonical_variance < (2 * 10**8 * delta * test_norm) ** 2
    assert D * delta**2 / K < (2 * s.Rational(1, 10**392) * delta) ** 2
    assert canonical_variance != D * delta**2 * test_norm**2 / K
    if delta == force.SMALL_DELTA:
        assert canonical_variance < (2 * s.Rational(1, 10**6) * test_norm) ** 2
        assert 4 * 10**10 * delta**2 == 4 * s.Rational(1, 10**18)


def test_canonical_mean_dictionary_and_no_metric_noise_conclusion():
    K = modes.KAPPA
    raw_field, raw_test = s.symbols("raw_field raw_test", real=True)
    d = s.Symbol("delta", positive=True)
    full = K * 4 * 10**10 * d * d * raw_field * raw_test
    assert (
        full.subs({raw_field: 3 / s.sqrt(K), raw_test: 7 / s.sqrt(K)})
        == 84 * 10**10 * d * d
    )
    assert force.data()["small_ball_boundary"].startswith("These are source-sector")
    assert "not a state-independent QEI" in covariance.data()["not_a_QEI"]
    assert audit.controls()["zero_reference_source_noise_not_zero_metric_noise"]
