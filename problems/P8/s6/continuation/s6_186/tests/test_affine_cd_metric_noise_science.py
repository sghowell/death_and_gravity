"""Independent full-stress, actual-preparation, pair and continuum fixtures."""

from functools import cache

import mpmath as mp
import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_cd_metric_noise import (
    audit,
    noise,
    reference,
    stress,
)
from p8_vector_hadamard import preparation, series
from p8_vector_state import wkb


def mpq(value):
    r = s.Rational(value)
    return mp.mpf(int(r.p)) / int(r.q)


@cache
def reference_functions(kind):
    return [
        s.lambdify((wkb.u, wkb.z), series.coefficient(kind, n), "mpmath")
        for n in range(1, 5)
    ]


def W8(time, k, kind):
    a = (1 + time * time) ** 2
    omega = mp.sqrt(1000**2 + k * k / a**2)
    z = (k * k / a**2) / omega**2
    ratio = 1 + sum(
        P(time, z) / omega ** (2 * n)
        for n, P in enumerate(reference_functions(kind), 1)
    )
    return omega * ratio


def ref_readouts(time, k, kind):
    a = (1 + time * time) ** 2
    H = 4 * time / (1 + time * time)
    omega = mp.sqrt(1000**2 + k * k / a**2)
    z = (k * k / a**2) / omega**2
    W = W8(time, k, kind)
    Wp = mp.diff(lambda u: W8(u, k, kind), time)
    d = H / 2 if kind == "transverse" else (mp.mpf("0.5") + z) * H
    f = 1 / mp.sqrt(2 * W)
    p = (-1j * W - Wp / (2 * W) - d) * f
    if kind == "transverse":
        Z = [p, 0, 0, 0, 1j * k / a * f, 0, 0, 1000 * f, 0, 0]
    else:
        Z = [
            0,
            0,
            1000 / omega * p,
            0,
            0,
            0,
            -1j * k / (a * omega) * p,
            0,
            0,
            omega * f,
        ]
    return a, omega, z, W, f, p, Z


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


def test_scope_gates_and_actual_counts():
    assert all(v is True for v in audit.gates().values())
    assert len(audit.residuals()) == 66 and audit.scalar_entry_count() == 81
    assert len(audit.gates()) == 44 and audit.rejected_inputs() == 108
    assert audit.validate_scope(audit.frontier(), audit.matching())


@pytest.mark.parametrize("a", range(4))
@pytest.mark.parametrize("b", range(4))
def test_full_stress_matrix_and_positive_energy(a, b):
    M = stress.data()["quadratic_matrices"][a, b]
    numeric = np.array(M, dtype=float)
    assert np.max(np.abs(np.linalg.eigvalsh(numeric))) <= 0.5 + 1e-14
    fields = s.Matrix(
        s.symbols("E0:3", real=True)
        + s.symbols("B0:3", real=True)
        + s.symbols("mA0:4", real=True)
    )
    assert (
        s.expand((fields.T * M * fields)[0] - stress.data()["physical_stress"][a, b])
        == 0
    )
    if a == b == 0:
        assert M == s.eye(10) / 2


@pytest.mark.parametrize("fixture", range(4))
def test_independent_covariant_action_metric_variation(fixture):
    eta = s.diag(1, -1, -1, -1)
    F = s.zeros(4)
    for a in range(4):
        for b in range(a + 1, 4):
            F[a, b] = s.Rational((a + 1) * (b + 2) + fixture, 7)
            F[b, a] = -F[a, b]
    V = s.Matrix([s.Rational(i + 1 + fixture, 5) for i in range(4)])
    H = s.Matrix(
        4, 4, lambda a, b: s.Rational(1 + min(a, b) + 2 * max(a, b) + fixture, 11)
    )
    e = s.Symbol("epsilon", real=True)
    inv = eta + e * H
    Fsquare = sum(
        inv[a, c] * inv[b, d] * F[a, b] * F[c, d]
        for a in range(4)
        for b in range(4)
        for c in range(4)
        for d in range(4)
    )
    lag = -Fsquare / 4 + (V.T * inv * V)[0] / 2
    density = (1 - e * s.trace(eta * H) / 2) * lag
    E = s.Matrix([F[0, i + 1] for i in range(3)])
    B = s.Matrix([F[2, 3], F[3, 1], F[1, 2]])
    values = [*E, *B, *V]
    names = (
        s.symbols("E0:3", real=True)
        + s.symbols("B0:3", real=True)
        + s.symbols("mA0:4", real=True)
    )
    actual = stress.data()["physical_stress"].subs(dict(zip(names, values)))
    expected = sum(H[a, b] * actual[a, b] for a in range(4) for b in range(4)) / 2
    assert s.expand(s.diff(density, e).subs(e, 0) - expected) == 0


@pytest.mark.parametrize("kind", ("transverse", "longitudinal"))
@pytest.mark.parametrize("nu", (1000, 3000, 10**7))
@pytest.mark.parametrize("center", (-s.Rational(1, 2), 0, s.Rational(1, 2)))
def test_full_complex_time_reference_and_constraint_amplitudes(kind, nu, center):
    with mp.workdps(80):
        k = mp.mpf(25) / 16 * mp.sqrt(mp.mpf(nu) ** 2 - 1000**2)
        t0 = mpq(center)
        outer = mpq(reference.OUTER)
        for j in range(6):
            t = t0 + outer * mp.exp(2j * mp.pi * (j + mp.mpf("0.3")) / 6)
            a = (1 + t * t) ** 2
            a0 = (1 + t0 * t0) ** 2
            omega = mp.sqrt(1000**2 + k * k / a**2)
            om0 = mp.sqrt(1000**2 + k * k / a0**2)
            z = (k * k / a**2) / omega**2
            W = W8(t, k, kind)
            assert abs(a0 * a0 / (a * a) - 1) < mp.mpf("0.001")
            assert abs(omega**2 / om0**2 - 1) < mp.mpf("0.001")
            assert abs(omega) > mp.mpf("0.99") * nu
            assert omega.real > mp.mpf("0.99") * nu
            assert abs(z) < mp.mpf("1.01")
            assert abs(W / omega - 1) < mp.mpf("0.001")
            assert W.real > nu / 2 and abs(W) < 3 * nu
            inner = t0 + (t - t0) / 2
            _, _, _, _, _, _, Z = ref_readouts(inner, k, kind)
            assert max(abs(x) for x in Z) < 1000 * mp.sqrt(nu)
            assert mp.sqrt(sum(abs(x) ** 2 for x in Z)) < 4000 * mp.sqrt(nu)
            l = mp.mpf("1.3") * k + mp.mpf("0.5")
            mu = mp.sqrt(1000**2 + l * l / (mp.mpf(25) / 16) ** 2)
            G = 1 / (W8(inner, k, kind) + W8(inner, l, kind))
            assert abs(G) < 2 / (nu + mu)


@pytest.mark.parametrize("kind", ("transverse", "longitudinal"))
@pytest.mark.parametrize("nu", (1000, 3000))
def test_actual_all_order_initial_mixing_not_reference_reset(kind, nu):
    with mp.workdps(100):
        packet = preparation.initial_data(kind, 1000, nu)
        k = mp.sqrt(mpq(packet["comoving_momentum_squared"]))
        t = mp.mpf("-0.5")
        W = mp.mpf(str(packet["all_order_frequency"].evalf(100)))
        rate = mp.mpf(str(packet["all_order_half_log_rate"].evalf(100)))
        f = 1 / mp.sqrt(2 * W)
        fd = (-1j * W - rate) * f
        refW = W8(t, k, kind)
        refWp = mp.diff(lambda u: W8(u, k, kind), t)
        rf = 1 / mp.sqrt(2 * refW)
        rfd = (-1j * refW - refWp / (2 * refW)) * rf
        basis = mp.matrix([[rf, mp.conj(rf)], [rfd, mp.conj(rfd)]])
        alpha, beta = mp.lu_solve(basis, mp.matrix([f, fd]))
        assert abs(alpha * rf + beta * mp.conj(rf) - f) < mp.mpf("1e-90")
        assert abs(alpha * rfd + beta * mp.conj(rfd) - fd) < mp.mpf("1e-85")
        assert abs(abs(alpha) ** 2 - abs(beta) ** 2 - 1) < mp.mpf("1e-85")
        assert (
            abs(alpha) < 2
            and abs(beta)
            < mpq(reference.data()["source_pinned_constants"]["B6"]) / nu**6
        )
        assert packet["active_higher_derivative_orders"] == ([] if nu == 1000 else [6])
        assert abs(beta) > 0


@pytest.mark.parametrize("fixture", range(6))
def test_arbitrary_momentum_polarization_rotation_and_stress_pair(fixture):
    rng = np.random.default_rng(7300 + fixture)
    n = rng.normal(size=3)
    n /= np.linalg.norm(n)
    v = rng.normal(size=3)
    e1 = v - n * np.dot(n, v)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    a = 1 + (fixture + 1) / 12
    k = 300 + 400 * fixture
    omega = np.sqrt(1000**2 + k * k / a**2)
    f = (1 + 0.2j) / np.sqrt(2 * omega)
    p = (-0.3 - 1j * omega) * f
    cols = []
    for ep in (e1, e2):
        cols.append(
            np.concatenate(
                [ep * p, 1j * k / a * np.cross(n, ep) * f, [0], 1000 * ep * f]
            )
        )
    cols.append(
        np.concatenate(
            [
                n * 1000 / omega * p,
                np.zeros(3),
                [-1j * k / (a * omega) * p],
                n * omega * f,
            ]
        )
    )
    for Z in cols:
        assert (
            abs(np.vdot(Z, Z).real / 2 - (abs(p) ** 2 + omega**2 * abs(f) ** 2) / 2)
            < 1e-8
        )
    F = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    F = (F + F.T) / 2
    M = sum(
        F[a, b] * np.array(stress.data()["quadratic_matrices"][a, b], dtype=complex)
        for a in range(4)
        for b in range(4)
    )
    assert np.linalg.norm(M, 2) <= 2 * np.linalg.norm(F) + 1e-12
    for left in cols:
        for right in cols:
            assert (
                abs(left @ M @ right)
                <= 2 * np.linalg.norm(F) * np.linalg.norm(left) * np.linalg.norm(right)
                + 1e-8
            )


@pytest.mark.parametrize("fixture", range(6))
def test_independent_Fock_pair_normalization_conjugation_and_centering(fixture):
    rng = np.random.default_rng(9800 + fixture)
    levels = 5
    annihilator = np.diag(np.sqrt(np.arange(1, levels)), 1)
    eye = np.eye(levels)
    ops = [np.kron(annihilator, eye), np.kron(eye, annihilator)]
    vacuum = np.zeros(levels**2, dtype=complex)
    vacuum[0] = 1
    Q = np.zeros((levels**2, levels**2), dtype=complex)
    pair = np.zeros((2, 2), dtype=complex)
    for step in range(3):
        U = rng.normal(size=(10, 2)) + 1j * rng.normal(size=(10, 2))
        F = rng.normal(size=(4, 4))
        F = (F + F.T) / 2
        M = sum(
            F[a, b] * np.array(stress.data()["quadratic_matrices"][a, b], dtype=float)
            for a in range(4)
            for b in range(4)
        )
        fields = [
            sum(
                U[i, r] * ops[r] + np.conj(U[i, r]) * ops[r].conjugate().T
                for r in range(2)
            )
            for i in range(10)
        ]
        weight = (step + 1) / 7
        Q += weight * sum(
            M[i, j] * fields[i] @ fields[j] for i in range(10) for j in range(10)
        )
        pair += weight * U.conjugate().T @ M @ U.conjugate()
    assert np.linalg.norm(Q - Q.conjugate().T) < 1e-10
    mean = np.vdot(vacuum, Q @ vacuum)
    centered = Q @ vacuum - mean * vacuum
    actual = np.vdot(centered, centered).real
    expected = 2 * np.sum(np.abs(pair) ** 2)
    assert abs(actual - expected) < 1e-8 * max(1, expected)
    assert abs(centered[2 * levels] - np.sqrt(2) * pair[0, 0]) < 1e-10
    assert abs(centered[2] - np.sqrt(2) * pair[1, 1]) < 1e-10
    assert abs(centered[levels + 1] - 2 * pair[0, 1]) < 1e-10
    assert expected > 0 and abs(actual - expected / 2) > 1e-6
    shifted = Q + (fixture + 3) * np.eye(levels**2)
    shifted_mean = np.vdot(vacuum, shifted @ vacuum)
    assert np.linalg.norm(shifted @ vacuum - shifted_mean * vacuum - centered) < 1e-10


@pytest.mark.parametrize("fixture", range(4))
def test_independent_third_pair_integration_by_parts(fixture):
    # A positive generic phase checks the operator identity, not a different CD parent.
    with mp.workdps(55):
        c = mp.mpf(3 + fixture)
        phase = lambda t: c * t + t * t / 2
        inv = lambda t: 1 / (c + t)
        amp = lambda t: (1 + 1j * t) / (1 + t * t)
        test = lambda t: (1 - 4 * t * t) ** 4

        def b(t):
            return amp(t) * test(t)

        def first(t):
            return mp.diff(lambda u: inv(u) * b(u), t)

        def second(t):
            return mp.diff(lambda u: inv(u) * first(u), t)

        def third(t):
            return mp.diff(lambda u: inv(u) * second(u), t)

        direct = mp.quad(
            lambda t: b(t) * mp.exp(1j * phase(t)), [-mp.mpf("0.5"), 0, mp.mpf("0.5")]
        )
        integrated = 1j**3 * mp.quad(
            lambda t: third(t) * mp.exp(1j * phase(t)),
            [-mp.mpf("0.5"), 0, mp.mpf("0.5")],
        )
        assert abs(direct - integrated) < mp.mpf("1e-45")


@pytest.mark.parametrize("power,expected", ((2, s.pi / 4), (5, 5 * s.pi / 256)))
def test_complete_infinite_radial_integrals(power, expected):
    with mp.workdps(60):
        actual = mp.quad(lambda y: y * y / (1 + y * y) ** power, [0, 1, mp.inf])
        assert abs(actual - mp.mpf(str(expected.evalf(60)))) < mp.mpf("1e-55")


@pytest.mark.parametrize("fixture", range(8))
def test_arbitrary_external_momentum_majorants(fixture):
    rng = np.random.default_rng(13400 + fixture)
    m = 1000
    A = 25 / 16
    for scale in (1, 1e3, 1e7):
        k = scale * rng.normal(size=3)
        P = scale * rng.normal(size=3)
        l = P - k
        nu = np.sqrt(m * m + np.dot(k, k) / A**2)
        mu = np.sqrt(m * m + np.dot(l, l) / A**2)
        L = 1 + np.linalg.norm(P) / (A * m)
        assert mu <= nu * L * (1 + 1e-14)
        assert nu <= mu * L * (1 + 1e-14)
        assert nu * mu / (nu + mu) ** 6 <= 1 / (4 * nu**4) * (1 + 1e-14)
        error = nu * mu * (nu**-6 + mu**-6) ** 2
        assert error <= 2 * L * (nu**-10 + mu**-10) * (1 + 1e-14)
        assert L <= 2 * (1 + np.dot(P, P))


@pytest.mark.parametrize("norm", (1, 2, s.Rational(3, 7), s.Rational(5, 3)))
def test_exact_gravity_normalizations_and_not_a_solution(norm):
    p = noise.data()
    total = p["complete_stress_variance_coefficient_upper"]
    K = s.Integer(10) ** 800
    assert total * norm**2 < (10**25 * norm) ** 2
    assert total * norm**2 / K**2 < (s.Rational(1, 10**775) * norm) ** 2
    assert total * norm**2 / K < (s.Rational(1, 10**375) * norm) ** 2
    assert (
        "not yet a fully constraint-reduced"
        in p["Einstein_normalized_metric_force_standard_deviation_display"]
    )
    assert audit.controls()["noise_input_not_a_metric_stability_certificate"]
    assert "conjugate(u_r" in p["exact_Wick_pair_formula"]
