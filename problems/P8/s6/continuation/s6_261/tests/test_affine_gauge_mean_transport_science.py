"""Independent quadrature, full Gaussian cubature, Fourier and matrix tests."""

import math

import numpy as np
import pytest
import sympy as s
from p8_affine import verify as certificate
from p8_vacuum_affine_coupled_gaussian_state import gaussian, phase
from p8_vacuum_affine_gauge_mean_transport import audit, crossing, geometry, modes
from scipy import integrate, interpolate, linalg


@pytest.mark.parametrize("packet", tuple(audit.packets()))
def test_all_exact_packets(packet):
    data = audit.packets()[packet]
    assert certificate.certify_residuals(data["checks"])
    assert all(bool(v) for v in data["gates"].values())


@pytest.mark.parametrize(
    "name,call,args",
    audit.bad_cases(),
    ids=lambda value: value if isinstance(value, str) else None,
)
def test_all_unsupported_inputs_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_full_unchanged_frontier_and_fixed_scope():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 117
    assert audit.matching()[:-1] == audit.previous.matching()
    assert all(audit.gates().values())
    assert "OPEN" in audit.observable()["original_problem"]
    assert audit.require_weight(s.Rational(2, 3)) == s.Rational(2, 3)
    assert audit.require_band(s.Rational(1, 3), 2) == (s.Rational(1, 3), 2)


@pytest.mark.parametrize("w", [-0.5, 2 / 3, 1.0, 1.4])
@pytest.mark.parametrize("amplitude", [0.07, 0.21, 0.4])
def test_exact_finite_periodic_map_by_independent_inverse_quadrature(w, amplitude):
    beta = (3 * w - 2) / (2 - w)
    y = np.linspace(0, 2 * math.pi, 16385)

    def v(z):
        return amplitude * np.cos(z) + amplitude * np.sin(2 * z) / 3

    def N(z):
        return 0.99 + 0.02 * np.cos(z + 0.3)

    def physical_C(z):
        return N(z) ** -2 - 0.5

    weighted = np.exp(-beta * v(y))
    K = np.trapezoid(weighted, y) / (2 * math.pi)
    x_of_y = integrate.cumulative_trapezoid(weighted / K, y, initial=0)
    x = np.linspace(0, 2 * math.pi, 8193)
    phi = interpolate.PchipInterpolator(x_of_y, y)(x)
    assert np.isfinite(phi).all() and np.min(np.diff(phi)) > 0
    jac = K * np.exp(beta * v(phi))
    new_v = v(phi) + np.log(jac) / 3
    actual_mean = np.trapezoid(new_v, x) / (2 * math.pi)
    expected_mean = (1 + beta / 3) * np.trapezoid(v(y) * weighted, y) / (
        2 * math.pi * K
    ) + math.log(K) / 3
    assert actual_mean == pytest.approx(expected_mean, rel=1e-6, abs=4e-8)
    actual_N = np.trapezoid(N(phi), x) / (2 * math.pi)
    expected_N = np.trapezoid(N(y) * weighted, y) / (2 * math.pi * K)
    assert actual_N == pytest.approx(expected_N, rel=1e-7, abs=4e-8)
    original_volume = np.trapezoid(physical_C(y) ** 1.5 * np.exp(3 * v(y)), y) / (
        2 * math.pi
    )
    pulled_volume = np.trapezoid(
        physical_C(phi) ** 1.5 * np.exp(3 * v(phi)) * jac, x
    ) / (2 * math.pi)
    assert pulled_volume == pytest.approx(original_volume, rel=4e-7, abs=4e-8)
    # Build the entire 3x3 metric, both transverse tensor profiles, and
    # its density inverse independently of the analytic x-density formula.
    for index in range(0, len(x), 512):
        yy = phi[index]
        tt = 0.2 * math.cos(yy)
        zz = 0.13 * math.sin(yy + 0.2)
        block = np.array([[math.exp(tt), zz], [zz, (1 + zz * zz) * math.exp(-tt)]])
        metric = linalg.block_diag(jac[index] ** 2, block) * math.exp(2 * v(yy))
        det = linalg.det(metric, check_finite=True)
        assert np.min(linalg.eigvalsh(metric, check_finite=True)) > 0
        Q = det ** (w / 2) * linalg.inv(metric, check_finite=True)
        assert Q[0, 0] == pytest.approx(K ** (w - 2), rel=2e-12, abs=2e-12)
        np.testing.assert_allclose(Q[0, 1:], 0, atol=1e-14)


@pytest.mark.parametrize("amplitude", [0.1, 0.3, 0.6])
def test_exact_mean_weight_slope_by_integrals(amplitude):
    def v(y):
        return amplitude * math.cos(y) + amplitude * math.sin(2 * y) / 3

    def n(y):
        return 1 + 0.03 * math.cos(y) + 0.02 * math.sin(3 * y)

    def avg(f):
        return integrate.quad(f, 0, 2 * math.pi, epsabs=2e-12, epsrel=2e-12)[0] / (
            2 * math.pi
        )

    def means(w):
        beta = (3 * w - 2) / (2 - w)
        K = avg(lambda y: math.exp(-beta * v(y)))
        return np.array(
            [
                (1 + beta / 3) * avg(lambda y: v(y) * math.exp(-beta * v(y))) / K
                + math.log(K) / 3,
                avg(lambda y: n(y) * math.exp(-beta * v(y))) / K,
            ]
        )

    step = 1e-5
    actual = (means(2 / 3 + step) - means(2 / 3 - step)) / (2 * step)
    expected = -2.25 * np.array(
        [
            avg(lambda y: v(y) ** 2) - avg(v) ** 2,
            avg(lambda y: v(y) * n(y)) - avg(v) * avg(n),
        ]
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-7, atol=2e-10)


def test_all_three_direction_products_against_exact_periodic_sampling():
    packet = modes.packet()
    N = 11
    points = np.stack(
        np.meshgrid(*([2 * math.pi * np.arange(N) / N] * 3), indexing="ij"), axis=-1
    ).reshape(-1, 3)

    def values(field):
        return (
            sum(
                complex(c) * np.exp(1j * (points @ np.array(k)))
                for k, c in field.items()
            )
            if field
            else np.zeros(len(points), complex)
        )

    v = values(packet["whole_three_direction_scalar_fixture"])
    n = values(packet["whole_nonindependent_lapse_fixture"])
    xi1 = np.array([values(f) for f in packet["whole_first_gauge_displacement"]])
    gv = np.array(
        [
            values(f)
            for f in modes.gradient(packet["whole_three_direction_scalar_fixture"])
        ]
    )
    gn = np.array(
        [
            values(f)
            for f in modes.gradient(packet["whole_nonindependent_lapse_fixture"])
        ]
    )
    np.testing.assert_allclose(
        np.mean(np.sum(xi1 * gv, axis=0)), -2.25 * np.mean(v * v), atol=1e-13
    )
    np.testing.assert_allclose(
        np.mean(np.sum(xi1 * gn, axis=0)), -2.25 * np.mean(v * n), atol=1e-13
    )
    assert np.mean(v * v).real == pytest.approx(
        float(packet["whole_finite_Cvv"]), abs=1e-13
    )
    assert np.mean(v * n).real == pytest.approx(
        float(packet["whole_finite_Cvn"]), abs=1e-13
    )
    for f in (
        *packet["whole_second_gauge_displacement"],
        packet["whole_second_log_volume_variation"],
    ):
        assert np.max(np.abs(values(f).imag)) < 1e-12


@pytest.mark.parametrize("seed", range(6))
@pytest.mark.parametrize("wave", [(1, 0, 0), (1, 2, -1), (0, 2, 3)])
@pytest.mark.parametrize("angle", [0.17, 1.23])
def test_fixed_entire_coupled_Weyl_moments_by_independent_Gaussian_cubature(
    seed, wave, angle
):
    S, diagonal, M = gaussian.fixture_preparation(seed)
    assert (S.T * M * S - diagonal).applyfunc(s.factor) == s.zeros(4)
    transport = gaussian.symplectic_fixture((seed + 3) % 12)
    V = crossing.transported_covariance(M, transport)
    expected_V = transport * S * S.T * transport.T / 2
    assert (V - expected_V).applyfunc(s.factor) == s.zeros(4)
    covariance = np.array(V, dtype=float)
    both = linalg.block_diag(covariance, covariance)
    L = linalg.cholesky(both, lower=True, check_finite=True)
    # Symmetric cubature with exact Gaussian moments through degree three.
    samples = np.concatenate((math.sqrt(8) * L.T, -math.sqrt(8) * L.T), axis=0)
    vr = np.array([1.0, 0.0, 0.0, 0.0])
    nr = np.array([0.31, 0.17, 0.0, -0.43])
    vc, vs = samples[:, :4] @ vr, samples[:, 4:] @ vr
    nc, ns = samples[:, :4] @ nr, samples[:, 4:] @ nr
    volume = 17.0
    norm = math.sqrt(2 / volume)
    v = norm * (vc * math.cos(angle) + vs * math.sin(angle))
    n = norm * (nc * math.cos(angle) + ns * math.sin(angle))
    k = np.array(wave, dtype=float)
    dv_wave = norm * (-vc * math.sin(angle) + vs * math.cos(angle))
    dn_wave = norm * (-nc * math.sin(angle) + ns * math.cos(angle))
    xi = -2.25 * dv_wave[:, None] * k[None, :] / (k @ k)
    grad_v, grad_n = dv_wave[:, None] * k[None, :], dn_wave[:, None] * k[None, :]
    dv2 = np.sum(xi * grad_v, axis=1)
    dn2 = np.sum(xi * grad_n, axis=1)
    Cvv = 2 / volume * (vr @ covariance @ vr)
    Cvn = 2 / volume * (vr @ covariance @ nr)
    assert np.mean(dv2) == pytest.approx(-2.25 * Cvv, rel=1e-12, abs=1e-12)
    assert np.mean(dn2) == pytest.approx(-2.25 * Cvn, rel=1e-12, abs=1e-12)
    dv1 = 0.75 * v
    contact = 3 * dv2 - 6 * dn2 + 9 * v * dv1 - 18 * n * dv1
    assert np.mean(contact) == pytest.approx(0, abs=1e-10)
    assert Cvv > 0 and np.max(np.abs(covariance - np.diag(np.diag(covariance)))) > 0


@pytest.mark.parametrize("theta", [s.Rational(0), s.Rational(1, 3), -s.Rational(2, 5)])
def test_reduced_lapse_commutator_by_independent_Fock_vacuum(theta):
    dim = 5
    annihilation = np.diag(np.sqrt(np.arange(1, dim)), 1)
    Q = (annihilation + annihilation.T) / math.sqrt(2)
    P = (annihilation - annihilation.T) / (1j * math.sqrt(2))
    I = np.eye(dim)
    fields = [np.kron(Q, I), np.kron(I, Q), np.kron(P, I), np.kron(I, P)]
    v = fields[0]
    J = 1.7
    n = (
        -0.8 * fields[0]
        + 0.3 * float(theta) * fields[1]
        + float(theta) * fields[2]
        - 0.05 * fields[3]
    ) / (2 * J)
    comm = v @ n - n @ v
    assert comm[0, 0] == pytest.approx(1j * float(theta) / (2 * J), abs=1e-13)
    if theta == 0:
        np.testing.assert_allclose(comm, 0, atol=1e-13)


@pytest.mark.parametrize("scale", [0.8, 1.0, 1.4])
@pytest.mark.parametrize("J", [1.51, 1.52, 1.7])
def test_bounce_symbols_against_full_noncommuting_positive_comparison_matrices(
    scale, J
):
    F, E, ell, T = 1199 / 800, -0.5, 0.1, 0.03
    speed = math.sqrt(F / J)
    Omega = np.array(phase.OMEGA, dtype=float)
    R = np.array([[E / math.sqrt(2 * J), 0], [-ell * E / math.sqrt(2 * J), 1]])
    invT = linalg.inv(R, check_finite=True).T
    errors = []
    for k in (64.0, 256.0, 1024.0):
        qq = np.diag([(k * speed / scale) ** 2, (k / scale) ** 2]) + np.array(
            [[1.0, 0.2], [0.2, 2.0]]
        )
        qp = np.array([[0.03, 0.05], [-0.07, 0.02]])
        pp = np.eye(2) + np.array([[0.1, 0.03], [0.03, 0.2]]) / k**2
        H = np.block([[qq, qp], [qp.T, pp]])
        assert np.min(linalg.eigvalsh(H, check_finite=True)) > 0
        determinant_root = math.sqrt(linalg.det(H, check_finite=True))
        V = (
            determinant_root * linalg.inv(H, check_finite=True) - Omega @ H @ Omega
        ) / (
            2
            * math.sqrt(-np.trace((Omega @ H) @ (Omega @ H)) / 2 + 2 * determinant_root)
        )
        q = k * k / scale**2
        vr = invT[0] / (2 * scale**1.5 * q)
        nr = (ell * E * invT[1] - (E + 3 * T / (2 * q)) * invT[0]) / (
            2 * J * scale**1.5
        )
        pp_cov = V[2:, 2:]
        measured = np.array(
            [k**3 * (vr @ pp_cov @ vr), k * (vr @ pp_cov @ nr), (nr @ pp_cov @ nr) / k]
        )
        expected = np.array(
            [
                (2 * J * speed / E**2 + ell**2) / 8,
                -speed / (4 * E * scale**2),
                speed / (4 * J * scale**4),
            ]
        )
        errors.append(np.max(np.abs(measured / expected - 1)))
    assert errors[-1] < 3e-6 and errors[1] < errors[0] / 8 and errors[2] < errors[1] / 8


@pytest.mark.parametrize(
    "bad", [s.zeros(4), s.eye(4) / 2, s.diag(1, 2, 1, 1), s.ones(4)]
)
def test_nonsymplectic_fixed_transport_rejected(bad):
    with pytest.raises(ValueError):
        crossing.transported_covariance(s.eye(4), bad)


def test_gauge_inverse_never_divides_constant_mode():
    with pytest.raises(ValueError, match="constant"):
        modes.inverse([{modes.ZERO: s.Integer(1)}, {}, {}])


def test_whole_mean_contacts_and_unchanged_quantum_boundary():
    packet = geometry.physical_volume_contact()
    assert packet["whole_linear_onepoint_volume_contact"] != 0
    assert packet["whole_quadratic_two_point_volume_contact"] != 0
    assert crossing.reference_rows()["whole_lapse_volume_equal_time_bracket"] != 0
    assert crossing.bounce_symbols()["whole_actual_bounce_F"] == s.Rational(1199, 800)
