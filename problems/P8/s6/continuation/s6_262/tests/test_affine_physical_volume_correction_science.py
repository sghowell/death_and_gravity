"""Independent actual-coefficient, metric, source, quadrature and Weyl checks."""

import math

import numpy as np
import pytest
import sympy as s
from p8_affine import verify as certificate
from p8_vacuum_affine_coupled_gaussian_state import gaussian
from p8_vacuum_affine_gauge_mean_transport import crossing, modes
from p8_vacuum_affine_nonlinear_auxiliary_measure import canonical
from p8_vacuum_affine_physical_volume_correction import actual, audit, contact, source
from scipy import integrate, interpolate, linalg


@pytest.mark.parametrize("packet", tuple(audit.packets()))
def test_all_exact_packets(packet):
    data = audit.packets()[packet]
    assert certificate.certify_residuals(data["checks"])
    assert all(bool(value) for value in data["gates"].values())


@pytest.mark.parametrize(
    "name,call,args",
    audit.bad_cases(),
    ids=lambda value: value if isinstance(value, str) else None,
)
def test_all_unsupported_inputs_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_unchanged_frontier_and_explicit_erratum_priority():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 118
    assert audit.matching()[:-1] == audit.previous.matching()
    assert all(audit.gates().values())
    assert "REFUTED" in audit.observable()["explicit_refutation"]
    assert "OPEN" in audit.observable()["original_problem"]
    assert audit.require_bounce_slope(s.Rational(3, 2)) == s.Rational(3, 2)
    assert source.C == canonical.Cchi**2 == canonical.M**-2
    assert s.simplify(source.C ** s.Rational(3, 2) - canonical.U) == 0


@pytest.mark.parametrize("seed", range(4))
@pytest.mark.parametrize("Rvalue", [0.64, 1.0, 1.44, 2.25])
@pytest.mark.parametrize("Nvalue", [0.81, 1.17])
def test_entire_shifted_metric_scalar_and_Maxwell_by_direct_four_contractions(
    seed, Rvalue, Nvalue
):
    rng = np.random.default_rng(320 + seed)
    L = rng.normal(size=(3, 3)) / 4 + np.eye(3)
    gamma = L @ L.T + np.eye(3) / 2
    shift = rng.normal(size=3) / 3
    C = Rvalue**-0.5
    g = np.block(
        [
            [
                np.array([[Nvalue**2 - C * shift @ gamma @ shift]]),
                (-C * shift @ gamma)[None, :],
            ],
            [(-C * gamma @ shift)[:, None], -C * gamma],
        ]
    )
    gi = linalg.inv(g, check_finite=True)
    gamma_inv = linalg.inv(gamma, check_finite=True)
    eig = linalg.eigvalsh(g, check_finite=True)
    assert np.count_nonzero(eig > 0) == 1
    hat_volume = math.sqrt(linalg.det(gamma, check_finite=True))
    volume = math.sqrt(-linalg.det(g, check_finite=True))
    assert volume / hat_volume == pytest.approx(Nvalue * Rvalue**-0.75, rel=2e-13)
    scalar = rng.normal(size=4)
    normal = scalar[0] - shift @ scalar[1:]
    expected_scalar = Rvalue**-0.75 * normal**2 / (2 * Nvalue)
    expected_scalar -= (
        Nvalue * Rvalue**-0.25 * (scalar[1:] @ gamma_inv @ scalar[1:]) / 2
    )
    assert volume * (scalar @ gi @ scalar) / (2 * hat_volume) == pytest.approx(
        expected_scalar, rel=2e-12, abs=2e-12
    )
    F = rng.normal(size=(4, 4))
    F -= F.T
    assert np.count_nonzero(np.triu(F, 1)) == 6
    direct_Maxwell = -volume * np.einsum("ac,bd,ab,cd", gi, gi, F, F) / (4 * hat_volume)
    electric = F[0, 1:] - shift @ F[1:, 1:]
    EE = electric @ gamma_inv @ electric
    BB = np.einsum("ik,jl,ij,kl", gamma_inv, gamma_inv, F[1:, 1:], F[1:, 1:])
    expected_Maxwell = (
        Rvalue**-0.25 * EE / (2 * Nvalue) - Nvalue * Rvalue**0.25 * BB / 4
    )
    assert direct_Maxwell == pytest.approx(expected_Maxwell, rel=3e-12, abs=3e-12)
    packet = actual.metric()
    values = {
        source.R: Rvalue,
        source.N: Nvalue,
        **dict(zip(packet["whole_three_shift_components"], shift, strict=True)),
        **dict(
            zip(
                (
                    packet["whole_hat_spatial_metric"][i, j]
                    for i, j in ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
                ),
                (
                    gamma[0, 0],
                    gamma[0, 1],
                    gamma[0, 2],
                    gamma[1, 1],
                    gamma[1, 2],
                    gamma[2, 2],
                ),
                strict=True,
            )
        ),
    }
    np.testing.assert_allclose(
        np.array(packet["whole_physical_metric"].subs(values), float),
        g,
        rtol=2e-13,
        atol=2e-13,
    )
    np.testing.assert_allclose(
        np.array(packet["whole_physical_inverse_metric"].subs(values), float),
        gi,
        rtol=2e-12,
        atol=2e-12,
    )


@pytest.mark.parametrize("slope", [-3.0, 0.0, 1.5, 4.0])
def test_metric_field_Jacobian_with_lapse_dependent_C_by_independent_differences(slope):
    N0, C0 = 1.13, 0.83
    base = np.array([N0, 0.1, -0.2, 0.3, 1.4, 0.12, -0.08, 1.1, 0.15, 1.6])
    pairs = [(i, j) for i in range(4) for j in range(i, 4)]

    def metric_coordinates(point, physical):
        N, shift = point[0], point[1:4]
        z = point[4:]
        gamma = np.array([[z[0], z[1], z[2]], [z[1], z[3], z[4]], [z[2], z[4], z[5]]])
        C = C0 * math.exp(slope * (N - N0)) if physical else 1.0
        g = np.block(
            [
                [
                    np.array([[N * N - C * shift @ gamma @ shift]]),
                    (-C * shift @ gamma)[None, :],
                ],
                [(-C * gamma @ shift)[:, None], -C * gamma],
            ]
        )
        return np.array([g[i, j] for i, j in pairs])

    def jac(physical):
        step = 2e-5
        return np.column_stack(
            [
                (
                    metric_coordinates(base + step * np.eye(10)[i], physical)
                    - metric_coordinates(base - step * np.eye(10)[i], physical)
                )
                / (2 * step)
                for i in range(10)
            ]
        )

    ratio = linalg.det(jac(True), check_finite=True) / linalg.det(
        jac(False), check_finite=True
    )
    assert ratio == pytest.approx(C0**9, rel=2e-8, abs=2e-10)


@pytest.mark.parametrize("uvalue", [-s.Rational(1, 3), s.Integer(0), s.Rational(2, 5)])
def test_six_clock_families_against_independent_germ_derivatives(uvalue):
    N = source.N
    h = (1 + uvalue**2) ** 3
    Rgerm = 1 + (N**-2 - 1) / h
    # A comparison remainder tests precisely the finite-jet boundary.
    # It is NOT substituted for any original off-clock physical function.
    comparison = Rgerm + (N - 1) ** 8 * s.exp(N - 1)
    functions = {
        "C": comparison ** (-s.Rational(1, 2)),
        "spatial_volume_U": comparison ** (-s.Rational(3, 4)),
        "spacetime_volume_NU": N * comparison ** (-s.Rational(3, 4)),
        "scalar_temporal_U_over_N": comparison ** (-s.Rational(3, 4)) / N,
        "scalar_spatial_NCchi": N * comparison ** (-s.Rational(1, 4)),
        "Maxwell_magnetic_NM": N * comparison ** s.Rational(1, 4),
    }
    rows = source.clock_jets()["whole_six_actual_ADM_clock_families"]
    for name, expr in functions.items():
        for order in range(5):
            expected = s.diff(expr, N, order).subs(N, 1)
            assert (
                s.factor(rows[name][order].subs(source.parent.u, uvalue) - expected)
                == 0
            )
    assert comparison.subs(N, 2) != Rgerm.subs(N, 2)


@pytest.mark.parametrize("Rvalue", [0.64, 1.0, 1.44, 2.25])
def test_wrong_positive_C_fails_independent_action_coefficients(Rvalue):
    wrong = Rvalue - 0.5
    correct = Rvalue**-0.5
    if Rvalue == 1.0:
        assert wrong == 0.5 and correct == 1.0
    assert wrong**1.5 != pytest.approx(Rvalue**-0.75, rel=1e-10)
    assert math.sqrt(wrong) != pytest.approx(Rvalue**-0.25, rel=1e-10)
    assert 1 / math.sqrt(wrong) != pytest.approx(Rvalue**0.25, rel=1e-10)


@pytest.mark.parametrize(
    "a,v,N", [(1.2, 0.13, 0.93), (0.8, -0.17, 1.07), (1.0, 0.0, 1.0)]
)
def test_whole_nonlinear_source_gradient_and_Hessian_by_independent_differences(
    a, v, N
):
    # A nonconstant positive-R fixture checks both R_N and R_NN terms.
    R = lambda n: math.exp(0.2 * n + 0.3 * n * n)
    F = lambda z: a**3 * math.exp(3 * z[0]) * R(z[1]) ** -0.75
    point, step = np.array([v, N]), 2e-4
    basis = np.eye(2)
    gradient = np.array(
        [(F(point + step * q) - F(point - step * q)) / (2 * step) for q in basis]
    )
    hessian = np.array(
        [
            [
                (
                    F(point + step * q + step * r)
                    - F(point + step * q - step * r)
                    - F(point - step * q + step * r)
                    + F(point - step * q - step * r)
                )
                / (4 * step**2)
                for r in basis
            ]
            for q in basis
        ]
    )
    slope = -0.75 * (0.2 + 0.6 * N)
    second = slope**2 - 0.45
    expected_gradient = F(point) * np.array([3, slope])
    expected_hessian = F(point) * np.array([[9, 3 * slope], [3 * slope, second]])
    np.testing.assert_allclose(gradient, expected_gradient, rtol=4e-7, atol=2e-8)
    np.testing.assert_allclose(hessian, expected_hessian, rtol=8e-7, atol=3e-8)


@pytest.mark.parametrize("w", [-0.5, 2 / 3, 1.0, 1.4])
def test_nonlinear_periodic_volume_with_correct_R_power_by_inverse_quadrature(w):
    beta = (3 * w - 2) / (2 - w)
    y = np.linspace(0, 2 * math.pi, 16385)
    v = lambda z: 0.21 * np.cos(z) + 0.07 * np.sin(2 * z)
    # Diagnostic positive scalar R, not a replacement for the fixed parent.
    R = lambda z: np.exp(0.1 * np.sin(z) + 0.2 * np.cos(2 * z))
    weight = np.exp(-beta * v(y))
    K = np.trapezoid(weight, y) / (2 * math.pi)
    x_of_y = integrate.cumulative_trapezoid(weight / K, y, initial=0)
    x = np.linspace(0, 2 * math.pi, 8193)
    phi = interpolate.PchipInterpolator(x_of_y, y)(x)
    jac = K * np.exp(beta * v(phi))
    original = np.trapezoid(R(y) ** -0.75 * np.exp(3 * v(y)), y)
    pulled = np.trapezoid(R(phi) ** -0.75 * np.exp(3 * v(phi)) * jac, x)
    assert pulled == pytest.approx(original, rel=3e-7, abs=3e-8)
    for index in range(0, len(x), 512):
        z = phi[index]
        tensor, cross = 0.2 * math.cos(z), 0.13 * math.sin(z + 0.2)
        block = np.array(
            [[math.exp(tensor), cross], [cross, (1 + cross**2) * math.exp(-tensor)]]
        )
        hat_metric = math.exp(2 * v(z)) * linalg.block_diag(jac[index] ** 2, block)
        physical = R(z) ** -0.5 * hat_metric
        density = math.sqrt(linalg.det(physical, check_finite=True))
        assert density == pytest.approx(
            R(z) ** -0.75 * math.exp(3 * v(z)) * jac[index], rel=2e-12
        )


@pytest.mark.parametrize("seed", range(6))
@pytest.mark.parametrize("angle", [0.17, 1.23])
def test_correct_contacts_and_old_error_by_whole_correlated_Weyl_cubature(seed, angle):
    S, diagonal, M = gaussian.fixture_preparation(seed)
    assert (S.T * M * S - diagonal).applyfunc(s.factor) == s.zeros(4)
    transport = gaussian.symplectic_fixture((seed + 3) % 12)
    covariance = np.array(crossing.transported_covariance(M, transport), float)
    L = linalg.cholesky(
        linalg.block_diag(covariance, covariance), lower=True, check_finite=True
    )
    samples = np.concatenate((math.sqrt(8) * L.T, -math.sqrt(8) * L.T), axis=0)
    vr, nr = np.array([1.0, 0.0, 0.0, 0.0]), np.array([0.31, 0.17, 0.0, -0.43])
    vc, vs = samples[:, :4] @ vr, samples[:, 4:] @ vr
    nc, ns = samples[:, :4] @ nr, samples[:, 4:] @ nr
    norm = math.sqrt(2 / 17)
    v = norm * (vc * math.cos(angle) + vs * math.sin(angle))
    n = norm * (nc * math.cos(angle) + ns * math.sin(angle))
    dv_wave = norm * (-vc * math.sin(angle) + vs * math.cos(angle))
    dn_wave = norm * (-nc * math.sin(angle) + ns * math.cos(angle))
    dv2, dn2, dv1 = -2.25 * dv_wave**2, -2.25 * dv_wave * dn_wave, 0.75 * v
    Cvv, Cvn = 2 / 17 * (vr @ covariance @ vr), 2 / 17 * (vr @ covariance @ nr)
    one = np.mean(3 * dv2 + 1.5 * dn2)
    two = np.mean(9 * v * dv1 + 4.5 * n * dv1)
    wrong_one = np.mean(3 * dv2 - 6 * dn2)
    wrong_two = np.mean(9 * v * dv1 - 18 * n * dv1)
    assert one == pytest.approx(-27 * Cvv / 4 - 27 * Cvn / 8, abs=2e-11)
    assert one + two == pytest.approx(0, abs=2e-11)
    assert wrong_one + wrong_two == pytest.approx(0, abs=2e-11)
    assert wrong_one - one == pytest.approx(135 * Cvn / 8, abs=2e-11)


def test_corrected_all_generated_Fourier_modes_by_unaliased_sampling():
    data = modes.packet()
    points = np.stack(
        np.meshgrid(*([2 * math.pi * np.arange(11) / 11] * 3), indexing="ij"), axis=-1
    ).reshape(-1, 3)

    def values(field):
        return sum(
            complex(c) * np.exp(1j * (points @ np.array(k))) for k, c in field.items()
        )

    v = values(data["whole_three_direction_scalar_fixture"])
    n = values(data["whole_nonindependent_lapse_fixture"])
    dv1 = values(data["whole_first_log_volume_variation"])
    dv2 = values(data["whole_second_log_volume_variation"])
    dn2 = values(data["whole_second_lapse_variation"])
    one, two = np.mean(3 * dv2 + 1.5 * dn2), np.mean(9 * v * dv1 + 4.5 * n * dv1)
    corrected = contact.corrected_contact()
    assert one == pytest.approx(
        complex(corrected["whole_corrected_three_direction_finite_onepoint_contact"]),
        abs=1e-12,
    )
    assert two == pytest.approx(
        complex(corrected["whole_corrected_three_direction_finite_twopoint_contact"]),
        abs=1e-12,
    )
    assert one + two == pytest.approx(0, abs=1e-12)
    assert abs(one) > 0 and abs(two) > 0


def test_explicit_refutation_retains_covariance_but_not_wrong_physical_label():
    erratum = contact.retained_symbols_and_erratum()
    old = crossing.bounce_symbols()
    assert erratum["retained_whole_bounce_reference_symbols"] == {
        key: value for key, value in old.items() if key not in ("checks", "gates")
    }
    details = erratum["explicit_erratum"]
    assert "R**(-1/2)" in details["correct_binding"]
    assert "false" in details["refuted_identification"]
    assert "not renewed physical certification" in details["immutability"]
    assert source.clock_jets()[
        "whole_actual_volume_lapse_slope_at_bounce"
    ] == s.Rational(3, 2)
