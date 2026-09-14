"""Independent non-aliased Fourier, cotangent and reference-symmetry diagnostics."""

import numpy as np
import pytest
import sympy as s
from numpy.testing import assert_allclose
from p8_vacuum_affine_nonlinear_spatial_reduction import (
    audit,
    cotangent,
    shape,
    source,
    translations,
)


@pytest.mark.parametrize(
    "packet", [source.data, shape.data, cotangent.data, translations.data]
)
def test_complete_packet_residuals_and_proof_gates(packet):
    data = packet()
    for name, value in data["checks"].items():
        assert all(
            x == 0
            for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
        ), name
    assert all(bool(v) for v in data["gates"].values())


def determinant(A):
    return (
        A[..., 0, 0] * (A[..., 1, 1] * A[..., 2, 2] - A[..., 1, 2] * A[..., 2, 1])
        - A[..., 0, 1] * (A[..., 1, 0] * A[..., 2, 2] - A[..., 1, 2] * A[..., 2, 0])
        + A[..., 0, 2] * (A[..., 1, 0] * A[..., 2, 1] - A[..., 1, 1] * A[..., 2, 0])
    )


def fourier_fixture(cutoff, amplitude):
    # Degree3 products have support<=3K, below the4K Nyquist surface.
    n = 8 * cutoff
    angle = np.arange(n) * 2 * np.pi / n
    xyz = np.meshgrid(angle, angle, angle, indexing="ij")
    axes = np.fft.fftfreq(n) * n
    k = np.stack(np.meshgrid(axes, axes, axes, indexing="ij"), axis=-1)
    square = np.sum(k * k, axis=-1)
    inv = np.divide(1.0, square, out=np.zeros_like(square), where=square > 0)
    P = np.eye(3) - k[..., :, None] * k[..., None, :] * inv[..., None, None]
    L = 2 + (square == 0)
    mask = np.max(np.abs(k), axis=-1) <= cutoff
    tau = np.zeros((n, n, n, 3, 3))
    for direction, (i, j) in enumerate([(1, 2), (0, 2), (0, 1)]):
        tau[..., i, j] = tau[..., j, i] = amplitude * np.cos(xyz[direction]) / 3

    def fft(x):
        return np.fft.fftn(x, axes=(0, 1, 2), norm="forward")

    def ifft(x):
        return np.fft.ifftn(x, axes=(0, 1, 2), norm="forward")

    def B(f):
        result = ifft(P * f[..., None, None])
        assert np.max(np.abs(result.imag)) < 1e-13
        return result.real

    def nonlinear(f):
        A = tau + B(f)
        tr = np.trace(A, axis1=-2, axis2=-1)
        return (
            -0.5 * fft(tr * tr - np.einsum("...ij,...ji->...", A, A)) / L
            - fft(determinant(A)) / L
        )

    f = np.zeros((n, n, n), complex)
    for _ in range(40):
        nxt = mask * nonlinear(f)
        change = np.sum(np.abs(nxt - f))
        f = nxt
        if change < 1e-18:
            break
    Q = np.eye(3) + tau + B(f)
    f2 = 0.5 * fft(np.einsum("...ij,...ji->...", tau, tau)) / L
    f3 = (fft(np.einsum("...ij,...ji->...", tau, B(f2))) - fft(determinant(tau))) / L
    return {
        "f": f,
        "Q": Q,
        "tau": tau,
        "B": B,
        "fft": fft,
        "k": k,
        "L": L,
        "map": nonlinear,
        "defect": np.sum(np.abs(f - nonlinear(f))),
        "f2": f2,
        "f3": f3,
        "amplitude": amplitude,
    }


@pytest.mark.parametrize("cutoff", [1, 2, 3, 4])
def test_full_nonaliased_projected_shape_with_generated_mean(cutoff):
    data = fourier_fixture(cutoff, 0.01)
    Q = data["Q"]
    assert abs(data["f"][0, 0, 0]) > 1e-7
    assert data["L"][0, 0, 0] == 3
    assert np.linalg.eigvalsh(Q).min() > 0.98
    divergence = np.sum(data["k"][..., None, :] * data["fft"](Q), axis=-1)
    assert np.max(np.abs(divergence)) < 1e-13
    assert data["defect"] < 0.001
    if cutoff == 4:
        assert data["defect"] < 1e-13
        assert np.max(np.abs(determinant(Q) - 1)) < 1e-11
    if cutoff == 1:
        # A projected fixed point is not the exact full determinant equation.
        assert data["defect"] > 1e-7
        assert np.max(np.abs(determinant(Q) - 1)) > 1e-7


@pytest.mark.parametrize("amplitude", [0.01, 0.007, 0.003, 0.001])
def test_full_cubic_shape_contacts_and_fourth_order_remainder(amplitude):
    data = fourier_fixture(4, amplitude)
    error = np.sum(np.abs(data["f"] - data["f2"] - data["f3"]))
    assert error < 60 * amplitude**4
    assert np.sum(np.abs(data["f2"])) < 0.75 * amplitude**2
    assert np.sum(np.abs(data["f3"])) < 4.125 * amplitude**3
    # Translation covariance checks the FULL nonlinear map, including the mean.
    shift = np.array([0.31, -0.17, 0.23])
    character = np.exp(1j * np.einsum("...i,i->...", data["k"], shift))
    assert character[0, 0, 0] == 1
    shifted = data["f"] * character
    shifted_Q = np.fft.ifftn(
        data["fft"](data["Q"]) * character[..., None, None],
        axes=(0, 1, 2),
        norm="forward",
    )
    assert np.max(np.abs(shifted_Q.imag)) < 1e-12
    assert np.max(np.abs(determinant(shifted_Q.real) - 1)) < 1e-10
    assert abs(shifted[0, 0, 0] - data["f"][0, 0, 0]) < 1e-16


@pytest.mark.parametrize("case", range(12))
def test_full_metric_shape_derivative_and_adjoint(case):
    rng = np.random.default_rng(2670 + case)
    C = rng.normal(size=(3, 3))
    gamma = np.eye(3) + 0.02 * (C + C.T)
    h = rng.normal(size=(3, 3))
    h = (h + h.T) / 2
    b = rng.normal(size=(3, 3))
    b = (b + b.T) / 2
    inv = np.linalg.inv(gamma)
    c = determinant(gamma) ** (1 / 3)
    dq = c * (np.trace(inv @ h) * inv / 3 - inv @ h @ inv)
    adj = c * (inv * np.trace(b @ inv) / 3 - inv @ b @ inv)
    assert_allclose(np.trace(b @ dq), np.trace(adj @ h), atol=1e-13)

    def Q(g):
        return determinant(g) ** (1 / 3) * np.linalg.inv(g)

    step = 1e-5
    finite = (
        -Q(gamma + 2 * step * h)
        + 8 * Q(gamma + step * h)
        - 8 * Q(gamma - step * h)
        + Q(gamma - 2 * step * h)
    ) / (12 * step)
    assert_allclose(dq, finite, atol=3e-10, rtol=3e-9)
    assert abs(np.trace(np.linalg.inv(Q(gamma)) @ dq)) < 1e-13


@pytest.mark.parametrize("case", range(12))
def test_complete_nonsymmetric_horizontal_and_dual_lift(case):
    rng = np.random.default_rng(26700 + case)
    n, g = 20, 3
    A = rng.normal(size=(g, n))
    V = A.T + 0.04 * rng.normal(size=(n, g))
    M = A @ V
    H = np.eye(n) - V @ np.linalg.solve(M, A)
    alpha = rng.normal(size=n)
    lifted = alpha - A.T @ np.linalg.solve(M.T, V.T @ alpha)
    assert_allclose(A @ H, 0, atol=4e-14)
    assert_allclose(H @ H, H, atol=4e-14)
    assert_allclose(lifted, H.T @ alpha, atol=4e-14)
    assert_allclose(V.T @ lifted, 0, atol=4e-14)
    tangent = H @ rng.normal(size=n)
    assert_allclose(lifted @ tangent, alpha @ tangent, atol=4e-14)
    C = rng.normal(size=(g, g))
    C = C - C.T
    block = np.block([[np.zeros((g, g)), M], [-M.T, C]])
    invM = np.linalg.inv(M)
    invblock = np.block([[invM.T @ C @ invM, -invM.T], [invM, np.zeros((g, g))]])
    assert_allclose(block @ invblock, np.eye(2 * g), atol=5e-14)


@pytest.mark.parametrize("n", [1, 2, 3, 8])
@pytest.mark.parametrize("case", range(6))
def test_full_coupled_pure_Gaussian_translation_action(n, case):
    rng = np.random.default_rng(267000 + 10 * n + case)
    S = np.block([[np.zeros((n, n)), np.eye(n)], [-np.eye(n), np.zeros((n, n))]])
    raw = rng.normal(size=(2 * n, 2 * n)) + 1j * rng.normal(size=(2 * n, 2 * n))
    raw = (raw + raw.T) / 2
    Z = 3 * np.eye(2 * n) + 0.02 * (raw + S.T @ raw @ S)
    A, B = Z.real, Z.imag
    assert np.linalg.eigvalsh(A).min() > 1
    assert_allclose(Z @ S, S @ Z, atol=1e-14)
    Ai = np.linalg.inv(A)
    covariance = 0.5 * np.block([[Ai, -Ai @ B], [-B @ Ai, A + B @ Ai @ B]])
    I = np.eye(2 * n)
    O = np.zeros_like(I)
    Omega = np.block([[O, I], [-I, O]])
    assert_allclose(covariance @ Omega @ covariance, Omega / 4, atol=3e-14)
    q = rng.normal(size=2 * n)
    theta = 0.271
    rotated = (np.cos(theta) * np.eye(2 * n) + np.sin(theta) * S) @ q
    assert abs((S @ q) @ Z @ q) < 1e-13
    assert_allclose(rotated @ Z @ rotated, q @ Z @ q, atol=1e-13)
    fullS = np.block([[S, O], [O, S]])
    assert_allclose(fullS @ covariance - covariance @ fullS, 0, atol=1e-13)


@pytest.mark.parametrize("sign", [-1, 1])
def test_rotation_invariant_mixture_need_not_have_zero_charge_support(sign):
    a, b = s.symbols("mixed_mode_a mixed_mode_b", real=True)
    polynomial = a + sign * s.I * b
    charge = -s.I * (b * s.diff(polynomial, a) - a * s.diff(polynomial, b))
    assert s.expand(charge + sign * polynomial) == 0
    # The radial Gaussian factor is invariant; the two orthogonal states
    # have opposite unit charges and their half-mixture has variance1.


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[x[0] for x in audit.bad_cases()]
)
def test_rejected_scope_and_domain_inputs(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_all_original_parameters_and_frontier_retained():
    assert audit.require_parameters(audit.parameters()) == audit.parameters()
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 123
    assert audit.require_norm(0) == 0
    assert audit.require_norm(shape.TAU_RADIUS) == shape.TAU_RADIUS
    assert audit.require_radius(s.Rational(3, 2)) == s.Rational(3, 2)
    assert audit.require_cutoff(4) == 4
