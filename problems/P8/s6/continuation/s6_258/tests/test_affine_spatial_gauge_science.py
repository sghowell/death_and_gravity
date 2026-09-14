"""Exact whole identities and independent nonlinear spatial-gauge diagnostics."""

import math

import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_spatial_gauge import audit, gauge, ghost, spatial
from scipy.linalg import det, eigvalsh, expm
from scipy.special import iv


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_whole_exact_identity(name):
    value = audit.residuals()[name]
    assert all(
        entry == 0
        for entry in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_every_domain_and_scope_rejection(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_whole_counts_fixed_parameters_and_unchanged_frontiers():
    assert (
        len(audit.residuals()),
        audit.scalar_entry_count(),
        len(audit.gates()),
        len(audit.controls()),
        audit.rejected_inputs(),
        len(audit.matching()),
    ) == (45, 1318, 40, 9, 312, 114)
    assert all(value is True for value in audit.gates().values())
    assert audit.require_parameters(audit.parameters()) == audit.parameters()
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()


@pytest.mark.parametrize("slot", range(13))
def test_independent_full_canonical_momentum_pairing(slot):
    data = spatial.generator()
    pi = data["whole_metric_and_metric_momenta"][1]
    pW = data["whole_vector_and_vector_momenta"][1]
    fields, ps = data["whole_two_matter_and_two_auxiliary_scalar_pairs"]
    x = data["whole_spatial_coordinates"]
    xi = data["whole_spatial_gauge_vector"]
    upper = [(i, j) for i in range(3) for j in range(i, 3)]
    variables = [pi[i, j] for i, j in upper] + list(pW) + list(ps)
    desired = (
        [(1 if i == j else 2) * data["whole_metric_Lie_action"][i, j] for i, j in upper]
        + list(data["whole_vector_Lie_action"])
        + [sum(xi[i] * s.diff(field, x[i]) for i in range(3)) for field in fields]
    )
    assert (
        s.expand(
            s.diff(data["whole_canonical_Lie_pairing_density"], variables[slot])
            - desired[slot]
        )
        == 0
    )


def test_boundary_flux_cannot_be_dropped_from_pointwise_generator():
    data = spatial.generator()
    flux = data["whole_spatial_generator_boundary_flux"]
    x = data["whole_spatial_coordinates"]
    remainder = s.expand(
        data["whole_canonical_Lie_pairing_density"]
        - data["whole_spatial_gauge_vector"].dot(
            data["whole_spatial_generator_density"]
        )
    )
    assert remainder != 0
    assert s.expand(remainder - sum(s.diff(flux[i], x[i]) for i in range(3))) == 0


def test_independent_spatially_varying_conformal_and_volume_frame():
    x, y, z = s.symbols("test_x test_y test_z", real=True)
    v, shape, logC = x * y + z, x * x + y * z, x * z + y * y
    gamma = s.diag(s.exp(2 * v + 2 * shape), s.exp(2 * v - 2 * shape), s.exp(2 * v))
    physical = s.exp(logC) * gamma
    Q = (s.exp(6 * v) ** s.Rational(1, 3) * gamma.inv()).applyfunc(s.simplify)
    Qphys = (s.exp(3 * logC + 6 * v) ** s.Rational(1, 3) * physical.inv()).applyfunc(
        s.simplify
    )
    assert Qphys == Q == s.diag(s.exp(-2 * shape), s.exp(2 * shape), 1)
    assert any(s.diff(logC, a) != 0 for a in (x, y, z))
    assert all(
        s.simplify(
            sum(s.diff(Qphys[i, j] - Q[i, j], a) for j, a in enumerate((x, y, z)))
        )
        == 0
        for i in range(3)
    )


@pytest.mark.parametrize("mode", ((1, 0, 0), (0, 2, -1), (1, -2, 3), (-3, 4, 2)))
def test_independent_nonzero_mode_TT_projection_and_volume_contact(mode):
    k = s.Matrix(mode)
    h = s.Matrix(
        [
            [s.Rational(2, 7), s.Rational(1, 11), s.Rational(-3, 13)],
            [s.Rational(1, 11), s.Rational(-1, 5), s.Rational(4, 17)],
            [s.Rational(-3, 13), s.Rational(4, 17), s.Rational(-3, 35)],
        ]
    )
    assert s.trace(h) == 0
    M = k.dot(k) * s.eye(3) + k * k.T / 3
    xi = M.inv() * (s.I * h * k)
    delta = s.I * (k * xi.T + xi * k.T - s.Rational(2, 3) * s.eye(3) * k.dot(xi))
    P = s.eye(3) - k * k.T / k.dot(k)
    TT = P * h * P - P * s.trace(P * h) / 2
    assert h + delta == TT
    assert TT * k == s.zeros(3, 1) and s.trace(TT) == 0
    assert s.I * k.dot(xi) / 3 == -(k.T * h * k)[0] / (4 * k.dot(k))
    assert M.inv() * (s.I * TT * k) == s.zeros(3, 1)
    assert audit.require_mode(mode) == mode


@pytest.mark.parametrize("seed", (11, 23, 47, 97))
def test_independent_nonsymmetric_principal_and_full_weak_coercivity(seed):
    rng = np.random.default_rng(seed)
    raw = rng.normal(size=(3, 3))
    symmetric = (raw + raw.T) / 2
    perturb = symmetric / (8 * np.linalg.norm(symmetric, ord=2))
    Q = np.eye(3) + perturb
    k = rng.normal(size=3)
    radial = k @ Q @ k
    principal = radial * np.eye(3) + np.outer(Q @ k, k) / 3
    symmetrized = np.linalg.solve(Q, principal)
    assert np.linalg.norm(symmetrized - symmetrized.T) < 1e-12
    assert min(eigvalsh(symmetrized, check_finite=True)) > 0
    assert abs(det(principal, check_finite=True) - (4 / 3) * radial**3) < 1e-9 * max(
        1, radial**3
    )
    G = rng.normal(size=(3, 3))
    energy = np.einsum("ij,jk,ik->", G, Q, G) + np.sum(G * Q) * np.trace(G) / 3
    assert energy >= 0.75 * np.sum(G * G) - 1e-12


def full_bessel_matrix(t, modes):
    """Independent entire exponential density; no Taylor replacement."""
    out = np.zeros((3 * len(modes), 3 * len(modes)))
    for i, pp in enumerate(modes):
        p = np.asarray(pp, dtype=float)
        for j, kk in enumerate(modes):
            k = np.asarray(kk, dtype=float)
            ell = p - k
            if ell[1] != 0 or ell[2] != 0:
                continue
            order = int(abs(ell[0]))
            Q = np.diag(
                [1.0 if order == 0 else 0.0, (-1) ** order * iv(order, t), iv(order, t)]
            )
            # This is the Fourier transform of the full density Lie derivative.
            out[3 * i : 3 * i + 3, 3 * j : 3 * j + 3] = (
                (p @ Q @ k) * np.eye(3)
                + np.outer(Q @ k, p)
                - np.outer(Q @ p, ell)
                - (2 / 3) * np.outer(Q @ p, k)
            )
    return out


@pytest.mark.parametrize("step", (0.001, 0.002, 0.003))
def test_independent_entire_Bessel_ghost_matrix_and_logdet_vertices(step):
    packet = ghost.finite_vertices()
    modes = packet["explicit_finite_nonzero_Fourier_modes"]
    base, plus, minus = (full_bessel_matrix(t, modes) for t in (0, step, -step))
    assert base.shape == (24, 24)
    assert (
        np.linalg.norm(
            base - np.asarray(packet["whole_three_ghost_reference_matrix"], dtype=float)
        )
        < 1e-13
    )
    first, second = (
        (plus - minus) / (2 * step),
        (plus + minus - 2 * base) / (step * step),
    )
    assert (
        np.linalg.norm(
            first - np.asarray(packet["whole_first_shape_ghost_matrix"], dtype=float)
        )
        < 3e-5
    )
    assert (
        np.linalg.norm(
            second - np.asarray(packet["whole_second_shape_ghost_matrix"], dtype=float)
        )
        < 3e-5
    )
    reference = det(base, check_finite=True)
    second_log = (
        math.log(det(plus, check_finite=True) / reference)
        + math.log(det(minus, check_finite=True) / reference)
    ) / (step * step)
    assert abs(second_log - 3) < 3e-5
    assert (
        abs(
            second_log
            - float(packet["incorrect_second_variation_without_shape_contact"])
        )
        > 0.1
    )


@pytest.mark.parametrize("cutoff", (1, 2))
@pytest.mark.parametrize("t", (-0.01, 0.01))
def test_full_three_dimensional_Galerkin_weak_coercivity(cutoff, t):
    modes = tuple(
        (a, b, c)
        for a in range(-cutoff, cutoff + 1)
        for b in range(-1, 2)
        for c in range(-1, 2)
        if (a, b, c) != (0, 0, 0)
    )
    matrix = full_bessel_matrix(t, modes)
    kinetic = np.repeat([sum(k * k for k in mode) for mode in modes], 3)
    normalized = matrix / np.sqrt(kinetic[:, None] * kinetic[None, :])
    # A real coefficient family makes this also the Hermitian part in complex modes.
    assert min(eigvalsh((normalized + normalized.T) / 2, check_finite=True)) >= 0.75
    assert matrix.shape[0] == 3 * len(modes)
    assert (0, 0, 0) not in modes


@pytest.mark.parametrize("L", (s.Rational(1, 7), 1, 10**30))
def test_retained_finite_torus_IR_scale(L):
    data = gauge.local_slice()
    symbol = next(iter(data["mean_zero_torus_L2_gap_lower_bound"].free_symbols))
    assert (
        data["mean_zero_torus_L2_gap_lower_bound"].subs(symbol, L)
        == s.Rational(3, 4) / L**2
    )
    assert (
        data["mean_zero_torus_inverse_L2_upper_bound"].subs(symbol, L)
        == s.Rational(4, 3) * L**2
    )
    assert audit.require_torus_radius(L) == L


def test_translation_kernel_is_surface_only_and_complement_is_not_algebra():
    Q = s.Matrix([[1, s.Rational(1, 7), 0], [s.Rational(1, 7), 2, 0], [0, 0, 3]])
    ell = s.Matrix([1, 2, 0])
    assert ghost.spectral_block(Q, (0, 0, 0), ell) == -(Q * ell) * ell.T
    assert ghost.spectral_block(Q, (0, 0, 0), ell) != s.zeros(3)
    assert ghost.spectral_block(Q, ell, -ell) == s.zeros(3)
    zero = gauge.zero_modes()
    assert zero["whole_nonzero_constant_gauge_bracket"] == -1
    assert zero["full_finite_SBP_Leibniz_defect"] != s.zeros(3, 1)


def test_independent_second_order_gauge_restoration_with_volume_contact():
    data = ghost.mixed_off_gauge_contact()
    ell = s.Matrix([1, 1, 0])
    Q = data["whole_mixed_exponential_shape_Fourier_contact"]
    xi = (ell.dot(ell) * s.eye(3) + ell * ell.T / 3).inv() * (-s.I * Q * ell)
    restored = Q - s.I * (
        ell * xi.T + xi * ell.T - s.Rational(2, 3) * s.eye(3) * ell.dot(xi)
    )
    assert restored * ell == s.zeros(3, 1)
    assert xi == data["whole_second_order_gauge_restoring_vector"]
    assert s.I * ell.dot(xi) / 3 == s.Rational(1, 32)
    assert Q * ell != s.zeros(3, 1)
    for k in ((1, -1, 1), (0, 2, 3), (-2, 1, 0)):
        k = s.Matrix(k)
        wanted = (k.T * restored * k)[0] * s.eye(3) + (restored * k) * k.T / 3
        assert ghost.spectral_block(restored, k, ell) == wanted


@pytest.mark.parametrize("step", (0.0005, 0.001))
def test_independent_full_exponential_mixed_Fourier_contact(step):
    data = ghost.mixed_off_gauge_contact()
    A, B = (
        np.asarray(item, dtype=float)
        for item in data["two_distinct_TT_shape_polarizations"]
    )
    coefficient = np.zeros((3, 3), dtype=complex)
    grid = 12
    for i in range(grid):
        x = 2 * math.pi * i / grid
        for j in range(grid):
            y = 2 * math.pi * j / grid
            a, b = A * math.cos(x), B * math.cos(y)
            mixed = (
                expm(-step * a - step * b)
                - expm(-step * a + step * b)
                - expm(step * a - step * b)
                + expm(step * a + step * b)
            ) / (4 * step * step)
            coefficient += mixed * np.exp(-1j * (x + y)) / (grid * grid)
    expected = np.asarray(
        data["whole_mixed_exponential_shape_Fourier_contact"], dtype=float
    )
    assert np.linalg.norm(coefficient - expected) < 2e-7
    assert np.linalg.norm(coefficient.imag) < 1e-10


def test_same_parent_fixed_reference_and_explicit_nonclosure():
    scope = audit.observable()
    assert "OPEN" in scope["original_problem"]
    assert "fixed quantum mean" in scope["scope_boundary"]
    assert audit.require_state(audit.STATE) == audit.STATE
    assert audit.require_gauge(audit.GAUGE) == audit.GAUGE
    assert audit.require_shape_distance(gauge.SHAPE_BOUND) == s.Rational(1, 8)
    assert "not simply div t=0" in spatial.frame()["fixed_reference_boundary"]
