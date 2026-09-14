"""Independent adjoint, full-covariance, normalization and cutoff diagnostics."""

from functools import cache

import mpmath as mp
import numpy as np
import pytest
import sympy as s
from numpy.testing import assert_allclose
from p8_vacuum_affine_nonlinear_spatial_reduction.cotangent import dq_adjoint
from p8_vacuum_affine_quantitative_phase_domain import (
    audit,
    geometry,
    invariants,
    reference,
    source,
)
from scipy.integrate import quad
from scipy.linalg import block_diag, expm
from scipy.special import expit, gammaincc, gammaln, iv


@pytest.mark.parametrize(
    "packet", [source.data, reference.data, geometry.data, invariants.data]
)
def test_complete_packets_and_written_proof_gates(packet):
    data = packet()
    for name, value in data["checks"].items():
        assert all(
            v == 0
            for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
        ), name
    assert all(bool(v) for v in data["gates"].values())


@cache
def full_Galerkin_ghost(cutoff, amplitude):
    # Exact Fourier coefficients of diag(1,exp(a cos x),exp(-a cos x)).
    # This is a diagnostic compression of the full smooth coefficient,
    # not a discrete diffeomorphism algebra or the full nonlinear theorem.
    modes = [
        (i, j, 0)
        for i in range(-cutoff, cutoff + 1)
        for j in range(-cutoff, cutoff + 1)
        if (i, j) != (0, 0)
    ]
    size = 3 * len(modes)
    M = np.zeros((size, size))
    full_adjoint = np.zeros_like(M)
    wrong_adjoint = np.zeros_like(M)
    flat = []
    for p in modes:
        p = np.asarray(p, dtype=float)
        flat.append(np.dot(p, p) * np.eye(3) + np.outer(p, p) / 3)
    for row, k in enumerate(modes):
        k = np.asarray(k, dtype=float)
        for column, p in enumerate(modes):
            p = np.asarray(p, dtype=float)
            difference = k - p
            if difference[1] != 0 or difference[2] != 0:
                continue
            n = int(difference[0])
            coefficient = np.diag(
                [
                    float(n == 0),
                    iv(abs(n), amplitude),
                    (-1.0) ** n * iv(abs(n), amplitude),
                ]
            )
            first = np.dot(p, coefficient @ p)
            block = first * np.eye(3) + np.outer(coefficient @ p, p) / 3
            # All derivative-Q contacts combine into the OUTPUT wavevector.
            adjoint = (
                np.dot(k, coefficient @ k) * np.eye(3)
                + np.outer(k, coefficient @ k) / 3
            )
            wrong = first * np.eye(3) + np.outer(p, coefficient @ p) / 3
            a, b = slice(3 * row, 3 * row + 3), slice(3 * column, 3 * column + 3)
            M[a, b], full_adjoint[a, b], wrong_adjoint[a, b] = block, adjoint, wrong
    return M, full_adjoint, wrong_adjoint, block_diag(*flat)


@pytest.mark.parametrize("cutoff", [1, 2, 3])
@pytest.mark.parametrize("amplitude", [5e-5, 1e-4, 2e-4])
def test_full_variable_coefficient_adjoint_and_Neumann_inverse(cutoff, amplitude):
    M, adjoint, wrong, flat = full_Galerkin_ghost(cutoff, amplitude)
    assert_allclose(adjoint, M.T, atol=5e-15)
    assert np.linalg.norm(wrong - M.T, 2) > amplitude / 10
    assert np.linalg.norm(M - M.T, 2) > amplitude / 10
    flat_inverse = np.linalg.inv(flat)
    perturbation = (adjoint - flat) @ flat_inverse
    # Conservative whole-coefficient A2 majorant from the exponential series.
    budget = 100 * np.expm1(4 * amplitude)
    assert budget < 0.1
    assert np.linalg.norm(perturbation, 2) < budget
    partial = np.eye(len(M))
    term = np.eye(len(M))
    for _ in range(8):
        term = -perturbation @ term
        partial += term
    inverse = flat_inverse @ partial
    assert_allclose(adjoint @ inverse, np.eye(len(M)), atol=2e-13)
    rng = np.random.default_rng(269 + cutoff)
    left = rng.normal(size=len(M)) + 1j * rng.normal(size=len(M))
    right = rng.normal(size=len(M)) + 1j * rng.normal(size=len(M))
    assert_allclose(
        np.vdot(left, M @ right), np.vdot(adjoint @ left, right), atol=2e-11
    )


@pytest.mark.parametrize("case", range(12))
def test_full_nondiagonal_cotangent_trace_and_shear(case):
    rng = np.random.default_rng(26900 + case)
    raw = rng.normal(size=(3, 3))
    gamma = np.eye(3) + 0.02 * (raw + raw.T)
    assert np.linalg.eigvalsh(gamma).min() > 0.7
    B = rng.normal(size=(3, 3))
    B = (B + B.T) / 2
    # The independent exact function keeps all metric inverse and trace terms.
    gs = s.Matrix([[s.Rational(str(x)) for x in row] for row in gamma])
    bs = s.Matrix([[s.Rational(str(x)) for x in row] for row in B])
    correction = dq_adjoint(gs, bs)
    assert s.factor(s.trace(correction * gs)) == 0
    numerical = np.array(correction.evalf().tolist(), dtype=float)
    Pi_v = 0.173
    base = Pi_v * np.linalg.inv(gamma) / 6
    lifted = base + numerical
    assert_allclose(2 * np.trace(lifted @ gamma), Pi_v, atol=4e-14)
    assert np.linalg.norm(lifted - base) > 1e-3
    tf = lifted - np.trace(lifted @ gamma) * np.linalg.inv(gamma) / 3
    assert_allclose(np.trace(tf @ gamma), 0, atol=4e-14)


@pytest.mark.parametrize("case", range(8))
def test_entire_correlated_whitening_row_identity(case):
    rng = np.random.default_rng(269000 + case)
    omega = np.block([[np.zeros((4, 4)), np.eye(4)], [-np.eye(4), np.zeros((4, 4))]])
    raw = rng.normal(size=(8, 8))
    symplectic = expm(0.03 * omega @ (raw + raw.T))
    covariance = symplectic @ symplectic.T / 2
    assert_allclose(covariance @ omega @ covariance, omega / 4, atol=2e-14)
    row = rng.normal(size=8)
    w = rng.normal(size=8)
    w /= np.linalg.norm(w)
    assert_allclose(
        np.dot(row @ symplectic, row @ symplectic),
        2 * row @ covariance @ row,
        atol=2e-14,
    )
    assert abs(row @ symplectic @ w) <= np.sqrt(2 * row @ covariance @ row) * (
        1 + 1e-14
    )
    # Physical kappa factors stay exact rather than underflowing in this fixture.
    assert (1 / s.sqrt(reference.KAPPA)) ** 2 == 1 / reference.KAPPA
    assert (
        reference.field_bounds()["Pi_v_A0"] * s.sqrt(reference.KAPPA)
        > reference.safe_bounds()["Pi_v_A0"]
    )


def smooth_cutoff(radial, core, support):
    x = (radial - core) / (support - core)
    if x <= 0:
        return 1.0
    if x >= 1:
        return 0.0
    return expit(1 / x - 1 / (1 - x))


@pytest.mark.parametrize("dimension", [1, 2, 8, 48])
@pytest.mark.parametrize("multiple", [1.0, 3.0])
def test_same_core_positive_volume_cutoffs_and_full_gamma_mean(dimension, multiple):
    core = dimension * multiple
    epsilon = 0.01

    def difference(radial):
        chi1 = smooth_cutoff(radial, core, 4 * core)
        chi2 = smooth_cutoff(radial, core, 3 * core)
        volume_difference = epsilon * (-np.expm1(-radial))
        density = np.exp((dimension - 1) * np.log(radial) - radial - gammaln(dimension))
        return (chi1 - chi2) * volume_difference * density

    actual, error = quad(
        difference, core, 4 * core, epsabs=1e-16, epsrel=1e-12, limit=200
    )
    tail = gammaincc(dimension, core)
    assert abs(actual) <= 2 * epsilon * tail + 10 * error
    assert abs(actual) <= 2 * epsilon * np.sqrt(tail) + 10 * error
    assert actual > -10 * error


def test_actual_large_radius_tail_without_floating_underflow():
    with mp.workdps(160):
        t = mp.mpf(5) * mp.mpf(10) ** 39
        polynomial = sum(t**k / mp.factorial(k) for k in range(48))
        log_tail = -t + mp.log(polynomial)
        assert mp.isfinite(log_tail) and log_tail < -(mp.mpf(10) ** 39)
        actual = mp.exp(log_tail)
        bound = mp.exp(-(mp.mpf(10) ** 39))
        assert 0 < actual < bound
    assert s.Rational(8, 3) ** 100 > 10**40


@pytest.mark.parametrize("name", tuple(invariants.bounds()))
def test_every_complete_density_invariant_is_strictly_inside_original_box(name):
    value = invariants.bounds()[name]
    assert value > 0
    assert value < invariants.IMAGE < source.auxiliary.DELTA / 2


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[x[0] for x in audit.bad_cases()]
)
def test_unsupported_phase_state_time_and_closure_inputs_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_all_original_frontiers_and_exact_family_retained():
    assert audit.require_parameters(audit.parameters()) == audit.parameters()
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 125
    assert audit.require_radius(0) == 0
    assert audit.require_radius(reference.SUPPORT) == reference.SUPPORT
    assert audit.require_time(0) == 0
    assert audit.require_momentum(reference.P) == reference.P
    assert audit.require_dimension(48) == 48
