"""Independent normalized Weyl kernels, full derivative and operator diagnostics."""

from functools import cache

import numpy as np
import pytest
import sympy as s
from numpy.polynomial.hermite import hermgauss
from numpy.testing import assert_allclose
from p8_vacuum_affine_weyl_operator_comparison import (
    audit,
    comparison,
    derivatives,
    kernel,
    source,
)
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm
from scipy.special import gammaln


@pytest.mark.parametrize(
    "packet", [source.data, kernel.data, derivatives.data, comparison.data]
)
def test_complete_packets_and_written_proof_gates(packet):
    data = packet()
    for name, value in data["checks"].items():
        assert all(
            x == 0
            for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
        ), name
    assert all(bool(g) for g in data["gates"].values())


def exact_Gaussian_Weyl_coherent_kernel(a, z, zprime):
    z, zprime = np.asarray(z), np.asarray(zprime)
    m = (z + zprime) / 2
    delta = z - zprime
    frequency = np.array([-delta[1], delta[0]])
    sigma = z[1] * zprime[0] - z[0] * zprime[1]
    linear = 2 * m + 1j * frequency
    # Dot, NOT conjugate-dot, in completing a complex Gaussian square.
    return np.exp(linear @ linear / (4 * (1 + a)) - m @ m + 0.5j * sigma) / (1 + a)


def defining_Weyl_kernel_quadrature(a, z, zprime, order):
    # Independent original x,y Weyl kernel after its exact xi Gaussian
    # integral. No coherent-kernel formula is used here.
    nodes, weights = hermgauss(order)
    x, y = np.meshgrid(nodes, nodes, indexing="ij")
    q, p = z
    qp, pp = zprime
    real = (
        -((x - q) ** 2) / 2
        - (y - qp) ** 2 / 2
        - a * (x + y) ** 2 / 4
        - (x - y) ** 2 / (4 * a)
        + x * x
        + y * y
    )
    phase = -p * (x - q / 2) + pp * (y - qp / 2)
    integrand = np.exp(real + 1j * phase) / (2 * np.pi * np.sqrt(a))
    return np.einsum("i,j,ij->", weights, weights, integrand)


@pytest.mark.parametrize("a", [0.3, 0.7, 1.2, 2.0])
@pytest.mark.parametrize("case", range(4))
def test_defining_Weyl_kernel_matches_independent_normalized_coherent_formula(a, case):
    rng = np.random.default_rng(2710 + case)
    z, zprime = rng.normal(size=(2, 2)) * 0.6
    expected = exact_Gaussian_Weyl_coherent_kernel(a, z, zprime)
    coarse = defining_Weyl_kernel_quadrature(a, z, zprime, 48)
    fine = defining_Weyl_kernel_quadrature(a, z, zprime, 72)
    assert_allclose(fine, expected, atol=2e-13, rtol=2e-12)
    assert abs(fine - expected) <= abs(coarse - expected) + 2e-13
    assert_allclose(
        exact_Gaussian_Weyl_coherent_kernel(a, zprime, z),
        np.conjugate(expected),
        atol=1e-14,
    )


@pytest.mark.parametrize("case", range(8))
def test_identity_symbol_frame_normalization_and_full_overlap_phase(case):
    rng = np.random.default_rng(27100 + case)
    z, zprime = rng.normal(size=(2, 2))
    delta = z - zprime
    sigma = z[1] * zprime[0] - z[0] * zprime[1]
    expected = np.exp(-delta @ delta / 4 - 0.5j * sigma)
    actual = exact_Gaussian_Weyl_coherent_kernel(0, z, zprime)
    assert_allclose(actual, expected, atol=2e-15)
    assert_allclose(exact_Gaussian_Weyl_coherent_kernel(0, z, z), 1, atol=2e-15)
    # The inconsistent printed extra factor2^d already fails at one pair.
    assert abs(2 * exact_Gaussian_Weyl_coherent_kernel(0, z, z) - 1) > 0.9


@pytest.mark.parametrize("a", [0.2, 0.8, 2.0])
@pytest.mark.parametrize("z", [np.array([0.0, 0.0]), np.array([0.3, -0.7])])
def test_actual_coherent_measure_Schur_integral_and_operator_bound(a, z):
    # Absolute kernel factors after Gaussian completion; integrate the
    # FULL real axes with coherent frame measure, not a finite phase grid.
    beta = (1 - a) / (2 * (1 + a))
    factors = []
    for component in z:
        value, error = quad(
            lambda v, component=component: np.exp(-v * v / 4 + beta * component * v),
            -np.inf,
            np.inf,
            epsabs=1e-12,
            epsrel=1e-12,
        )
        assert error < 1e-9
        factors.append(value)
    integral = np.prod(factors) * np.exp(-z @ z / 4) / (2 * np.pi * (1 + a))
    exact = 2 / (1 + a) * np.exp(-a * (z @ z) / (1 + a) ** 2)
    assert_allclose(integral, exact, atol=3e-12)
    assert integral <= 2 / (1 + a) + 3e-12
    assert 1 / (1 + a) < 2 / (1 + a)


@pytest.mark.parametrize("order", range(193))
def test_every_summed_product_Schur_multiindex_weight(order):
    weights = kernel.identities()["weights"]
    assert weights[order].is_Integer and weights[order] > 0
    assert len(weights) == 193
    if order == 0:
        assert weights[order] == 5**96
    if order == 192:
        assert weights[order] == 1
    assert sum(weights) == 10**96


@pytest.mark.parametrize("order", range(197))
def test_every_full_radial_partition_and_high_Cauchy_product_bound(order):
    value = derivatives.radial_partition(order)
    assert value.is_Rational and value > 0
    assert value <= 2048**order * s.factorial(order) ** 3
    assert derivatives.coefficient(order) == 2056**order * s.factorial(order) ** 3
    assert derivatives.amplitude_jet(order, source.VERROR) > 0
    assert 2 + s.sqrt(96) / 8 < 4


@pytest.mark.parametrize("order", range(196))
def test_every_exact_decreasing_high_phase_derivative_ratio(order):
    ratio = derivatives.growth(order + 1) / derivatives.growth(order)
    assert ratio == 2056 * (order + 1) ** 3 / source.RADIUS
    assert 0 < ratio < s.Rational(1, 10**9)


def test_all_high_derivatives_are_phase_not_time_profile_derivatives():
    assert source.MAX_ORDER == 196
    assert 2 + s.sqrt(96) / 4 > 4
    assert 2 + s.sqrt(96) / 8 < 4
    ratio = s.Rational(1, 32)
    assert 9 * ratio + 18 * ratio / (1 - ratio) == s.Rational(855, 992) < 1
    assert derivatives.amplitude_jet(4, source.VERROR) > derivatives.amplitude_jet(
        196, source.VERROR
    )
    assert 2**48 * 10**96 < kernel.CONSTANT


@pytest.mark.parametrize("a", [1.2, 2.0, 3.0])
def test_positive_Gaussian_symbol_is_not_generically_positive_Weyl_operator(a):
    # Exact Fock eigenvalues of OpW(exp(-a(q^2+p^2))).
    # The classical symbol is positive everywhere but the first excited
    # Weyl eigenvalue is negative if a>1.
    eigen0 = 1 / (1 + a)
    eigen1 = (1 - a) / (1 + a) ** 2
    integral, error = quad(
        lambda t: (2 * t - 1) * np.exp(-(1 + a) * t),
        0,
        np.inf,
        epsabs=1e-13,
        epsrel=1e-13,
    )
    assert eigen0 > 0 and eigen1 < 0
    assert_allclose(integral, eigen1, atol=5e-13)
    assert abs(eigen1) > 10 * error


@cache
def gaussian_derivative_ceiling(order):
    x = s.Symbol("normalized_coordinate", real=True)
    poly = s.expand(s.diff(s.exp(-x * x), x, order) / s.exp(-x * x))
    return sum(abs(c) * s.factorial(m[0]) for m, c in s.Poly(poly, x).terms())


@pytest.mark.parametrize("radius", [30, 80, 200])
def test_full_Weyl_and_calibrated_coherent_Gaussian_spectrum_against_Schur_error(
    radius,
):
    a = 1 / radius**2
    amplitude = 0.02
    # Coherent radial integral uses s=(q^2+p^2)/2 and the entire
    # exp(-s)s^n/n! measure, with the complete first calibration.
    for occupation in (0, 1, 2, 10, 100, 1000, 10000, 100000):
        Weyl = (1 - a) ** occupation / (1 + a) ** (occupation + 1)
        coherent = (1 + a) / (1 + 2 * a) ** (occupation + 1) - 2 * a * a * (
            occupation + 1
        ) / (1 + 2 * a) ** (occupation + 2)
        if occupation <= 2:

            def integrand(t, occupation=occupation):
                density = (
                    np.exp(occupation * np.log(t) - t - gammaln(occupation + 1))
                    if t
                    else float(occupation == 0)
                )
                return (1 + a - 2 * a * a * t) * np.exp(-2 * a * t) * density

            actual, error = quad(integrand, 0, np.inf, epsabs=2e-13, epsrel=2e-13)
            assert_allclose(actual, coherent, atol=10 * error + 3e-13)
        weights = [5, 4, 1]
        bound = 0.0
        for i in range(3):
            for j in range(3):
                full = (
                    gaussian_derivative_ceiling(i + 4) * gaussian_derivative_ceiling(j)
                    + 2
                    * gaussian_derivative_ceiling(i + 2)
                    * gaussian_derivative_ceiling(j + 2)
                    + gaussian_derivative_ceiling(i)
                    * gaussian_derivative_ceiling(j + 4)
                )
                bound += (
                    weights[i] * weights[j] * float(full) / radius ** (i + j + 4) / 32
                )
        bound *= np.pi * amplitude / 2
        assert amplitude * abs(coherent - Weyl) <= bound + 1e-13
        assert 1 + amplitude * Weyl > 0.99
        assert 1 + amplitude * coherent > 0.99
    assert (
        abs((1 + a) / (1 + 2 * a) - 2 * a * a / (1 + 2 * a) ** 2 - 1 / (1 + a))
        > a * a / 2
    )


@pytest.mark.parametrize("case", range(8))
def test_entire_two_ordering_UNITARY_operator_norm_and_changed_readout(case):
    # Independent noncommuting bounded-operator fixture. Not a claim that
    # these matrices are a finite CCR representation or the actual P8 H.
    rng = np.random.default_rng(271000 + case)
    matrices = []
    for _ in range(3):
        raw = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        matrices.append((raw + raw.conj().T) / 10)
    H0, H1, delta = matrices
    H0 += 0.43 * np.eye(4)
    delta *= 0.04
    T = 0.03
    assert np.linalg.norm(H0 @ H1 - H1 @ H0) > 0.001

    def evolve(extra):
        def rhs(t, flat):
            return (
                -1j * (H0 + t * H1 + np.cos(t) * extra) @ flat.reshape(4, 4)
            ).ravel()

        solution = solve_ivp(
            rhs,
            (0, T),
            np.eye(4, dtype=complex).ravel(),
            method="DOP853",
            rtol=2e-13,
            atol=2e-14,
        )
        assert solution.success
        return solution.y[:, -1].reshape(4, 4)

    UC, UW = evolve(np.zeros((4, 4))), evolve(delta)
    assert_allclose(UW.conj().T @ UW, np.eye(4), atol=2e-13)
    assert_allclose(UC.conj().T @ UC, np.eye(4), atol=2e-13)
    errorH = np.linalg.norm(delta, 2)
    assert np.linalg.norm(UW - UC, 2) <= T * errorH + 3e-13
    seed = rng.normal(size=4) + 1j * rng.normal(size=4)
    seed /= np.linalg.norm(seed)
    FC = np.diag([1.001, 1.003, 0.998, 1.002])
    change = np.diag([0.0001, -0.0002, 0.0003, -0.0001])
    FW = FC + change
    assert np.linalg.eigvalsh(FW).min() > 0.99
    errorF = np.linalg.norm(change, 2)
    psiC, psiW = UC @ seed, UW @ seed
    actual = abs(np.vdot(psiW, FW @ psiW) - np.vdot(psiC, FC @ psiC))
    assert actual <= errorF + 4 * T * errorH + 3e-13
    free = expm(-1j * 0.17 * H1)
    assert_allclose(
        np.vdot(free @ psiW, (free @ FW @ free.conj().T) @ (free @ psiW)),
        np.vdot(psiW, FW @ psiW),
        atol=2e-13,
    )


def test_all_actual_operator_bounds_are_positive_exact_and_separate_from_older_claims():
    b = comparison.bounds()
    assert all(value.is_Rational is True and value > 0 for value in b.values())
    assert b["whole_Weyl_volume_operator_ordering_error"] < s.Rational(1, 10**200)
    assert b["whole_Weyl_volume_positive_operator_floor"] > s.Rational(1, 2)
    assert b["whole_two_ordering_unitary_OPERATOR_norm_error"] < s.Rational(1, 10**944)
    assert b["whole_two_Weyl_cutoff_state_error"] > s.Rational(1, 10**1970)
    assert b["whole_two_Weyl_cutoff_readout_mean_error"] > s.Rational(1, 10**1230)
    assert b["whole_evolved_Weyl_coherent_outside_core_probability"] > 0


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[x[0] for x in audit.bad_cases()]
)
def test_unsupported_high_derivative_amplitude_operator_and_closure_inputs(
    name, call, args
):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_same_source_state_time_cutoffs_and_original_frontiers_retained():
    assert audit.require_parameters(audit.parameters()) == audit.parameters()
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 127
    for order in (0, 4, 192, 196):
        assert audit.require_order(order) == order
    for amplitude in (source.HBOUND, source.VERROR):
        assert audit.require_amplitude(amplitude) == amplitude
    for width in (1, 2):
        assert audit.require_cutoff(width) == width
    for time in (-source.TIME, 0, source.TIME):
        assert audit.require_time(time) == time
