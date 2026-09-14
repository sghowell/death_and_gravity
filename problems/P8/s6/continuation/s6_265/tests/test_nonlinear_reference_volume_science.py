"""Independent full-function, spectral, covariance and integrability diagnostics."""

import math
from itertools import pairwise

import numpy as np
import pytest
import sympy as s
from p8_affine import verify as certificate
from p8_vacuum_affine_coupled_gaussian_state import gaussian as original_gaussian
from p8_vacuum_affine_coupled_gaussian_state import phase
from p8_vacuum_affine_nonlinear_reference_volume import (
    audit,
    gaussian,
    local,
    physical,
    source,
)
from scipy import integrate, linalg, special


@pytest.mark.parametrize(
    "name",
    (
        "whole_original_R_and_global_positive_lapse_coefficient_bounds",
        "whole_full_original_Hamiltonian_and_lapse_weight_bindings",
        "whole_complete_finite_strip_inverse_power_derivative_bounds",
        "whole_strongly_commuting_bounce_rows_and_nonzero_covariance",
        "whole_complete_Gaussian_fourth_remainders_and_tail_bounds",
        "whole_positive_lapse_nonlinear_reference_composites",
        "whole_singular_individual_coefficient_integrability_boundary",
    ),
)
def test_all_complete_exact_packets(name):
    data = audit.packets()[name]
    assert certificate.certify_residuals(data["checks"])
    assert all(bool(value) for value in data["gates"].values())


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=lambda x: x if isinstance(x, str) else None
)
def test_every_unsupported_input_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("alpha", source.ALPHAS)
def test_only_full_stated_negative_power_functions(alpha):
    assert audit.require_power(alpha) == alpha


def test_bounce_spectral_time_is_exactly_zero_not_tiny_interval():
    assert audit.require_joint_time(0) == 0
    with pytest.raises(ValueError):
        audit.require_joint_time(source.TIME / 100)


def test_unchanged_frontier_and_explicit_no_go_boundary():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 121
    assert audit.matching()[:-1] == audit.previous.matching()
    assert all(audit.gates().values())
    assert "OPEN" in audit.observable()["original_problem"]
    assert "not" in audit.observable()["strict_probability_boundary"].lower()
    assert "infinite" in audit.observable()["negative_result"].lower()


def complete_fixture_R(u, N, order=8, localizer=17):
    # Deliberately finite parameter fixtures of the complete TWO-switch
    # structure, not the actual m=1024,A=1e420 parent or a new P8 state.
    if N <= 0:
        raise ValueError("Positive lapse only")
    X = N**-2
    T = 1.0 if X == 1 else special.expit(order * (math.log(X) - math.log(abs(1 - X))))
    bump = order * X * X * math.exp(-order * X * X)
    B = T + (1 - T) * bump
    heavy = order * X * X * (1 - X) ** 8 * math.exp(-localizer * X * X)
    return 1 + B * (X - 1) / (1 + u * u) ** 3 + heavy


@pytest.mark.parametrize("order,localizer", [(8, 17), (16, 29), (32, 67)])
@pytest.mark.parametrize("u", [-0.3, 0, 0.4])
def test_full_positive_source_and_large_X_coefficients_with_nonzero_switch_fixtures(
    order, localizer, u
):
    h = (1 + u * u) ** 3
    for N in np.geomspace(0.005, 20, 91):
        R = complete_fixture_R(u, N, order, localizer)
        assert R > 0.5
        if N <= 1:
            assert 1 / (2 * h * N * N) <= R <= 2 / (N * N)
        else:
            assert R <= 2
        for alpha in (0.25, 0.5, 0.75):
            assert 0 < R**-alpha < 2
    # This is an actual nonzero heavy term in the diagnostic, not a
    # floating-point underflow counted as exact source elimination.
    X = 0.2
    delta = order * X * X * (1 - X) ** 8 * math.exp(-localizer * X * X)
    assert delta > 1e-5
    far = complete_fixture_R(u, 1e-5, order, localizer) * 1e-10
    assert far == pytest.approx(1 / (2 * h), rel=2e-8)


@pytest.mark.parametrize("case", range(12))
def test_pure_full_covariance_inverse_and_lower_operator_bound_independently(case):
    _, _, M = original_gaussian.fixture_preparation(case)
    V = np.array(original_gaussian.ground_covariance(M), float)
    omega = np.array(phase.OMEGA, float)
    np.testing.assert_allclose(
        linalg.inv(V, check_finite=True), -4 * omega @ V @ omega, rtol=2e-8, atol=2e-8
    )
    values = linalg.eigvalsh(V, check_finite=True)
    assert values[0] >= 1 / (4 * values[-1]) * (1 - 1e-8)
    P = 3.0 + case
    J = 1.52
    T = 0.01
    rows = np.array(
        [
            [0, 0, 1 / (2 * P * P), 0],
            [0, 0, 1 / (4 * J) - 3 * T / (4 * J * P * P), -1 / (40 * J)],
        ]
    )
    covariance = rows @ (P * V) @ rows.T
    assert linalg.eigvalsh(covariance, check_finite=True)[0] > 0
    assert covariance[1, 1] >= P * rows[1, 3] ** 2 / (4 * values[-1]) * (1 - 1e-9)
    np.testing.assert_allclose(
        rows @ (P * omega) @ rows.T, np.zeros((2, 2)), atol=1e-14
    )
    # Purity is necessary for the exact inverse identity, not silently
    # promoted to every mixed positive covariance.
    mixed = V + np.eye(4)
    assert (
        np.linalg.norm(linalg.inv(mixed, check_finite=True) + 4 * omega @ mixed @ omega)
        > 1
    )


@pytest.mark.parametrize("alpha", source.ALPHAS)
def test_complete_lapse_center_jets_with_independent_finite_symbolic_composition(alpha):
    N = s.Symbol("independent_positive_N", positive=True)
    h = s.Symbol("independent_positive_h", positive=True)
    q = 1 + (N**-2 - 1) / h
    exact = [s.factor(s.diff(q**-alpha, N, j).subs(N, 1)) for j in range(5)]
    for u in (-s.Rational(1, 10), 0, s.Rational(1, 10)):
        expected = [value.subs(h, (1 + u * u) ** 3) for value in exact]
        actual = [
            value.subs(source.u, u)
            for value in local.derivatives()[
                "whole_complete_center_jets_zero_through_four"
            ][str(alpha)]
        ]
        assert all(s.factor(a - b) == 0 for a, b in zip(actual, expected))


def joint_volume_mean(scale, correlation=0.3):
    vv = 0.36 * scale * scale
    nn = scale * scale
    vn = correlation * scale * scale
    mean = 3 * vn
    sd = math.sqrt(nn)
    lower = (-1 - mean) / sd

    def integrand(z):
        N = 1 + mean + sd * z
        if N <= 0:
            return 0.0
        return (
            math.exp(-z * z / 2)
            / math.sqrt(2 * math.pi)
            * complete_fixture_R(0, N) ** -0.75
        )

    result = integrate.quad(
        integrand, max(lower, -13), 13, epsabs=2e-13, epsrel=2e-13, limit=200
    )[0]
    return math.exp(4.5 * vv) * result, 1 + 4.5 * vv + 4.5 * vn + 0.375 * nn


@pytest.mark.parametrize("correlation", [-0.3, 0.0, 0.3])
def test_entire_commuting_volume_Gaussian_integral_has_fourth_order_remainder(
    correlation,
):
    errors = []
    for scale in (0.08, 0.04, 0.02, 0.01):
        full, jet = joint_volume_mean(scale, correlation)
        error = abs(full - jet)
        errors.append(error)
        assert full > 0 and error < 100 * scale**4
    assert all(a / b > 6 for a, b in pairwise(errors))
    if correlation:
        actual, _ = joint_volume_mean(0.02, correlation)
        erased, _ = joint_volume_mean(0.02, 0)
        assert abs(actual - erased) > 0.0004


@pytest.mark.parametrize("alpha", [0.25, 0.5, 0.75])
@pytest.mark.parametrize("u", [-0.1, 0, 0.1])
def test_single_full_nonlinear_coefficient_Gaussian_integral_away_from_bounce(alpha, u):
    h = (1 + u * u) ** 3
    second = 4 * alpha * (alpha + 1) / (h * h) - 6 * alpha / h
    for sigma in (0.04, 0.02):

        def integrand(z, sigma=sigma):
            N = 1 + sigma * z
            return (
                math.exp(-z * z / 2)
                / math.sqrt(2 * math.pi)
                * complete_fixture_R(u, N) ** -alpha
                if N > 0
                else 0.0
            )

        value = integrate.quad(
            integrand, -13, 13, epsabs=2e-13, epsrel=2e-13, limit=200
        )[0]
        assert 0 < value < 2
        assert abs(value - (1 + second * sigma * sigma / 2)) < 100 * sigma**4


@pytest.mark.parametrize("case", range(4))
def test_full_cubic_Gaussian_polynomial_by_correlated_Hermite_cubature(case):
    covariance = np.array(
        [[0.04 + case / 100, (-1) ** case * 0.013], [(-1) ** case * 0.013, 0.06]]
    )
    L = linalg.cholesky(covariance, lower=True, check_finite=True)
    nodes, weights = np.polynomial.hermite.hermgauss(7)
    polynomial = gaussian.composites()["whole_bounce_complete_cubic_volume_polynomial"]
    function = s.lambdify((gaussian.v, gaussian.n), polynomial, "numpy")
    value = sum(
        weights[i]
        * weights[j]
        * function(*(math.sqrt(2) * L @ np.array([x, y])))
        / math.pi
        for i, x in enumerate(nodes)
        for j, y in enumerate(nodes)
    )
    wanted = (
        1 + 4.5 * covariance[0, 0] + 4.5 * covariance[0, 1] + 0.375 * covariance[1, 1]
    )
    assert value == pytest.approx(wanted, rel=2e-13)
    assert abs(value - (1 + 4.5 * covariance[0, 0] + 0.375 * covariance[1, 1])) > 0.01


@pytest.mark.parametrize("sigma", [0.4, 0.5, 0.7])
def test_full_singular_kinetic_and_canonical_matter_square_integrals(sigma):
    density = lambda N: (
        math.exp(-((N - 1) ** 2) / (2 * sigma * sigma))
        / (sigma * math.sqrt(2 * math.pi))
    )
    target = 2 * (0.5) ** 0.25 * density(0)
    kinetic = []
    matter2 = []
    for cutoff in (1e-3, 1e-5, 1e-7):

        def kd(y):
            N = y**-2
            return 2 * complete_fixture_R(0, N) ** 0.25 / N * density(N) / y**3

        value = integrate.quad(
            kd, 1, cutoff**-0.5, epsabs=1e-10, epsrel=2e-11, limit=200
        )[0]
        kinetic.append(value)

        def q2(t):
            N = math.exp(-t)
            Q = N * complete_fixture_R(0, N) ** 0.75
            return Q * Q * density(N) * N

        matter2.append(
            integrate.quad(
                q2, 0, -math.log(cutoff), epsabs=1e-11, epsrel=2e-11, limit=200
            )[0]
        )
    assert all(b > a for a, b in pairwise(kinetic))
    assert abs(math.sqrt(1e-7) * kinetic[-1] - target) < 0.003
    expected = 0.5**1.5 * density(0) * math.log(100)
    assert abs(matter2[-1] - matter2[-2] - expected) < 1e-4
    assert density(0) > 0


@pytest.mark.parametrize("sigma", [0.4, 0.5, 0.7])
def test_complete_canonical_matter_first_moment_remains_finite(sigma):
    density = lambda N: (
        math.exp(-((N - 1) ** 2) / (2 * sigma * sigma))
        / (sigma * math.sqrt(2 * math.pi))
    )

    def near(y):
        if y == 0:
            return 2 * (0.5) ** 0.75 * density(0)
        N = y * y
        return N * complete_fixture_R(0, N) ** 0.75 * density(N) * 2 * y

    value = integrate.quad(near, 0, 1, epsabs=1e-11, epsrel=2e-11, limit=200)[0]
    value += integrate.quad(
        lambda N: N * complete_fixture_R(0, N) ** 0.75 * density(N),
        1,
        12,
        epsabs=1e-11,
        epsrel=2e-11,
        limit=200,
    )[0]
    assert math.isfinite(value) and value > 0
    bound = 4 / (sigma * math.sqrt(2 * math.pi)) + 2 * (1 + sigma)
    assert value < bound


def test_actual_tiny_tail_is_not_computed_as_floating_zero():
    data = physical.integrability()
    sigma = next(
        z
        for z in data["whole_positive_gaussian_density_at_zero_lapse"].free_symbols
        if str(z) == "positive_lapse_standard_deviation"
    )
    actual = data["whole_positive_gaussian_density_at_zero_lapse"].subs(
        sigma, s.Rational(1, 10**250)
    )
    assert actual.is_positive is True and actual != 0
    assert (
        physical.bounce()["whole_current_uniform_single_lapse_variance_lower_enclosure"]
        > source.VAR_LOWER
    )
    assert (
        gaussian.remainder()["whole_fourth_order_bulk_remainder_bound"]
        + gaussian.remainder()["whole_positive_series_tail_bound"]
        < source.VOLUME_ERROR
    )
