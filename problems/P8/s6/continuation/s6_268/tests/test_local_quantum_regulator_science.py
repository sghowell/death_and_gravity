"""Independent full-covariance quadrature, phase tails and time-ordered dynamics."""

from functools import cache
from itertools import product
from math import factorial

import numpy as np
import pytest
import sympy as s
from numpy.testing import assert_allclose
from p8_vacuum_affine_local_quantum_regulator import (
    audit,
    dynamics,
    localization,
    ordering,
    source,
)
from scipy.integrate import solve_ivp
from scipy.special import gammainc, gammaincc, gammaln


@pytest.mark.parametrize(
    "packet", [source.data, ordering.data, localization.data, dynamics.data]
)
def test_whole_packet_exact_residuals_and_gates(packet):
    data = packet()
    for name, value in data["checks"].items():
        assert all(
            v == 0
            for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
        ), name
    assert all(value is True for value in data["gates"].values())


def hermite_grid(dimension, order):
    nodes, weights = np.polynomial.hermite.hermgauss(order)
    indices = np.array(list(product(range(order), repeat=dimension)))
    return nodes[indices], np.prod(weights[indices], axis=1) / np.pi ** (dimension / 2)


@pytest.mark.parametrize("case", range(8))
@pytest.mark.parametrize("amount", [1, 2])
def test_entire_correlated_Gaussian_heat_against_independent_quadrature(case, amount):
    variables, V, _, _ = ordering.correlated_fixture()
    rng = np.random.default_rng(26800 + case)
    coefficients = rng.integers(-3, 4, (3, 4))
    linear = [sum(int(row[i]) * variables[i] for i in range(4)) for row in coefficients]
    polynomial = s.expand(
        linear[0] ** 4 + linear[1] ** 3 + linear[0] * linear[1] * linear[2] + 7
    )
    covariance = amount * np.array(V.tolist(), dtype=float)
    grid, weights = hermite_grid(4, 4)
    samples = np.sqrt(2) * grid @ np.linalg.cholesky(covariance).T
    actual = np.dot(weights, s.lambdify(variables, polynomial, "numpy")(*samples.T))
    expected = float(
        ordering.polynomial_heat(polynomial, variables, V, amount).subs(
            dict.fromkeys(variables, 0)
        )
    )
    assert_allclose(actual, expected, rtol=3e-14, atol=2e-11)


@pytest.mark.parametrize("degree", range(8))
def test_whole_polynomial_calibration_has_exact_heat_integral_remainder(degree):
    x, V, _, _ = ordering.correlated_fixture()
    f = (x[0] + 2 * x[1] - x[2] + 3 * x[3]) ** degree
    D = lambda a: ordering.heat_generator(a, x, V)
    matched = ordering.polynomial_heat(f - D(f), x, V)
    t = s.Symbol("heat_time", real=True)
    integral = s.integrate(t * ordering.polynomial_heat(D(D(f)), x, V, t), (t, 0, 1))
    assert s.expand(matched - f + integral) == 0
    if degree <= 3:
        assert s.expand(matched - f) == 0


@pytest.mark.parametrize("dimension", [1, 2, 8, 16, 256])
@pytest.mark.parametrize(
    "radial", [s.Rational(1, 4), s.Integer(3), s.Integer(20), s.Integer(300)]
)
@pytest.mark.parametrize("occupation", [0, 1, 4])
def test_full_dimension_tail_against_incomplete_gamma(dimension, radial, occupation):
    t = s.Symbol("half_square_radius", real=True)
    degree = dimension + occupation
    tail = s.exp(-radial) * localization.tail_polynomial(degree, t).eval(radial)
    actual = float(s.N(tail, 60))
    expected = gammaincc(degree, float(radial))
    assert_allclose(actual, expected, atol=3e-14, rtol=3e-13)
    # Underflow of tiny positive CDFs is not treated as an exact zero eigenvalue.
    assert 0 <= gammainc(degree, float(radial)) <= 1
    next_tail = s.exp(-radial) * localization.tail_polynomial(degree + 1, t).eval(
        radial
    )
    gap = float(s.N(next_tail - tail, 60))
    stable_gap = np.exp(
        -float(radial) + degree * np.log(float(radial)) - gammaln(degree + 1)
    )
    assert_allclose(gap, stable_gap, atol=3e-14, rtol=5e-12)


@pytest.mark.parametrize("dimension", [1, 2, 8, 16, 256])
@pytest.mark.parametrize("multiple", [1.5, 2.0, 4.0])
def test_dimension_dependent_Chernoff_bound(dimension, multiple):
    t = multiple * dimension
    tail = gammaincc(dimension, t)
    bound = np.exp(-t + dimension + dimension * np.log(t / dimension))
    assert tail <= bound * (1 + 1e-13) + 1e-300


@cache
def numeric_fixture():
    fix = dynamics.fixture()
    return {
        name: np.array(fix[name].tolist(), dtype=complex)
        for name in ("H0", "H1", "F", "J", "psi")
    }


@pytest.mark.parametrize("which", ["H0", "H1", "F", "identity"])
def test_entire_coherent_matrix_integral_against_Gaussian_quadrature(which):
    # Scale the Gaussian into the quadrature weight; all remaining factors
    # are polynomials, so this independently integrates every matrix entry.
    beta = {"H0": 0.5, "H1": 2 / 3, "F": 0.5, "identity": 0.0}[which]
    grid, weights = hermite_grid(4, 7)
    grid /= np.sqrt(1 + beta)
    alpha = grid[:, :2] + 1j * grid[:, 2:]
    q, p = np.sqrt(2) * alpha.real, np.sqrt(2) * alpha.imag
    q2, p2 = np.sum(q * q, axis=1), np.sum(p * p, axis=1)
    qp = np.sum(q * p, axis=1)
    values = {
        "H0": q2 - p2 + (q2 + p2) ** 2 / 10,
        "H1": qp + qp**2 / 20,
        "F": (q2 - p2) ** 2 + 1,
        "identity": np.ones(len(grid)),
    }[which]
    basis = dynamics.fixture()["basis"]
    monomials = np.array(
        [
            alpha[:, 0] ** m[0]
            * alpha[:, 1] ** m[1]
            / np.sqrt(factorial(m[0]) * factorial(m[1]))
            for m in basis
        ]
    )
    actual = (monomials * (weights * values / (1 + beta) ** 2)) @ monomials.conj().T
    expected = np.eye(len(basis)) if which == "identity" else numeric_fixture()[which]
    assert_allclose(actual, expected, atol=2e-14, rtol=2e-13)


@cache
def numerical_Dyson_terms(last):
    fix = numeric_fixture()
    terms = [[np.eye(10, dtype=complex)]]
    for _ in range(last):
        previous = terms[-1]
        new = [np.zeros((10, 10), dtype=complex) for _ in range(len(previous) + 2)]
        for degree, coefficient in enumerate(previous):
            new[degree + 1] += -1j * (fix["H0"] @ coefficient) / (degree + 1)
            new[degree + 2] += -1j * (fix["H1"] @ coefficient) / (degree + 2)
        terms.append(new)
    return terms


def evaluate_Dyson(time, order):
    total = np.zeros((10, 10), dtype=complex)
    for coefficients in numerical_Dyson_terms(order):
        polynomial = np.zeros((10, 10), dtype=complex)
        for coefficient in reversed(coefficients):
            polynomial = polynomial * time + coefficient
        total += polynomial
    return total


@cache
def independent_propagator(time):
    fix = numeric_fixture()
    result = solve_ivp(
        lambda t, y: (-1j * (fix["H0"] + t * fix["H1"]) @ y.reshape(10, 10)).ravel(),
        (0, time),
        np.eye(10, dtype=complex).ravel(),
        method="DOP853",
        rtol=2e-13,
        atol=2e-15,
    )
    assert result.success
    return result.y[:, -1].reshape(10, 10)


@pytest.mark.parametrize("time", [0.1, 0.3, 0.6])
@pytest.mark.parametrize("order", [4, 8, 12])
def test_full_noncommuting_time_order_Dyson_bound_and_readout(time, order):
    fix = numeric_fixture()
    U = independent_propagator(time)
    approximate = evaluate_Dyson(time, order)
    norm_integral = (
        time * np.linalg.norm(fix["H0"], "fro")
        + time * time * np.linalg.norm(fix["H1"], "fro") / 2
    )
    bound = norm_integral ** (order + 1) / factorial(order + 1)
    actual = np.linalg.norm(U - approximate, 2)
    assert actual <= bound + 5e-13  # Independent floating ODE tolerance is explicit.
    assert_allclose(U.conj().T @ U, np.eye(10), atol=3e-13)
    assert_allclose(U @ fix["J"], fix["J"] @ U, atol=3e-13)
    state, approximate_state = U @ fix["psi"], approximate @ fix["psi"]
    assert_allclose(fix["J"] @ state, 0, atol=3e-13)
    exact_mean = (state.conj().T @ fix["F"] @ state)[0, 0]
    approximate_mean = (approximate_state.conj().T @ fix["F"] @ approximate_state)[0, 0]
    assert abs(exact_mean.imag) < 2e-13 and exact_mean.real > 0
    assert (
        abs(exact_mean - approximate_mean)
        <= np.linalg.norm(fix["F"], 2) * (2 * bound + bound**2) + 3e-12
    )


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[x[0] for x in audit.bad_cases()]
)
def test_unsupported_domain_state_ordering_and_closure_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_all_original_parameters_and_frontier_retained():
    assert audit.require_parameters(audit.parameters()) == audit.parameters()
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 124
    assert audit.require_radius(s.Rational(3, 2)) == s.Rational(3, 2)
    assert audit.require_dimension(16) == 16
