"""Independent complete Hamiltonian, nonlinear-root and response diagnostics."""

import numpy as np
import pytest
import sympy as s
from numpy.testing import assert_allclose
from p8_vacuum_affine_finite_window_growth.intervals import I, evaluate
from p8_vacuum_affine_nonlinear_lapse_branch import branch, response
from p8_vacuum_affine_nonlinear_lapse_branch import diagnostics as num
from p8_vacuum_affine_nonlinear_lapse_branch import source as q


@pytest.mark.parametrize(
    "packet",
    [q.data, branch.enclosure, branch.auxiliaries, response.jets, response.physical],
)
def test_whole_packet_exact_residuals_and_bounds(packet):
    data = packet()
    assert all(bool(value) for value in data["gates"].values())
    for name, value in data["checks"].items():
        assert all(
            x == 0
            for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
        ), name


@pytest.mark.parametrize("case", range(5))
@pytest.mark.parametrize("scale", [1e-6, 1e-5, 3e-5, 7e-5])
def test_full_nonlinear_lapse_and_temporal_roots(case, scale):
    values = num.vector(case, scale)
    N = num.root(case, values)
    data = num.fixture(case)
    assert 0.999 < N < 1.001
    assert abs(data["C"](N, *values)) < 5e-13
    T = float(data["T"](N, *values))
    R = data["source"](N)[0]
    U, M, r = R ** (-0.75), R**0.25, R - 1
    p, G = values[:2]
    aa = -M / 3
    K = (p + U * r * T) / (2 * aa + U * r * r)
    assert abs(U * (T - r * K) + G) < 1e-14
    assert 1 - 3 * r * r / (2 * R) > 0.99


@pytest.mark.parametrize("case", range(5))
@pytest.mark.parametrize("component", range(12))
def test_every_full_invariant_response_against_independent_root_solves(case, component):
    values = num.vector(case)
    data = num.fixture(case)
    N = num.root(case, values)
    step = 2e-6
    left = values.copy()
    left[component] -= step
    right = values.copy()
    right[component] += step
    Cdiff = (data["C"](N, *right) - data["C"](N, *left)) / (2 * step)
    predicted = -Cdiff / data["C_N"](N, *values)
    observed = (num.root(case, right) - num.root(case, left)) / (2 * step)
    assert_allclose(observed, predicted, atol=2e-8, rtol=2e-6)


@pytest.mark.parametrize("case", range(5))
def test_literal_full_Hamiltonian_and_nonzero_primitive_against_numeric_variation(case):
    values = num.vector(case, 3e-5)
    data = num.fixture(case)
    for N in [0.9998, num.root(case, values), 1.0002]:
        expected = N * data["normal"](N, *values) + num.primitive(case, N)
        assert_allclose(
            num.independent_hamiltonian(case, N, values),
            expected,
            atol=2e-14,
            rtol=2e-14,
        )
        step = 1e-4
        h = lambda x: num.independent_hamiltonian(case, x, values)
        derivative = (
            -h(N + 2 * step) + 8 * h(N + step) - 8 * h(N - step) + h(N - 2 * step)
        ) / (12 * step)
        assert_allclose(derivative, data["C"](N, *values), atol=2e-10, rtol=2e-8)
    assert abs(num.primitive(case, 1.0002)) > 1e-9


@pytest.mark.parametrize("case", range(5))
def test_whole_contraction_against_independent_bracket_root(case):
    values = num.vector(case)
    data = num.fixture(case)
    target = num.root(case, values)
    N = 1.0
    errors = []
    for _ in range(6):
        errors.append(abs(N - target))
        N += data["C"](N, *values) / 3
    assert errors[-1] < 2e-13
    assert errors[1] < errors[0] / 20


@pytest.mark.parametrize("case", range(5))
def test_complete_auxiliary_Schur_matrix_and_inverse(case):
    values = num.vector(case)
    N = num.root(case, values)
    data = num.fixture(case)
    R = data["source"](N)[0]
    gamma = 1 - 3 * (R - 1) ** 2 / (2 * R)
    temporal = -N * R ** (-0.75) / gamma
    slope = data["T_N"](N, *values)
    cn = data["C_N"](N, *values)
    matrix = np.array(
        [
            [cn + temporal * slope * slope, -temporal * slope],
            [-temporal * slope, temporal],
        ],
        float,
    )
    inverse = np.array(
        [[1 / cn, slope / cn], [slope / cn, 1 / temporal + slope * slope / cn]], float
    )
    assert np.linalg.eigvalsh(matrix)[-1] < -0.5
    assert matrix[0, 0] * matrix[1, 1] - matrix[0, 1] ** 2 > 1.5
    assert_allclose(matrix @ inverse, np.eye(2), atol=2e-15)
    # Independent derivative of the complete reduced Hamiltonian.
    step = 1e-4
    h = lambda x: num.independent_hamiltonian(case, x, values)
    second = (
        -h(N + 2 * step)
        + 16 * h(N + step)
        - 30 * h(N)
        + 16 * h(N - step)
        - h(N - 2 * step)
    ) / (12 * step**2)
    assert_allclose(second, cn, atol=1e-6, rtol=2e-6)


@pytest.mark.parametrize("dimension", [2, 3, 4, 5])
def test_full_nonlinear_chart_chain_keeps_every_second_contact(dimension):
    variables = s.symbols("q0:" + str(dimension), real=True)
    invariant = s.Matrix(
        [
            variables[i]
            + variables[(i + 1) % dimension] ** 2
            + variables[i] * variables[(i + 2) % dimension]
            for i in range(dimension)
        ]
    )
    z = s.symbols("z0:" + str(dimension), real=True)
    function = sum(
        (i + 1) * z[i] ** 3 + z[i] * z[(i + 1) % dimension] for i in range(dimension)
    )
    mapping = dict(zip(z, invariant, strict=True))
    gradient = s.Matrix(
        [s.diff(function, x).subs(mapping, simultaneous=True) for x in z]
    )
    hessian = s.hessian(function, z).subs(mapping, simultaneous=True)
    mapped1, mapped2 = response.pullback(
        gradient,
        hessian,
        invariant.jacobian(variables),
        [s.hessian(x, variables) for x in invariant],
    )
    literal = function.subs(mapping, simultaneous=True)
    assert all(
        s.expand(x) == 0
        for x in mapped1 - s.Matrix([s.diff(literal, x) for x in variables])
    )
    assert all(s.expand(x) == 0 for x in mapped2 - s.hessian(literal, variables))
    incomplete = (
        invariant.jacobian(variables).T * hessian * invariant.jacobian(variables)
    )
    assert any(s.expand(x) != 0 for x in mapped2 - incomplete)


@pytest.mark.parametrize("case", range(8))
def test_all_third_implicit_contacts_against_exact_implicit_curve(case):
    # Construct C(x,t) with a known nonlinear root, not by using the tested formula.
    t, x = s.symbols("t x", real=True)
    exact = (case + 1) * t + (case + 2) * t * t / 3 - (case + 3) * t**3 / 5
    C = (x - exact) * (2 + (case + 1) * t + 3 * x + t * x + x * x)
    at = {t: 0, x: 0}
    derivative = lambda i, j: s.diff(C, x, i, t, j).subs(at)
    result = response.third_direction(
        derivative(1, 0),
        derivative(0, 1),
        derivative(1, 1),
        derivative(0, 2),
        derivative(2, 0),
        derivative(3, 0),
        derivative(2, 1),
        derivative(1, 2),
        derivative(0, 3),
    )
    assert all(
        s.factor(value - s.diff(exact, t, k).subs(t, 0)) == 0
        for k, value in enumerate(result, 1)
    )


def test_actual_outward_roots_are_centered_not_independent_subtractions():
    packet = branch.enclosure()
    assert packet["whole_centered_MVT_constraint_residual_bound"] / 3 < q.INNER
    assert packet["whole_source_only_center_residual_bound"] / 3 < q.CENTER_RADIUS
    assert q.CENTER_RADIUS < q.INNER < q.fifth.EPSILON
    # Explicitly demonstrate the inherited fixed-width root interval floor.
    tiny = s.Rational(1, 10**380)
    iv = evaluate(
        s.Symbol("positive_probe") ** s.Rational(1, 4),
        {s.Symbol("positive_probe"): I(1 - tiny, 1 + tiny)},
    )
    assert iv.hi - iv.lo > tiny


def test_full_constraint_not_the_naive_linear_Gaussian_lapse():
    assert s.diff(q.CONSTRAINT, q.ph, 2) != 0
    assert s.diff(q.CONSTRAINT, q.eta, 2) != 0
    assert s.diff(q.CONSTRAINT, q.p, q.G) != 0
    assert s.diff(q.CONSTRAINT, q.G, 2) != 0
    assert s.diff(q.CONSTRAINT, q.N, q.curv) != 0


@pytest.mark.parametrize("case", range(5))
def test_diagnostic_sources_keep_both_switches_and_finite_heavy_force(case):
    data = num.fixture(case)
    values = num.vector(case, 5e-5)
    first = values.copy()
    first[4] += 3e-4
    second = values.copy()
    second[4] -= 3e-4
    assert abs(data["C"](1, *first) - data["C"](1, *second)) > 1e-9
    assert float(data["source"](0.999)[0]) > 0
    assert float(data["source"](1.001)[0]) > 0


from p8_vacuum_affine_nonlinear_lapse_branch import audit


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_invalid_inputs_and_false_scope_promotions_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("component", range(12))
@pytest.mark.parametrize("sign", [-1, 1])
def test_exact_closed_invariant_box_faces_allowed(component, sign):
    values = {str(z): 0 for z in q.COORDS}
    values[str(q.COORDS[component])] = sign * q.DELTA
    assert audit.require_coordinates(values)[str(q.COORDS[component])] == sign * q.DELTA


@pytest.mark.parametrize("observable", audit.OBSERVABLES)
def test_only_whole_nonlinear_and_physical_chain_observables_allowed(observable):
    assert audit.require_observable(observable) == observable


def test_original_parameters_state_measure_and_frontier_preserved():
    assert audit.require_parameters(audit.parameters()) == audit.parameters()
    assert audit.require_state(audit.STATE) == audit.STATE
    assert audit.require_measure(audit.MEASURE) == audit.MEASURE
    assert audit.require_time(0) == 0
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 122
    assert audit.rejected_inputs() == 291
