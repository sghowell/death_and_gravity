"""Whole exact identities and independent finite connection/measure diagnostics."""

import cmath
import math

import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_complement_measure import audit, gaussian, jets, split
from scipy.linalg import det, eigvalsh, solve


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
    ) == (74, 29960, 42, 9, 326, 115)
    assert all(value is True for value in audit.gates().values())
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()
    assert audit.require_parameters(audit.parameters()) == audit.parameters()


@pytest.mark.parametrize("p", (s.Rational(3, 8), s.Rational(1, 2), s.Rational(27, 50)))
def test_independent_full_56_and_60_determinants_and_inverse(p):
    data = split.matrices()
    A = np.asarray(data["A"].subs(split.P, p), dtype=float)
    K = np.asarray(data["kernel"], dtype=float)
    L = np.asarray(data["lift"].subs(split.P, p), dtype=float)
    full = np.asarray(data["parent"]["M_new"].subs(split.P, p), dtype=float)
    inverse = solve(full, np.eye(60), assume_a="sym", check_finite=True)
    independent = (
        K @ solve(A, K.T, assume_a="sym", check_finite=True)
        + L @ np.diag([1.0, -1.0, -1.0, -1.0]) @ L.T
    )
    assert np.linalg.norm(inverse - independent) < 1e-10
    expected = float(
        split.complement()["whole_56_direction_complement_determinant"].subs(split.P, p)
    )
    assert abs(det(A, check_finite=True) / expected - 1) < 1e-10
    assert audit.require_p(p) == p
    eigenvalues = eigvalsh(A, check_finite=True)
    assert (sum(eigenvalues > 0), sum(eigenvalues < 0)) == (26, 30)


@pytest.mark.parametrize("block_index", range(9))
def test_independent_exact_reference_inertia_congruence(block_index):
    data = split.matrices()
    block = data["block_matrices"][block_index].subs(split.P, s.Rational(1, 2))
    transform, diagonal, pivots, inertia = split.rational_congruence(block)
    assert transform.T * block * transform == diagonal
    assert abs(transform.det()) == 1
    assert s.prod(p.det() for p in pivots) == block.det()
    eigenvalues = eigvalsh(np.asarray(block, dtype=float), check_finite=True)
    assert inertia == (sum(eigenvalues > 0), sum(eigenvalues < 0))


@pytest.mark.parametrize("p", (s.Rational(3, 8), s.Rational(1, 2), s.Rational(27, 50)))
def test_independent_all_64_point_inverse_and_projective_measure(p):
    data = split.matrices()
    J = np.asarray(data["full"].subs(split.P, p), dtype=float)
    inverse = solve(J, np.eye(64), check_finite=True)
    expected = np.asarray(data["full_inverse"].subs(split.P, p), dtype=float)
    assert np.linalg.norm(inverse - expected) < 1e-11
    assert abs(det(J, check_finite=True) + 256) < 1e-9
    F, G = (
        np.asarray(data["gauge"], dtype=float),
        np.asarray(data["projective"], dtype=float),
    )
    assert np.array_equal(F @ G, 4 * np.eye(4))
    assert det(F @ G, check_finite=True) == 256
    assert (
        abs(det(J, check_finite=True) * det(inverse.T, check_finite=True) - 1) < 1e-12
    )


@pytest.mark.parametrize("p", (s.Rational(3, 8), s.Rational(1, 2), s.Rational(27, 50)))
@pytest.mark.parametrize("scale", (0.7, 3.0))
def test_full_56_convergent_Gaussian_limit_and_continued_phase(p, scale):
    data = split.matrices()
    A = np.asarray(data["A"].subs(split.P, p), dtype=float)
    eigenvalues = eigvalsh(A, check_finite=True)
    epsilon = 1e-9
    log_integral = (
        28 * math.log(2 * math.pi)
        - sum(cmath.log(epsilon - 1j * scale * item) for item in eigenvalues) / 2
    )
    expected_log_amplitude = (
        28 * math.log(2 * math.pi / scale) - math.log(det(A, check_finite=True)) / 2
    )
    ratio = cmath.exp(log_integral - expected_log_amplitude)
    assert abs(ratio + 1) < 2e-6
    assert abs(ratio - 1) > 1.9
    normalized = ratio * (-1)
    assert abs(normalized - 1) < 2e-6
    # The final determinant loses the eigenvalue-by-eigenvalue continuation phase.
    endpoint = (-1j) ** 56 * det(scale * A, check_finite=True)
    assert endpoint.real > 0 and abs(endpoint.imag) < 1e-15


@pytest.mark.parametrize("p", (s.Rational(3, 8), s.Rational(1, 2), s.Rational(27, 50)))
def test_full_56_source_contact_from_convergent_Gaussian(p):
    data = split.matrices()
    A = np.asarray(data["A"].subs(split.P, p), dtype=float)
    source = np.sin(np.arange(56) + 0.3) / 50
    scale = 1.7
    B = 1e-10 * np.eye(56) - 1j * scale * A
    regular_inverse = solve(B, np.eye(56), assume_a="sym", check_finite=True)
    exact_inverse = np.asarray(data["inverse_A"].subs(split.P, p), dtype=float) / scale
    assert np.linalg.norm(regular_inverse - 1j * exact_inverse) < 2e-6
    gaussian_source = cmath.exp(-0.5 * source @ regular_inverse @ source)
    expected = cmath.exp(-0.5j * source @ exact_inverse @ source)
    assert abs(gaussian_source - expected) < 1e-8
    assert np.linalg.norm(1j * regular_inverse @ source + exact_inverse @ source) < 1e-7


@pytest.mark.parametrize("p", (s.Rational(3, 8), s.Rational(1, 2)))
def test_independent_complete_64_source_centered_action(p):
    data = split.matrices()
    packet = split.source_map()
    full = s.Matrix(packet["whole_source_centered_64_updated_Hessian"])
    forcing = s.Matrix(packet["whole_source_centered_64_updated_linear_term"])
    center = s.Matrix(packet["whole_actual_64_connection_center_at_retained_trace"])
    symbols = sorted(
        full.free_symbols | forcing.free_symbols | center.free_symbols, key=str
    )
    substitutions = {symbol: s.Rational(i + 1, 29) for i, symbol in enumerate(symbols)}
    substitutions[split.P] = p
    substitutions[split.old.S] = s.Rational(6, 5)
    full, forcing, center = (
        value.subs(substitutions) for value in (full, forcing, center)
    )
    E = data["embedding"] * data["kernel"]
    z = s.Matrix([s.Rational((i % 7) - 3, 31) for i in range(56)])
    displaced = center + E * z
    original_change = (displaced.T * full * displaced - center.T * full * center)[
        0
    ] / 2 + forcing.dot(displaced - center)
    reduced = (z.T * data["A"].subs(split.P, p) * z)[0] / 2
    assert original_change == reduced
    assert data["gauge"] * displaced == s.zeros(4, 1)


@pytest.mark.parametrize("seed", (7, 13, 29))
def test_projective_source_and_full_contact_without_retained_Proca(seed):
    data = split.matrices()
    raw = s.Matrix([s.Rational((i * seed) % 11 - 5, 17) for i in range(64)])
    G = data["projective"]
    source = raw - G * (G.T * raw) / 4
    assert audit.require_projective_source(source) == source
    E = data["embedding"] * data["kernel"]
    covariance = gaussian.source_contacts()[
        "whole_64_connection_source_contact_kernel"
    ].subs(split.P, s.Rational(1, 2))
    assert data["gauge"] * covariance == s.zeros(4, 64)
    result = covariance * source
    quotient = data["full_inverse"].subs(split.P, s.Rational(1, 2)) * result
    assert quotient[56:, :] == s.zeros(8, 1)
    assert (
        s.simplify(
            source.dot(result)
            - (E.T * source).dot(
                data["inverse_A"].subs(split.P, s.Rational(1, 2)) * (E.T * source)
            )
        )
        == 0
    )


@pytest.mark.parametrize("step", (1e-4, 3e-4))
def test_full_56_logdet_first_second_R_variations(step):
    data = split.matrices()

    def logdet(R):
        p = math.sqrt(R) / 2
        A = np.asarray(data["A"].subs(split.P, p), dtype=float)
        return math.log(det(A, check_finite=True))

    base, plus, minus = (logdet(R) for R in (1, 1 + step, 1 - step))
    assert abs((plus - minus) / (2 * step) - 116 / 3) < 1e-5
    assert abs((plus + minus - 2 * base) / step**2 + 1228 / 27) < 2e-4


@pytest.mark.parametrize("value", (1, 2, 13, 10**20))
def test_finite_cell_and_positive_scale_scope(value):
    assert audit.require_cells(value) == value
    assert audit.require_cell_scale(s.Rational(value, 7)) == s.Rational(value, 7)
    assert audit.require_measure(audit.MEASURE) == audit.MEASURE


@pytest.mark.parametrize("q,v", ((0, 0), (s.Rational(1, 3), s.Rational(2, 5)), (-1, 2)))
def test_independent_jet_constraint_density_and_physical_Hamiltonian(q, v):
    data = jets.primary_secondary()
    variables = data["whole_finite_jet_base_auxiliary_variables"]
    sub = {variables[0]: q, variables[1]: v}
    A = data["whole_finite_indefinite_auxiliary_Hessian"].subs(sub)
    bracket = data["whole_six_constraint_bracket_on_auxiliary_surface"].subs(sub)
    assert bracket.det() == A.det() ** 2
    assert (
        s.factor(
            data["whole_reduced_physical_Hamiltonian"]
            - (variables[4] ** 2 + variables[0] ** 2) / 2
        )
        == 0
    )
    assert A.det() < 0


def test_jet_lift_does_not_assume_original_velocity_Legendre_regular():
    data = jets.primary_secondary()
    q, v = data["whole_finite_jet_base_auxiliary_variables"][:2]
    singular = data["whole_nonempty_unextended_Legendre_singularity"]
    original = data["whole_original_unextended_velocity_Hessian_at_auxiliary_solution"]
    assert s.simplify(original.subs({q: 0, v: s.sqrt(singular)})) == 0
    assert (
        data["whole_finite_indefinite_auxiliary_Hessian"]
        .subs({q: 0, v: s.sqrt(singular)})
        .det()
        != 0
    )


def test_entire_source_pullback_changes_remaining_velocity_constraint():
    packet = jets.source_dependent_jet()
    assert all(
        s.factor(entry) == 0
        for value in packet["checks"].values()
        for entry in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    reduced = packet["whole_source_dependent_retained_jet_Lagrangian"]
    source2 = next(
        symbol
        for symbol in reduced.free_symbols
        if str(symbol) == "original_connection_source2"
    )
    velocity = next(
        symbol
        for symbol in reduced.free_symbols
        if str(symbol) == "independent_velocity"
    )
    zeros = {
        symbol: 0
        for symbol in reduced.free_symbols
        if str(symbol).startswith("original_connection_source")
    }
    assert s.factor(s.diff(reduced, source2, velocity, 2).subs(zeros)) == 2


def test_independent_nonlinear_point_ordering_difference():
    data = jets.ordering_boundary()
    whole = data["whole_unitarily_pulled_back_free_Hamiltonian"]
    naive = data["incorrect_unchanged_naive_Weyl_symbol_quantization"]
    difference = s.simplify(whole - naive)
    assert difference == data["whole_nonzero_ordering_contact"] and difference != 0


def test_same_parent_and_explicit_nonclosure():
    assert "OPEN" in audit.observable()["original_problem"]
    assert "fixed quantum mean" in audit.observable()["scope_boundary"]
    assert audit.require_state(audit.STATE) == audit.STATE


@pytest.mark.parametrize("parameter", (0.2, 0.7, -0.3))
def test_independent_full_56_basis_density_and_source_covariance(parameter):
    data = split.matrices()
    A = np.asarray(data["A"].subs(split.P, s.Rational(1, 2)), dtype=float)
    B = np.eye(56)
    B[:2, :2] = np.array([[1 + parameter**2, parameter], [parameter, 2 + parameter**2]])
    transformed = B.T @ A @ B
    jacobian = det(B, check_finite=True)
    assert (
        abs(
            det(transformed, check_finite=True)
            / (jacobian**2 * det(A, check_finite=True))
            - 1
        )
        < 1e-11
    )
    old_inverse = solve(A, np.eye(56), assume_a="sym", check_finite=True)
    new_inverse = solve(transformed, np.eye(56), assume_a="sym", check_finite=True)
    assert np.linalg.norm(B @ new_inverse @ B.T - old_inverse) < 1e-10
    assert (
        abs(
            math.sqrt(det(transformed, check_finite=True))
            / (abs(jacobian) * math.sqrt(det(A, check_finite=True)))
            - 1
        )
        < 1e-11
    )
    eigenvalues = eigvalsh(transformed, check_finite=True)
    assert (sum(eigenvalues > 0), sum(eigenvalues < 0)) == (26, 30)
