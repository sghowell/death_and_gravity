"""Independent full-time, complex-adjoint, cutoff and bounded-unitary diagnostics."""

from functools import cache

import mpmath as mp
import numpy as np
import pytest
import sympy as s
from numpy.testing import assert_allclose
from p8_vacuum_affine_quantitative_local_time import (
    audit,
    branch,
    quantum,
    reference,
    source,
)
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag, expm
from scipy.special import iv


@pytest.mark.parametrize(
    "packet", [source.data, branch.data, reference.data, quantum.data]
)
def test_complete_packets_and_written_proof_gates(packet):
    data = packet()
    for name, value in data["checks"].items():
        entries = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(entry == 0 for entry in entries), name
    assert all(bool(value) for value in data["gates"].values())


@pytest.mark.parametrize("orders", [(i, j) for i in range(6) for j in range(6 - i)])
def test_all_full_mixed_real_time_complex_lapse_source_rows(orders):
    row = source.complex_source_bounds()["rows"][orders]
    assert 0 < row[0] < 10**30 and 0 < row[1] < 10**30
    assert 0 < row[2] < 1
    i, j = orders
    assert i + j <= 5
    assert (
        source.TIME + source.complex_source_bounds()["radius"]
        < source.complex_source_bounds()["distance"]
    )


@pytest.mark.parametrize("name", tuple(branch.derivative_bounds()))
def test_all_187_complete_auxiliary_majorants(name):
    value = branch.derivative_bounds()[name]
    assert value.is_Rational is True and 0 <= value < branch.MAJORANT
    assert s.ceiling(value) >= value


@pytest.mark.parametrize("i,j", [(0, 1), (1, 1), (0, 2), (2, 1), (1, 2)])
def test_whole_primitive_N_elimination_and_pure_time_integral_are_distinct(i, j):
    derivative = s.diff(source.primitive, source.u, i, source.N, j)
    assert source.eliminate_N_primitive(derivative) == s.diff(
        source.IN, source.u, i, source.N, j - 1
    )
    pure_time = s.diff(source.primitive, source.u, i + 1)
    assert source.eliminate_N_primitive(pure_time) == pure_time
    assert branch.magnitude(pure_time) == source.LAPSE_RADIUS * branch.magnitude(
        s.diff(source.IN, source.u, i + 1)
    )


def test_bounce_primitive_time_contact_cannot_be_deleted():
    whole = source.CONSTRAINT
    missing = whole.subs(s.Derivative(source.primitive, source.u), 0)
    # In C the first Iu/N contact simplifies to (Iu)_N; that exact N
    # derivative is eliminated already. Deleting IN before differentiating
    # also loses genuine second time-source contacts.
    wrong_H = source.HAMILTONIAN.subs(
        {source.primitive: 0, s.Derivative(source.primitive, source.u): 0}
    )
    wrong_C = source.eliminate_N_primitive(s.diff(wrong_H, source.N))
    assert s.simplify(whole - wrong_C) != 0
    assert branch.bounce_identity() == 0
    assert s.simplify(missing - whole) == 0


@cache
def complex_Galerkin_ghost(cutoff, amplitude):
    # Full smooth coefficient diag(1,exp(a cos x),exp(-a cos x)).
    # Its exact Bessel coefficients are compressed for this diagnostic only.
    modes = [
        (i, j, 0)
        for i in range(-cutoff, cutoff + 1)
        for j in range(-cutoff, cutoff + 1)
        if (i, j) != (0, 0)
    ]
    M = np.zeros((3 * len(modes), 3 * len(modes)), dtype=complex)
    formal = np.zeros_like(M)
    missing = np.zeros_like(M)
    flat = []
    for mode in modes:
        p = np.asarray(mode, dtype=float)
        flat.append(np.dot(p, p) * np.eye(3) + np.outer(p, p) / 3)
    for row, mode in enumerate(modes):
        k = np.asarray(mode, dtype=float)
        for column, input_mode in enumerate(modes):
            p = np.asarray(input_mode, dtype=float)
            n = k - p
            if n[1] or n[2]:
                continue
            order = int(n[0])
            Q = np.diag(
                [
                    float(order == 0),
                    iv(abs(order), amplitude),
                    (-1.0) ** order * iv(abs(order), amplitude),
                ]
            )
            a, b = slice(3 * row, 3 * row + 3), slice(3 * column, 3 * column + 3)
            M[a, b] = np.dot(p, Q @ p) * np.eye(3) + np.outer(Q @ p, p) / 3
            formal[a, b] = np.dot(k, Q @ k) * np.eye(3) + np.outer(k, Q @ k) / 3
            missing[a, b] = np.dot(p, Q @ p) * np.eye(3) + np.outer(p, Q @ p) / 3
    return M, formal, missing, block_diag(*flat)


@pytest.mark.parametrize("cutoff", [1, 2])
@pytest.mark.parametrize("amplitude", [5e-5 + 7e-5j, -3e-5 + 9e-5j, 1e-4 - 1e-4j])
def test_holomorphic_formal_adjoint_not_complex_conjugate_adjoint(cutoff, amplitude):
    M, formal, missing, flat = complex_Galerkin_ghost(cutoff, amplitude)
    # The coefficient is even in space, so Fourier mode reversal cancels
    # here; a general Fourier bilinear adjoint additionally reverses modes.
    assert_allclose(formal, M.T, atol=1e-14)
    assert np.linalg.norm(formal - M.conj().T, 2) > abs(amplitude.imag) / 10
    assert np.linalg.norm(missing - formal, 2) > abs(amplitude) / 10
    rng = np.random.default_rng(270 + cutoff)
    left = rng.normal(size=len(M)) + 1j * rng.normal(size=len(M))
    right = rng.normal(size=len(M)) + 1j * rng.normal(size=len(M))
    assert_allclose(left @ (M @ right), (formal @ left) @ right, atol=3e-11)
    inverse0 = np.linalg.inv(flat)
    perturbation = (formal - flat) @ inverse0
    budget = 100 * np.expm1(4 * abs(amplitude))
    assert np.linalg.norm(perturbation, 2) < budget < 0.1
    partial, term = np.eye(len(M), dtype=complex), np.eye(len(M), dtype=complex)
    for _ in range(8):
        term = -perturbation @ term
        partial += term
    assert_allclose(formal @ (inverse0 @ partial), np.eye(len(M)), atol=4e-13)


@pytest.mark.parametrize("case", range(8))
def test_correlated_pure_covariance_whitening_inverse_and_heat_operator(case):
    rng = np.random.default_rng(27000 + case)
    omega = np.block([[np.zeros((4, 4)), np.eye(4)], [-np.eye(4), np.zeros((4, 4))]])
    raw = rng.normal(size=(8, 8))
    S = expm(0.1 * omega @ (raw + raw.T))
    V = S @ S.T / 2
    inverse = np.linalg.inv(S)
    assert_allclose(S @ omega @ S.T, omega, atol=1e-13)
    assert_allclose(V @ omega @ V, omega / 4, atol=1e-13)
    assert_allclose(np.linalg.norm(inverse, 2), np.linalg.norm(S, 2), atol=1e-13)
    assert_allclose(inverse @ V @ inverse.T / 2, np.eye(8) / 4, atol=1e-13)
    assert np.linalg.norm(V - np.diag(np.diag(V))) > 0.01


@pytest.mark.parametrize("name", tuple(source.COORDS))
def test_every_actual_off_time_invariant_and_original_domain_margin(name):
    values = reference.phase_domain()["invariants"]
    assert 0 < values[name] < source.IMAGE
    assert source.IMAGE < source.bounce.DELTA / 2


def transition(y):
    if y <= 0:
        return mp.mpf(1)
    if y >= 1:
        return mp.mpf(0)
    left, right = mp.exp(-1 / y), mp.exp(-1 / (1 - y))
    return right / (left + right)


@pytest.mark.parametrize(
    "point",
    [
        s.Rational(1, 100),
        s.Rational(1, 10),
        s.Rational(1, 3),
        s.Rational(1, 2),
        s.Rational(4, 5),
        s.Rational(99, 100),
    ],
)
@pytest.mark.parametrize("order", range(5))
def test_independent_full_transition_derivatives_against_exact_ceilings(point, order):
    with mp.workdps(80):
        y = mp.mpf(int(s.numer(point))) / int(s.denom(point))
        value = mp.diff(transition, y, order)
        ceiling = int(quantum.cutoff_bounds()["transition"][order])
        assert mp.isfinite(value) and abs(value) <= ceiling
        denominator = mp.exp(-1 / y) + mp.exp(-1 / (1 - y))
        assert denominator > mp.mpf(1) / 9


@pytest.mark.parametrize("width", [1, 2])
@pytest.mark.parametrize(
    "point", [s.Rational(1, 10), s.Rational(2, 5), s.Rational(4, 5)]
)
@pytest.mark.parametrize("order", range(5))
def test_full_radial_cutoff_derivatives_and_support(width, point, order):
    with mp.workdps(80):
        y = mp.mpf(int(s.numer(point))) / int(s.denom(point))
        radial = mp.sqrt(1 + width * y)
        value = mp.diff(lambda r: transition((r * r - 1) / width), radial, order)
        assert abs(value) <= int(quantum.cutoff_bounds()["radial"][order])
        assert transition((mp.mpf("0.9") ** 2 - 1) / width) == 1
        assert transition((mp.mpf(2) ** 2 - 1) / width) == 0
    assert 1 + width < 4


@pytest.mark.parametrize("order", range(5))
def test_individual_Cauchy_polydisc_and_every_Leibniz_contact(order):
    # Only the coordinates in this derivative are varied. At total order4
    # there are at most4, even though the physical phase dimension is96.
    active = max(1, order)
    assert 2 + s.sqrt(active) / 4 < 4
    assert s.sqrt(96) / 4 + 2 > 4  # Reject the incorrect all-coordinate premise.
    raw = quantum.cutoff_bounds()
    expected = sum(
        s.binomial(order, k)
        * raw["radial"][k]
        * s.factorial(order - k)
        * 4 ** (order - k)
        for k in range(order + 1)
    )
    assert raw["product"][order] == expected
    if order >= 2:
        endpoints = raw["radial"][order] + s.factorial(order) * 4**order
        assert expected > endpoints


def test_mixed_radial_and_heat_contacts_are_not_omittable():
    x, y = s.symbols("x y", real=True)
    phi = s.Function("transition")
    a = x * x + y * y - 1
    literal = s.diff(phi(a), x, 2, y, 2)
    t = s.Symbol("argument")
    expected = (
        16 * x * x * y * y * s.diff(phi(t), t, 4)
        + 8 * (x * x + y * y) * s.diff(phi(t), t, 3)
        + 4 * s.diff(phi(t), t, 2)
    ).subs(t, a)
    assert s.simplify(literal - expected) == 0
    assert (
        s.simplify(literal - 16 * x * x * y * y * s.diff(phi(t), t, 4).subs(t, a)) != 0
    )
    chi = x * x * y + y * y
    f = x * y + x**3
    lap = lambda value: (s.diff(value, x, 2) + s.diff(value, y, 2)) / 4
    cross = lap(chi * f) - chi * lap(f) - f * lap(chi)
    assert (
        s.expand(
            cross - (s.diff(chi, x) * s.diff(f, x) + s.diff(chi, y) * s.diff(f, y)) / 2
        )
        == 0
    )
    assert s.expand(cross) != 0


@pytest.mark.parametrize("sign", [-1, 1])
def test_actual_heavy_frozen_coefficient_diagnostic_resolves_tiny_time_and_wrong_balance(
    sign,
):
    # Frozen bounce-coefficient KG oscillator only: this is NOT claimed to
    # be a numerical evolution of the full time-dependent P8 Hamiltonian.
    with mp.workdps(4200):
        P = mp.mpf(10) ** 64
        mass2 = mp.mpf(10) ** 200 / 512 + 2
        omega = mp.sqrt(mass2 + P * P)
        balance = mp.sqrt(mass2 + P * P / 16)
        T = sign * mp.mpf(10) ** -2000
        cosine, sine = mp.cos(omega * T), mp.sin(omega * T)
        assert 0 < 1 - cosine
        assert (
            mp.mpf("0.49") * (omega * T) ** 2
            < 1 - cosine
            < mp.mpf("0.51") * (omega * T) ** 2
        )
        wrong = mp.cos(balance * T)
        difference = wrong - cosine
        expected = (15 * P * P / 16) * T * T / 2
        assert mp.mpf("0.99") * expected < difference < mp.mpf("1.01") * expected
        matrix = mp.matrix([[cosine, sine / omega], [-omega * sine, cosine]])
        assert abs(mp.det(matrix) - 1) < mp.mpf(10) ** -4100
        assert abs(matrix[1, 0]) < 2 * omega**2 * abs(T)
        assert mp.sign(matrix[0, 1]) == sign
        assert mp.mpf(10) ** 98 < balance < omega < mp.mpf(10) ** 100


@pytest.mark.parametrize("case", range(8))
def test_full_noncommuting_Duhamel_leakage_and_two_readout_comparison(case):
    # Finite-dimensional bounded-operator fixture, not a replacement for
    # the infinite-dimensional CCR or an actual P8 mode evolution.
    rng = np.random.default_rng(270000 + case)
    raw0 = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    raw1 = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    H0, H1 = (raw0 + raw0.conj().T) / 8, (raw1 + raw1.conj().T) / 8
    H0 += 0.37 * np.eye(3)  # Keep the full scalar background phase.
    delta = np.zeros((3, 3), dtype=complex)
    delta[1:, 1:] = np.array([[0.07, 0.02j], [-0.02j, -0.04]])
    assert np.linalg.norm(H0 @ H1 - H1 @ H0) > 0.001
    tail, T = 0.0004, 0.01
    seed = np.array([np.sqrt(1 - tail), np.sqrt(tail), 0], dtype=complex)
    identity = np.eye(3, dtype=complex)

    def solve(extra):
        def rhs(t, flat):
            U = flat.reshape(3, 3)
            return (-1j * (H0 + t * H1 + extra) @ U).ravel()

        solution = solve_ivp(
            rhs, (0, T), identity.ravel(), method="DOP853", rtol=2e-13, atol=2e-14
        )
        assert solution.success
        U = solution.y[:, -1].reshape(3, 3)
        assert_allclose(U.conj().T @ U, identity, atol=1e-12)
        return U @ seed

    psi1, psi2 = solve(delta), solve(np.zeros((3, 3), dtype=complex))
    B = np.linalg.norm(H0, 2) + T * np.linalg.norm(H1, 2) + np.linalg.norm(delta, 2)
    M = np.linalg.norm(delta, 2)
    state_bound = M * np.sqrt(tail) * T + M * B * T * T / 2
    assert np.linalg.norm(psi1 - psi2) <= state_bound + 2e-12
    assert np.linalg.norm(psi1 - seed) <= B * T + 2e-12
    outside = np.diag([0.0, 1.0, 1.0])
    assert np.vdot(psi1, outside @ psi1).real <= tail + 2 * B * T + 2e-12
    F2 = np.diag([1.01, 0.99, 1.02])
    difference = np.diag([0.0, 0.006, -0.004])
    F1 = F2 + difference
    assert np.linalg.eigvalsh(F1).min() > 0.9
    assert np.linalg.eigvalsh(F2).min() > 0.9
    actual = abs(np.vdot(psi1, F1 @ psi1) - np.vdot(psi2, F2 @ psi2))
    volume_bound = np.linalg.norm(difference, 2)
    complete = 4 * state_bound + volume_bound * tail + 2 * volume_bound * B * T
    assert actual <= complete + 3e-12
    # Identical physical-picture transport cancels on BOTH state and readout.
    reference_U = expm(-1j * 0.19 * H1)
    physical = reference_U @ psi1
    observable = reference_U @ F1 @ reference_U.conj().T
    assert_allclose(
        np.vdot(physical, observable @ physical), np.vdot(psi1, F1 @ psi1), atol=1e-12
    )


def test_actual_quantum_bounds_remain_exact_and_do_not_underflow():
    data = quantum.dynamical_bounds()
    assert all(value.is_Rational is True and value > 0 for value in data.values())
    assert data["whole_same_seed_unitary_state_error"] == s.Rational(1, 10**990)
    assert data["whole_two_explicit_regulator_state_difference"] < s.Rational(
        1, 10**1970
    )
    assert data["whole_two_readout_mean_difference"] < s.Rational(1, 10**1230)
    cut = quantum.cutoff_bounds()
    assert cut["Dvolume"] > 0 and cut["D2volume"] > 0
    assert 1 - quantum.VERROR - cut["Dvolume"] > s.Rational(1, 2)
    with mp.workdps(100):
        log_tail = -(mp.mpf(10) ** 39)
        assert 0 < mp.exp(log_tail) < mp.mpf(10) ** -4000


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[x[0] for x in audit.bad_cases()]
)
def test_unsupported_source_time_cutoff_and_closure_inputs_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_all_original_frontiers_full_family_and_two_declared_cutoffs_retained():
    assert audit.require_parameters(audit.parameters()) == audit.parameters()
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 126
    assert audit.require_radius(0) == 0
    assert audit.require_radius(2 * reference.CORE) == 2 * reference.CORE
    for T in (-source.TIME, 0, source.TIME):
        assert audit.require_time(T) == T
    for width in (1, 2):
        assert audit.require_cutoff(width) == width
    assert audit.require_dimension(48) == 48
    assert audit.require_momentum(reference.field.P) == reference.field.P
