"""Independent constrained-action, shear vertices, state and exact-response fixtures."""

from functools import cache
from itertools import permutations

import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_shear_response import (
    audit,
    control,
    response,
    vertices,
)
from p8_vector_hadamard import preparation
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag, expm

J = np.block([[np.zeros((3, 3)), np.eye(3)], [-np.eye(3), np.zeros((3, 3))]])


def geometry(time):
    return (1 + time * time) ** 2, 4 * time / (1 + time * time)


def cross(k):
    x, y, z = k
    return np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])


def metric_matrix(time, k, gamma, mass=1000.0):
    a, _ = geometry(time)
    E = expm(gamma)
    Q = expm(-gamma)
    C = cross(k)
    K = E / a + np.outer(k, k) / (a**3 * mass**2)
    V = a * mass**2 * Q + C.T @ E @ C / a
    return block_diag(V, K)


def first_matrix(time, k, G, mass=1000.0):
    a, _ = geometry(time)
    C = cross(k)
    return block_diag(-a * mass**2 * G + C.T @ G @ C / a, G / a)


def second_matrix(time, k, G, D, mass=1000.0):
    a, _ = geometry(time)
    C = cross(k)
    P = (G @ D + D @ G) / 2
    return block_diag(a * mass**2 * P + C.T @ P @ C / a, P / a)


def square_roots(M):
    eigen, vectors = np.linalg.eigh(M)
    assert min(eigen) > 0
    return (vectors * np.sqrt(eigen)) @ vectors.T, (
        vectors / np.sqrt(eigen)
    ) @ vectors.T


def shear(rng):
    G = rng.normal(size=(3, 3))
    G = (G + G.T) / 2
    G -= np.eye(3) * np.trace(G) / 3
    return G / np.linalg.norm(G, 2)


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_exact_residual(name):
    value = audit.residuals()[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    ), name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_unsupported_scope_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("name", list(audit.packets()))
def test_continuous_proof_packets(name):
    assert all(audit.packets()[name].get("gates", {}).values())


def test_scope_counts_and_gates():
    assert len(audit.residuals()) == 43 and audit.scalar_entry_count() == 348
    assert len(audit.gates()) == 25 and audit.rejected_inputs() == 110
    assert len(audit.controls()) == 9 and all(audit.gates().values())
    assert audit.validate_scope(audit.frontier(), audit.matching())


@pytest.mark.parametrize("fixture", range(8))
def test_independent_anisotropic_Legendre_transform_and_constraint(fixture):
    rng = np.random.default_rng(11000 + fixture)
    gamma = 0.2 * shear(rng)
    k = rng.normal(size=3) * (100 + fixture * 300)
    A = (rng.normal(size=3) + 1j * rng.normal(size=3)) / 1000
    pi = rng.normal(size=3) + 1j * rng.normal(size=3)
    time = -0.4 + fixture / 10
    a, _ = geometry(time)
    E, Q = expm(gamma), expm(-gamma)
    inverse = Q / a**2
    A0 = -1j * np.dot(k, pi) / (a**3 * 1000**2)
    velocity = E @ pi / a + 1j * k * A0
    electric = velocity - 1j * k * A0
    F = 1j * (np.outer(k, A) - np.outer(A, k))
    magnetic = sum(
        inverse[i, l] * inverse[j, n] * np.conj(F[i, j]) * F[l, n]
        for i in range(3)
        for j in range(3)
        for l in range(3)
        for n in range(3)
    )
    L = a**3 * (
        np.vdot(electric, inverse @ electric) / 2
        - magnetic / 4
        + 1000**2 * abs(A0) ** 2 / 2
        - 1000**2 * np.vdot(A, inverse @ A) / 2
    )
    independent = np.vdot(pi, velocity).real - L.real
    Z = np.concatenate([A, pi])
    actual = np.vdot(Z, metric_matrix(time, k, gamma) @ Z).real / 2
    assert np.isclose(independent, actual, rtol=1e-11, atol=1e-10)
    assert np.isclose(np.linalg.det(E), 1, rtol=1e-12)
    K = metric_matrix(time, k, gamma)[3:, 3:]
    V = metric_matrix(time, k, gamma)[:3, :3]
    omega2 = 1000**2 + np.dot(k, Q @ k) / a**2
    assert np.linalg.norm(K @ V - omega2 * np.eye(3)) < 1e-10 * omega2


@pytest.mark.parametrize("fixture", range(10))
def test_full_exponential_and_metric_vertex_relative_majorants(fixture):
    rng = np.random.default_rng(12000 + fixture)
    time = -0.4 + fixture * 0.08
    k = rng.normal(size=3) * (10 + fixture * 500)
    G, D = shear(rng), shear(rng)
    base = metric_matrix(time, k, np.zeros((3, 3)))
    _, inv = square_roots(base)
    d = 0.1 + fixture / 100
    ratio = inv @ metric_matrix(time, k, d * G) @ inv
    eigen = np.linalg.eigvalsh(ratio)
    assert min(eigen) >= np.exp(-d) - 1e-10 and max(eigen) <= np.exp(d) + 1e-10
    assert np.linalg.norm(inv @ first_matrix(time, k, G) @ inv, 2) <= 1 + 1e-10
    assert np.linalg.norm(inv @ second_matrix(time, k, G, D) @ inv, 2) <= 1 + 1e-10
    step = 1e-5
    numeric = (metric_matrix(time, k, step * G) - metric_matrix(time, k, -step * G)) / (
        2 * step
    )
    assert np.linalg.norm(inv @ (numeric - first_matrix(time, k, G)) @ inv) < 1e-8


@pytest.mark.parametrize("degree", (1, 2, 3, 4))
def test_noncommuting_all_order_exponential_vertex_bound(degree):
    rng = np.random.default_rng(13000 + degree)
    directions = [shear(rng) for _ in range(degree)]
    product = sum(
        (
            np.linalg.multi_dot(order)
            if degree > 2
            else order[0] @ order[1]
            if degree == 2
            else order[0]
        )
        for order in permutations(directions)
    ) / s.factorial(degree)
    product = np.asarray(product, dtype=float)
    assert np.allclose(product, product.T)
    assert np.linalg.norm(product, 2) <= 1 + 1e-12


@pytest.mark.parametrize("fixture", range(5))
def test_independent_mixed_second_derivative_and_nonzero_contact(fixture):
    rng = np.random.default_rng(14000 + fixture)
    k = rng.normal(size=3) * 600
    G, D = shear(rng), shear(rng)
    time = 0.1
    base = metric_matrix(time, k, np.zeros((3, 3)))
    _, inv = square_roots(base)
    step = 2e-4
    numeric = (
        metric_matrix(time, k, step * (G + D))
        - metric_matrix(time, k, step * (G - D))
        - metric_matrix(time, k, step * (-G + D))
        + metric_matrix(time, k, -step * (G + D))
    ) / (4 * step**2)
    exact = second_matrix(time, k, G, D)
    assert np.linalg.norm(inv @ (numeric - exact) @ inv) < 2e-7
    contact = second_matrix(time, k, G, G)
    # A positive finite-mode covariance makes this diagonal second-shear contact nonzero.
    assert np.trace(inv @ contact @ inv) > 0.1


def test_generic_shear_mixes_transverse_and_longitudinal_modes():
    k = np.array([0.0, 0.0, 1000.0])
    G = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    M = first_matrix(0.0, k, G)
    assert abs(M[0, 2]) > 0 and abs(M[3, 5]) > 0
    assert "mixes" in vertices.data()["state_boundary"]


def pulse(time):
    return np.exp(1 - 1 / (1 - (time / 0.4) ** 2)) if abs(time) < 0.4 else 0.0


@cache
def prepared(nu):
    a, H = geometry(-0.5)
    n = np.array([2.0, 3.0, 5.0])
    n /= np.linalg.norm(n)
    k = a * np.sqrt(nu**2 - 1000**2) * n
    e1 = np.cross(n, [0.0, 0.0, 1.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    modes = []
    for kind, direction in (
        ("transverse", e1),
        ("transverse", e2),
        ("longitudinal", n),
    ):
        packet = preparation.initial_data(kind, 1000, nu)
        W = float(packet["all_order_frequency"].evalf(40))
        rate = float(packet["all_order_half_log_rate"].evalf(40))
        z = 1 - 1000**2 / nu**2
        d = H / 2 if kind == "transverse" else H * (0.5 + z)
        g = np.sqrt(a) if kind == "transverse" else np.sqrt(a) * 1000 / nu
        f = 1 / np.sqrt(2 * W)
        p = (-1j * W - rate - d) * f
        modes.append(np.concatenate([direction * f / g, direction * g * p]))
    F = np.array(modes).T
    M = metric_matrix(-0.5, k, np.zeros((3, 3)))
    root, inverse = square_roots(M)
    return k, F, root, inverse


@cache
def actual_response(nu):
    k, F, root, inverse = prepared(nu)
    G = np.diag([1.0, -1.0, 0.0])
    epsilon = 1e-8
    initial = root @ F / np.sqrt(nu)
    cross_direction = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    plane = np.diag([1.0, 1.0, 0.0])
    Ck = cross(k)

    def matrices(time):
        a, _ = geometry(time)
        f = pulse(time)
        angle = 3 * time if nu == 1400 else 0.0
        direction = np.cos(angle) * G + np.sin(angle) * cross_direction
        even = 2 * np.sinh(epsilon * f / 2) ** 2 * plane
        E = np.eye(3) + np.sinh(epsilon * f) * direction + even
        Q = np.eye(3) - np.sinh(epsilon * f) * direction + even
        K0 = np.eye(3) / a + np.outer(k, k) / (a**3 * 1000**2)
        V0 = a * 1000**2 * np.eye(3) + Ck.T @ Ck / a
        K = E / a + np.outer(k, k) / (a**3 * 1000**2)
        V = a * 1000**2 * Q + Ck.T @ E @ Ck / a
        base = root @ J @ block_diag(V0, K0) @ inverse
        exact = root @ J @ block_diag(V, K) @ inverse
        first = root @ J @ first_matrix(time, k, direction) @ inverse
        return base, exact, first

    def equation(time, flat):
        Y = flat.reshape(6, 12)
        S, F1, Fe = Y[:, :6], Y[:, 6:9], Y[:, 9:]
        base, exact, first = matrices(time)
        return np.column_stack(
            [base @ S, base @ F1 + pulse(time) * first @ S @ initial, exact @ Fe]
        ).ravel()

    Y0 = np.column_stack([np.eye(6, dtype=complex), np.zeros((6, 3), complex), initial])
    maximum = np.sqrt(1000**2 + np.dot(k, k))
    result = solve_ivp(
        equation,
        (-0.5, 0.0),
        Y0.ravel(),
        method="DOP853",
        rtol=5e-12,
        atol=2e-13,
        max_step=1 / (4 * maximum),
    )
    assert result.success
    Y = result.y[:, -1].reshape(6, 12)
    return k, G, epsilon, initial, Y, root, inverse


@pytest.mark.parametrize("nu", (1000, 1200, 1400))
def test_actual_all_order_initial_state_CCR_and_energy(nu):
    _, F, root, _ = prepared(nu)
    W = F @ F.conj().T
    assert np.linalg.norm(W - W.T - 1j * J) < 1e-10
    assert min(np.linalg.eigvalsh(W)) > -1e-8
    trace = np.trace(root @ W.real @ root)
    assert trace < 108 * nu
    assert "all-order" in response.data()["actual_initial_state"]


@pytest.mark.parametrize("nu", (1000, 1200, 1400))
def test_actual_finite_amplitude_covariance_remainder(nu):
    k, _, epsilon, initial, Y, _, _ = actual_response(nu)
    S, F1, Fe = Y[:, :6], Y[:, 6:9], Y[:, 9:]
    F1I = np.linalg.solve(S, F1)
    FeI = np.linalg.solve(S, Fe)
    Cinit = (initial @ initial.conj().T).real
    first = (F1I @ initial.conj().T + initial @ F1I.conj().T).real
    remainder = (FeI @ FeI.conj().T).real - Cinit - epsilon * first
    norm = np.sum(np.abs(np.linalg.eigvalsh((remainder + remainder.T) / 2)))
    maximum = np.sqrt(1000**2 + np.dot(k, k))
    A = 16 * maximum * np.expm1(epsilon)
    R = 8 * maximum * epsilon**2 * np.exp(epsilon) + A * A * np.exp(A) / 2
    bound = np.trace(Cinit) * (2 * R + np.expm1(A) ** 2)
    assert maximum <= 2000 and norm < bound
    assert bound / np.trace(Cinit) < 1e-6
    assert np.linalg.norm(S, 2) < 4 and np.linalg.norm(np.linalg.inv(S), 2) < 4


@pytest.mark.parametrize("nu", (1000, 1200, 1400))
def test_actual_current_tangent_includes_metric_contact(nu):
    k, G, epsilon, initial, Y, _, inverse = actual_response(nu)
    S, F1, Fe = Y[:, :6], Y[:, 6:9], Y[:, 9:]
    F0 = inverse @ S @ initial * np.sqrt(nu)
    dF = inverse @ F1 * np.sqrt(nu)
    actual = inverse @ Fe * np.sqrt(nu)
    C0 = (F0 @ F0.conj().T).real
    dC = (dF @ F0.conj().T + F0 @ dF.conj().T).real
    Ce = (actual @ actual.conj().T).real
    first = first_matrix(0, k, G)
    contact = second_matrix(0, k, G, G)
    exact_contact = -np.trace(contact @ C0) / 2
    tangent = -np.trace(first @ dC) / 2 + exact_contact
    a, _ = geometry(0)
    E, Q = expm(epsilon * G), expm(-epsilon * G)
    Ck = cross(k)
    readout = block_diag(-a * 1000**2 * Q @ G + Ck.T @ E @ G @ Ck / a, E @ G / a)
    finite = (-np.trace(readout @ Ce) / 2 + np.trace(first @ C0) / 2) / epsilon
    assert abs(exact_contact) > 1
    assert abs(finite - tangent) < 0.005 * abs(exact_contact)
    assert abs(finite - (-np.trace(first @ dC) / 2)) > 0.9 * abs(exact_contact)


@pytest.mark.parametrize("fixture", range(3))
def test_independent_three_oscillator_Fock_Wick_and_Kubo_identity(fixture):
    # The vacuum products reach at most occupation2, so three levels per factor suffice here.
    ann = np.diag(np.sqrt(np.arange(1, 3)), 1)
    eye = np.eye(3)
    annihilators = [
        np.kron(np.kron(ann, eye), eye),
        np.kron(np.kron(eye, ann), eye),
        np.kron(np.kron(eye, eye), ann),
    ]
    q = [(a + a.T) / np.sqrt(2) for a in annihilators]
    p = [-1j * (a - a.T) / np.sqrt(2) for a in annihilators]
    rng = np.random.default_rng(15000 + fixture)
    C = rng.normal(size=(3, 3))
    C = (C + C.T) / 2
    S = np.block([[np.eye(3), np.zeros((3, 3))], [C, np.eye(3)]])
    Z = [sum(S[i, j] * op for j, op in enumerate(q + p)) for i in range(6)]
    W = np.array([[(a @ b)[0, 0] for b in Z] for a in Z])
    A = rng.normal(size=(6, 6))
    A = (A + A.T) / 2
    B = rng.normal(size=(6, 6))
    B = (B + B.T) / 2
    HA = sum(A[i, j] * Z[i] @ Z[j] / 2 for i in range(6) for j in range(6))
    HB = sum(B[i, j] * Z[i] @ Z[j] / 2 for i in range(6) for j in range(6))
    connected = (HA @ HB)[0, 0] - HA[0, 0] * HB[0, 0]
    assert np.allclose(connected, np.trace(A @ W @ B @ W.T) / 2, atol=1e-10)
    comm = (HA @ HB - HB @ HA)[0, 0]
    assert np.allclose(
        comm, 1j * np.trace((A @ J @ B - B @ J @ A) @ W.real) / 2, atol=1e-10
    )
    covariance_tangent = J @ B @ W.real - W.real @ B @ J
    assert np.allclose(1j * comm, -np.trace(A @ covariance_tangent) / 2, atol=1e-10)


def test_no_integrable_UV_claim_or_cutoff_from_mode_benchmark():
    packet = control.data()
    assert "not an integrable" in packet["ultraviolet_boundary"]
    assert "not a physical cutoff" in packet["ultraviolet_boundary"]
    assert packet["explicit_benchmark"][
        "covariance_relative_remainder_upper"
    ] < s.Rational(1, 10**6)
    assert "continuum renormalization" in response.data()["Kubo_identity"]
