"""Independent inhomogeneous Hamiltonian, Fourier, contact and transport checks."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_spatial_current import (
    audit,
    bounds,
    hamiltonian,
    response,
    vertices,
)


def mp_matrix(M):
    return mp.matrix([[mp.mpf(str(s.N(x, 100))) for x in row] for row in M.tolist()])


def trace(M):
    return sum(M[i, i] for i in range(M.rows))


def kron(A, B):
    return mp.matrix(
        [
            [
                A[i // B.rows, j // B.cols] * B[i % B.rows, j % B.cols]
                for j in range(A.cols * B.cols)
            ]
            for i in range(A.rows * B.rows)
        ]
    )


def block(*matrices):
    size = sum(M.rows for M in matrices)
    out = mp.zeros(size)
    j = 0
    for M in matrices:
        out[j : j + M.rows, j : j + M.cols] = M
        j += M.rows
    return out


def cross(k):
    return mp.matrix([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])


def opnorm(M):
    _, singular, _ = mp.svd(M)
    return max(abs(x) for x in singular)


@cache
def grid():
    with mp.workdps(100):
        N = 5
        positions = [2 * mp.pi * j / N for j in range(N)]
        momenta = [0, 1, 2, -2, -1]
        U = mp.matrix(
            [[mp.exp(-mp.j * k * x) / mp.sqrt(N) for x in positions] for k in momenta]
        )
        derivative = U.H * mp.diag([mp.j * k for k in momenta]) * U
        derivative = mp.matrix(
            [[mp.re(derivative[i, j]) for j in range(N)] for i in range(N)]
        )
        ez = mp.matrix([0, 0, 1])
        grad = kron(derivative, ez)
        curl = kron(derivative, cross(ez))
        F0 = mp.matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 100
        F1 = mp.matrix([[2, -1, 1], [-1, 0, 2], [1, 2, -2]]) / 100
        F2 = mp.matrix([[0, 1, 2], [1, 1, -1], [2, -1, -1]]) / 100
        G = [F0 + mp.cos(x) * F1 + mp.sin(x) * F2 for x in positions]
        D = [F1 + mp.sin(x) * F0 - mp.cos(x) * F2 for x in positions]
        return N, positions, momenta, U, grad, curl, G, D


def lattice_H(e, b):
    N, _positions, _momenta, _U, grad, curl, G, D = grid()
    a = mp.mpf(25) / 16
    m = mp.mpf(1000)
    E = block(*(mp.expm(e * G[j] + b * D[j]) for j in range(N)))
    Q = block(*(mp.expm(-e * G[j] - b * D[j]) for j in range(N)))
    K = E / a + grad * grad.T / (a**3 * m * m)
    V = a * m * m * Q + curl * E * curl / a
    return V, K, E, Q


@pytest.mark.parametrize("epsilon", (-mp.mpf("0.01"), mp.mpf(0), mp.mpf("0.01")))
def test_independent_spatial_Legendre_transform_and_temporal_constraint(epsilon):
    with mp.workdps(85):
        N, _, _, _, grad, curl, _, _ = grid()
        a = mp.mpf(25) / 16
        m = mp.mpf(1000)
        V, K, E, Q = lattice_H(epsilon, mp.mpf(0))
        A = mp.matrix([mp.mpf((j % 7) - 3) / 23 for j in range(3 * N)])
        pi = mp.matrix([mp.mpf((j % 5) - 2) / 19 for j in range(3 * N)])
        divpi = -grad.T * pi
        A0 = -divpi / (a**3 * m * m)
        velocity = E * pi / a + grad * A0
        electric = velocity - grad * A0
        L = (
            a * (electric.T * Q * electric)[0] / 2
            + a**3 * m * m * (A0.T * A0)[0] / 2
            - a * m * m * (A.T * Q * A)[0] / 2
            - ((curl * A).T * E * (curl * A))[0] / (2 * a)
        )
        literal = (pi.T * velocity)[0] - L
        expected = ((A.T * V * A)[0] + (pi.T * K * pi)[0]) / 2
        assert abs(literal - expected) < mp.mpf("1e-60")
        assert mp.norm(-divpi - a**3 * m * m * A0) < mp.mpf("1e-70")
        assert mp.norm(curl - curl.T) < mp.mpf("1e-70")
        assert mp.norm(curl) > 1 and mp.norm(grad) > 1
        assert mp.norm(E * Q - mp.eye(3 * N)) < mp.mpf("1e-70")


@cache
def differentiated_grid(kind):
    with mp.workdps(95):
        _N, _, _, U, _, _, _, _ = grid()
        transform = kron(U, mp.eye(3))
        derivative = (0, 1) if kind == "first" else (1, 1)
        raw = [
            mp.diff(
                lambda e, b, index=index: lattice_H(e, b)[index],
                (mp.mpf(0), mp.mpf(0)),
                derivative,
            )
            for index in (0, 1)
        ]
        return [transform * M * transform.H for M in raw]


@pytest.mark.parametrize("ki", (0, 1, -1))
@pytest.mark.parametrize("qi", (0, 1, -1))
@pytest.mark.parametrize("kind", ("first", "contact"))
def test_independent_real_lattice_DFT_two_momentum_vertices(ki, qi, kind):
    with mp.workdps(80):
        N, positions, momenta, _, _, _, G, D = grid()
        kslot, qslot = momenta.index(ki), momenta.index(qi)
        p = ki - qi
        local = (
            D
            if kind == "first"
            else [(D[j] * G[j] + G[j] * D[j]) / 2 for j in range(N)]
        )
        coefficient = (
            sum(
                (mp.exp(-mp.j * p * positions[j]) * local[j] for j in range(N)),
                mp.zeros(3),
            )
            / N
        )
        a = mp.mpf(25) / 16
        m = mp.mpf(1000)
        k, q = mp.matrix([0, 0, ki]), mp.matrix([0, 0, qi])
        predicted = [
            (-1 if kind == "first" else 1) * a * m * m * coefficient
            + cross(k).T * coefficient * cross(q) / a,
            coefficient / a,
        ]
        direct = differentiated_grid(kind)
        for index in (0, 1):
            actual = direct[index][3 * kslot : 3 * kslot + 3, 3 * qslot : 3 * qslot + 3]
            assert mp.norm(actual - predicted[index]) < mp.mpf("1e-55") * max(
                1, mp.norm(predicted[index])
            )
        assert any(mp.norm(G[j] * D[j] - D[j] * G[j]) > 0 for j in range(N))


def balanced(k):
    a = mp.mpf(25) / 16
    m = mp.mpf(1000)
    k = mp.matrix(k)
    k2 = (k.T * k)[0]
    omega = mp.sqrt(m * m + k2 / a**2)
    P = k * k.T / k2 if k2 else mp.zeros(3)
    B = (mp.eye(3) - P + (omega / m) * P) / mp.sqrt(a)
    inverse = mp.sqrt(a) * (mp.eye(3) - P + (m / omega) * P)
    return omega, block(B / mp.sqrt(omega), mp.sqrt(omega) * inverse)


@pytest.mark.parametrize("which", (0, 1, 2, 3, 4, 5))
@pytest.mark.parametrize("contact", (False, True))
def test_all_momentum_energy_relative_vertex_and_contact_norm(which, contact):
    with mp.workdps(90):
        cases = [
            ([0, 0, 0], [0, 0, 0]),
            ([0, 0, 0], [1, 2, 3]),
            ([1000, 0, 0], [0, 1000, 0]),
            ([1, 2, 3], [-2, 1, 2]),
            ([10**20, 0, 0], [0, 1, 0]),
            ([10**20, 2 * 10**20, 3 * 10**20], [-2 * 10**20, 10**20, 2 * 10**20]),
        ]
        kr, qr = cases[which]
        k, q = s.Matrix(kr), s.Matrix(qr)
        D = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 10
        G = s.Matrix([[2, -1, 1], [-1, 0, 2], [1, 2, -2]]) / 10
        H = (D * G + G * D) / 2 if contact else D
        module = vertices.metric_contact if contact else vertices.metric_vertex
        raw = module(k, q, H, a=s.Rational(25, 16), m=s.Integer(1000))
        wk, Tk = balanced(kr)
        wq, Tq = balanced(qr)
        actual = Tk.T * mp_matrix(raw) * Tq / mp.sqrt(wk * wq)
        expected = opnorm(mp_matrix(D)) * (opnorm(mp_matrix(G)) if contact else 1)
        assert opnorm(actual) < expected + mp.mpf("1e-35")
        for momentum, T, w in ((k, Tk, wk), (q, Tq, wq)):
            F = mp_matrix(
                hamiltonian.features(momentum, a=s.Rational(25, 16), m=s.Integer(1000))
            )
            R = F * T / mp.sqrt(w)
            assert mp.norm(R.T * R - mp.eye(6)) < mp.mpf("1e-35")


@pytest.mark.parametrize("momentum", (1, 3, 1000))
def test_high_external_transfer_has_zero_first_vertex_but_nonzero_contact(momentum):
    with mp.workdps(70):
        # Keep only the real constant vector oscillator. Integrate the
        # actual metric's second derivative, not a projected first vertex.
        D = mp.diag([1, -1, 0])
        m = mp.mpf(1000)
        E = lambda e, x: mp.expm(e * mp.cos(momentum * x) * D)
        first = (
            mp.quad(
                lambda x: mp.diff(lambda e: E(e, x)[0, 0], mp.mpf(0)),
                [0, 2 * mp.pi / momentum],
            )
            * momentum
            / (2 * mp.pi)
        )
        second = (
            mp.quad(
                lambda x: mp.diff(lambda e: E(e, x)[0, 0], mp.mpf(0), 2),
                [0, 2 * mp.pi / momentum],
            )
            * momentum
            / (2 * mp.pi)
        )
        assert abs(first) < mp.mpf("1e-55")
        assert abs(second - mp.mpf(1) / 2) < mp.mpf("1e-55")
        C = block(mp.eye(3) / (2 * m), m * mp.eye(3) / 2)
        contact = block(m * m * D * D / 2, D * D / 2)
        assert -trace(contact * C) / 2 == -m * trace(D * D) / 4
        assert trace(contact * C) > 0


@pytest.mark.parametrize("name", ("hamiltonian", "vertices", "response", "bounds"))
def test_exact_core_packet(name):
    p = {
        "hamiltonian": hamiltonian,
        "vertices": vertices,
        "response": response,
        "bounds": bounds,
    }[name].data()
    for value in p["checks"].values():
        entries = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.cancel(x) == 0 for x in entries)
    assert all(p["gates"].values())


@pytest.mark.parametrize("which", (0, 1, 2))
@pytest.mark.parametrize("duration", ("0.1", "0.5", "1"))
def test_independent_actual_two_momentum_coupled_evolution_and_contact(which, duration):
    with mp.workdps(80):
        cases = [
            ([1, 2, 3], [-2, 1, 2]),
            ([1000, 0, 0], [0, 2000, 0]),
            ([10**16, 0, 0], [0, 1, 0]),
        ]
        kr, qr = cases[which]
        k, q = s.Matrix(kr), s.Matrix(qr)
        D = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 10
        G = s.Matrix([[2, -1, 1], [-1, 0, 2], [1, 2, -2]]) / 10
        wk, Tk = balanced(kr)
        wq, Tq = balanced(qr)
        m = s.Integer(1000)
        a = s.Rational(25, 16)
        HG = Tk.T * mp_matrix(vertices.metric_vertex(k, q, G, a, m)) * Tq
        HD = Tk.T * mp_matrix(vertices.metric_vertex(k, q, D, a, m)) * Tq
        Hlocal = (D * G + G * D) / 2
        NK = Tk.T * mp_matrix(vertices.metric_contact(k, k, Hlocal, a, m)) * Tk
        NQ = Tq.T * mp_matrix(vertices.metric_contact(q, q, Hlocal, a, m)) * Tq
        J = mp.zeros(6)
        J[:3, 3:] = mp.eye(3)
        J[3:, :3] = -mp.eye(3)
        Jall = block(J, J)
        generator = block(wk * J, wq * J)
        perturbation = mp.zeros(12)
        perturbation[:6, 6:] = HG
        perturbation[6:, :6] = HG.T
        readout = mp.zeros(12)
        readout[:6, 6:] = HD
        readout[6:, :6] = HD.T
        contact = block(NK, NQ)
        r = mp.matrix([1, -2, 1, 3, -1, 2]) / 20
        v = mp.matrix([2, 1, -1, 1, 2, -3]) / 20
        C0 = block(mp.eye(6) / 2 + r * r.T, mp.eye(6) * mp.mpf(3) / 4 + v * v.T)
        total = mp.mpf(duration) / max(wk, wq)

        def free(t):
            return block(
                mp.cos(wk * t) * mp.eye(6) + mp.sin(wk * t) * J,
                mp.cos(wq * t) * mp.eye(6) + mp.sin(wq * t) * J,
            )

        def actual(e):
            U = mp.expm((generator + e * Jall * perturbation) * total)
            return -trace((readout + e * contact) * U * C0 * U.T) / 2

        literal = mp.diff(actual, mp.mpf(0))
        points, weights = mp.gauss_quadrature(24, "legendre")
        tangent = mp.zeros(12)
        for x, w in zip(points, weights):
            time = (x + 1) * total / 2
            U = free(time)
            C = U * C0 * U.T
            source = Jall * perturbation * C + C * perturbation * Jall.T
            transport = free(total - time)
            tangent += w * total / 2 * (transport * source * transport.T)
        final_cov = free(total) * C0 * free(total).T
        predicted = -trace(readout * tangent) / 2 - trace(contact * final_cov) / 2
        assert abs(literal - predicted) < mp.mpf("1e-40") * max(1, abs(literal))
        assert abs(trace(contact * final_cov)) > mp.mpf("1e-30")
        assert wk != wq or which == 0
        assert mp.norm(tangent[:6, 6:]) > 0
        assert mp.norm(tangent[6:, :6] - tangent[:6, 6:].T) < mp.mpf("1e-40") * max(
            1, mp.norm(tangent)
        )


@pytest.mark.parametrize("K", (1000, 10**6, 10**16))
def test_complete_band_measure_and_two_separate_current_displays(K):
    c = bounds.constants()
    assert 4 * bounds.READOUT * s.Integer(K) ** 5 / 54 < s.Integer(10) ** 24 * K**5
    assert 2 * bounds.CONTACT * s.Integer(K) ** 4 / 54 < s.Integer(10) ** 13 * K**4
    assert bounds.require_band(K) == K
    assert (
        c["full_covariance_pair_source_and_transport"]
        == 2 * bounds.COV * bounds.PROP * 2
    )
    assert c["full_six_dimensional_pair_readout"] == 3 * 2 * bounds.PAIR


@pytest.mark.parametrize(
    "bad",
    (
        True,
        False,
        1.0,
        s.Float(1000),
        "1000",
        None,
        s.oo,
        -s.oo,
        s.zoo,
        s.I,
        s.nan,
        s.Symbol("K"),
        -1000,
        0,
        999,
    ),
)
def test_band_is_a_finite_exact_computational_parameter(bad):
    with pytest.raises((ValueError, TypeError)):
        bounds.require_band(bad)


def test_spatial_vertex_cannot_identify_the_two_internal_momenta():
    k, q = s.Matrix([1, 2, 3]), s.Matrix([3, -1, 2])
    D = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 10
    actual = vertices.metric_vertex(k, q, D, a=s.Integer(1), m=s.Integer(1000))
    wrong = vertices.metric_vertex(k, k, D, a=s.Integer(1), m=s.Integer(1000))
    assert actual != wrong
    assert actual[:3, :3] != actual[:3, :3].T
    assert (
        actual == vertices.metric_vertex(q, k, D, a=s.Integer(1), m=s.Integer(1000)).T
    )


def test_bounds_do_not_claim_a_renormalized_infinite_spatial_tail():
    assert (
        "not a covariantly renormalized full spatial response"
        in bounds.data()["UV_boundary"]
    )
    assert "not a physical cutoff" in bounds.data()["UV_boundary"]
    assert "two-propagator" in bounds.data()["preparation_and_norm"]
    assert "homogeneous subtraction" in bounds.data()["UV_boundary"]


@pytest.mark.parametrize(
    "name,value",
    audit.residuals().items(),
    ids=lambda x: x if isinstance(x, str) else None,
)
def test_all_exact_residuals(name, value):
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(s.cancel(x) == 0 for x in entries), name


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=lambda x: x if isinstance(x, str) else None
)
def test_all_scope_mutations_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_complete_counts_and_unclosed_frontier():
    assert len(audit.residuals()) == 25 and audit.scalar_entry_count() == 272
    assert len(audit.gates()) == 39 and all(audit.gates().values())
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 117
    assert len(audit.frontier()) == 9 and audit.matching()[-1] == audit.ITEM
