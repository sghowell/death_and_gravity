"""Independent full-matrix, actual-state and infinite-tail science fixtures."""

from functools import cache
from math import comb

import mpmath as mp
import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_matrix_adiabatic import audit, frame, initial, jets, riccati
from p8_vector_hadamard import cutoffs, preparation, series
from p8_vector_state import wkb
from scipy.linalg import expm


@pytest.mark.parametrize(
    "name,value",
    audit.residuals().items(),
    ids=lambda x: x if isinstance(x, str) else None,
)
def test_exact_identity(name, value):
    if isinstance(value, s.MatrixBase):
        assert all(s.cancel(x) == 0 for x in value), name
    else:
        assert s.cancel(value) == 0, name


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=lambda x: x if isinstance(x, str) else None
)
def test_unsupported_scope(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("name", tuple(audit.packets()))
def test_packet_gates(name):
    assert all(audit.packets()[name]["gates"].values())


def test_counts_and_unclosed_frontier():
    assert len(audit.residuals()) == 80
    assert audit.scalar_entry_count() == 257
    assert len(audit.gates()) == 53
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 112
    assert len(audit.frontier()) == 9
    assert audit.matching()[-1] == audit.ITEM


@pytest.mark.parametrize("power", (1, 2, 3))
@pytest.mark.parametrize("order", range(13))
def test_independent_literal_scale_derivatives(power, order):
    actual = s.diff((1 + jets.t**2) ** (-2 * power), jets.t, order) * (
        1 + jets.t**2
    ) ** (2 * power + order)
    assert s.cancel(actual - jets.scale_polynomial(power, order)) == 0


@pytest.mark.parametrize(
    "magnitude,time",
    (
        (0, s.Rational(-1, 4)),
        (1000, 0),
        (10**6, s.Rational(1, 5)),
        (10**20, s.Rational(1, 4)),
        (10**50, s.Rational(-1, 5)),
    ),
)
def test_actual_ill_conditioned_Hamiltonian_frame(magnitude, time):
    with mp.workdps(160):
        t = mp.mpf(str(s.N(time, 160)))
        k = mp.mpf(magnitude) * mp.matrix([mp.mpf(2) / 3, mp.mpf(1) / 3, mp.mpf(2) / 3])
        G0 = mp.matrix([[2, 1, -1], [1, -1, 2], [-1, 2, -1]]) / 10000
        G1 = mp.matrix([[1, -2, 1], [-2, 1, 1], [1, 1, -2]]) / 10000
        G2 = mp.matrix([[1, 1, 0], [1, -2, 1], [0, 1, 1]]) / 10000

        def gamma(x):
            return G0 + (x - t) * G1 + (x - t) ** 2 * G2 / 2

        def K(x):
            a = (1 + x * x) ** 2
            return mp.expm(gamma(x)) / a + (k * k.T) / (a**3 * 1000**2)

        K0 = K(t)
        Kd = mp.diff(K, t)
        vals, O = mp.eigsy(K0)
        roots = [mp.sqrt(x) for x in vals]
        B = O * mp.diag(roots) * O.T
        Q = B**-1 * Kd * B**-1
        rotated = O.T * Kd * O
        Bd = (
            O
            * mp.matrix(
                [
                    [rotated[i, j] / (roots[i] + roots[j]) for j in range(3)]
                    for i in range(3)
                ]
            )
            * O.T
        )
        L = B**-1 * Bd
        assert mp.norm(B * Bd + Bd * B - Kd) / max(1, mp.norm(Kd)) < mp.mpf("1e-50")
        assert mp.norm(L) <= mp.norm(Q) * (1 + mp.mpf("1e-40"))
        assert mp.norm(L) < 10
        a = (1 + t * t) ** 2
        E = mp.expm(gamma(t))
        cross = mp.matrix([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
        V = a * 1000**2 * mp.expm(-gamma(t)) + cross.T * E * cross / a
        omega2 = 1000**2 + (k.T * mp.expm(-gamma(t)) * k)[0] / a**2
        assert mp.norm(K0 * V - omega2 * mp.eye(3)) / omega2 < mp.mpf("1e-50")
        real_qp = mp.matrix(6)
        physical = mp.matrix(6)
        transform = mp.matrix(6)
        transform_dot = mp.matrix(6)
        Binv = B**-1
        for i in range(3):
            for j in range(3):
                physical[i, j + 3] = K0[i, j]
                physical[i + 3, j] = -V[i, j]
                transform[i, j] = Binv[i, j]
                transform[i + 3, j + 3] = B[i, j]
                transform_dot[i, j] = -(Binv * Bd * Binv)[i, j]
                transform_dot[i + 3, j + 3] = Bd[i, j]
                real_qp[i, j] = -L[i, j]
                real_qp[i + 3, j + 3] = L[j, i]
                real_qp[i, j + 3] = int(i == j)
                real_qp[i + 3, j] = -omega2 * int(i == j)
        actual = (transform_dot + transform * physical) * transform**-1
        assert mp.norm(actual - real_qp) / max(1, mp.norm(real_qp)) < mp.mpf("1e-50")


@pytest.mark.parametrize("condition", (1, 10**5, 10**30))
@pytest.mark.parametrize("order", (1, 2, 3, 4))
def test_independent_higher_square_root_balance(condition, order):
    t = s.Symbol("local_time", real=True)
    B0 = s.diag(1, 3, condition)
    G = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 31
    D = s.Matrix([[0, 1, 2], [1, -1, -1], [2, -1, 1]]) / 37
    B = B0 + t * G + t * t * D / 2
    derivatives = [s.diff(B, t, n).subs(t, 0) for n in range(order + 1)]
    M = [B0.inv() * d for d in derivatives]
    Q = B0.inv() * s.diff(B * B, t, order).subs(t, 0) * B0.inv()
    right = Q - sum(
        (comb(order, j) * M[j] * M[order - j].T for j in range(1, order)), s.zeros(3)
    )
    assert M[order] + M[order].T == right
    exact = s.Matrix(3, 3, lambda i, j: B0[j, j] * right[i, j] / (B0[i, i] + B0[j, j]))
    assert exact == M[order]
    assert sum(abs(x) ** 2 for x in M[order]) <= sum(abs(x) ** 2 for x in right)


def _fixture(seed=0):
    rng = np.random.default_rng(seed)
    G = rng.normal(size=(3, 3))
    D = rng.normal(size=(3, 3))
    R = (G - G.T) / 100
    S = (D + D.T) / 100
    r = (G + G.T + 1j * (D + D.T)) / 500
    return R, S, r


def _cov(r):
    F = np.vstack((np.eye(3) + r, -1j * (np.eye(3) - r))) / np.sqrt(2)
    return np.real(F @ np.linalg.inv(np.eye(3) - r.conj().T @ r) @ F.conj().T)


@pytest.mark.parametrize("seed", range(8))
def test_full_noncommuting_symplectic_covariance(seed):
    R, S, r0 = _fixture(seed)
    omega = 1000 + 137 * seed
    unit = np.eye(3)
    complexflow = np.block([[-1j * omega * unit + R, S], [S, 1j * omega * unit + R]])
    T = np.block([[unit, 1j * unit], [unit, -1j * unit]]) / np.sqrt(2)
    realflow = np.real(T.conj().T @ complexflow @ T)
    dt = 0.001 * (seed + 1)
    Z = expm(complexflow * dt) @ np.vstack((unit, r0))
    r = Z[3:] @ np.linalg.inv(Z[:3])
    realS = expm(realflow * dt)
    assert np.linalg.norm(r - r.T) < 1e-13
    assert np.linalg.norm(_cov(r) - realS @ _cov(r0) @ realS.T) < 2e-12
    assert np.linalg.norm(R @ S - S @ R) > 1e-7


@pytest.mark.parametrize("seed", range(8))
def test_full_covariance_uniform_Lipschitz(seed):
    _, _, r = _fixture(seed)
    _, _, q = _fixture(seed + 31)
    assert max(np.linalg.norm(r, 2), np.linalg.norm(q, 2)) < 0.1
    actual = np.linalg.norm(_cov(r) - _cov(q), 2)
    assert actual < 16 * np.linalg.norm(r - q, 2)
    exact = frame.covariance(s.Matrix(r))
    assert np.linalg.norm(np.asarray(exact, dtype=float) - _cov(r)) < 1e-13


@cache
def _finite_reference(order):
    R = s.Matrix([[0, 1, -2], [-1, 0, 3], [2, -3, 0]]) / 37
    G = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 41
    D = s.Matrix([[2, -1, 1], [-1, 0, 2], [1, 2, -2]]) / 43
    omega = [s.Integer(1000)] + [
        s.Rational((-1) ** j, j + 1) for j in range(1, order + 1)
    ]
    rotations = [R * s.Rational((-1) ** j, j + 1) for j in range(order + 1)]
    squeezes = [
        (G if j % 2 == 0 else D) * s.Rational((-1) ** j, j + 1)
        for j in range(order + 1)
    ]
    return (
        omega,
        rotations,
        squeezes,
        riccati.reference_jets(omega, rotations, squeezes, order),
    )


@pytest.mark.parametrize("order", (1, 2, 3, 4, 6, 8, 10))
def test_independent_full_ordered_Riccati_residual(order):
    omega, R, S, ref = _finite_reference(order)
    total = sum((ref[n][0] for n in ref), s.zeros(3))
    derivative = sum((ref[n][1] for n in ref), s.zeros(3))
    literal = (
        derivative
        - 2 * s.I * omega[0] * total
        - R[0] * total
        + total * R[0]
        - S[0]
        + total * S[0] * total
    )
    tail = ref[order][1] - R[0] * ref[order][0] + ref[order][0] * R[0]
    tail += sum(
        (ref[j][0] * S[0] * ref[l][0] for j in ref for l in ref if j + l >= order),
        s.zeros(3),
    )
    assert all(s.cancel(x) == 0 for x in literal - tail)
    assert all(s.cancel(x) == 0 for x in total - total.T)
    assert S[0] * S[1] != S[1] * S[0]


@cache
def _doubled_reference():
    omega, R, S, _ref = _finite_reference(10)
    return riccati.reference_jets([2 * x for x in omega], R, S, 10)


@pytest.mark.parametrize("coefficient", range(1, 11))
def test_frequency_order_without_commuting_matrices(coefficient):
    ref = _finite_reference(10)[3]
    doubled = _doubled_reference()
    assert all(
        s.cancel(x) == 0
        for x in doubled[coefficient][0] * 2**coefficient - ref[coefficient][0]
    )


@cache
def _base_polynomials(kind):
    H = wkb.background()["H"]
    lam = wkb.background()["lambda"]
    d = H / 2 if kind == "transverse" else H * (s.Rational(1, 2) + wkb.z)
    frequency = [s.Integer(1)]
    squeeze = [d + lam / 2]
    for _ in range(10):
        frequency.append(s.factor(series.D0(frequency[-1]) + lam * frequency[-1]))
        squeeze.append(series.D0(squeeze[-1]))
    return frequency, squeeze


def _truncated_product(a, b, N):
    return [sum(a[j] * b[n - j] for j in range(n + 1)) for n in range(N + 1)]


def _truncated_inverse(a, N):
    out = [1 / a[0]]
    for n in range(1, N + 1):
        out.append(-out[0] * sum(a[j] * out[n - j] for j in range(1, n + 1)))
    return out


@cache
def _actual_initial_series(kind, nu):
    point = {wkb.u: -s.Rational(1, 2), wkb.z: 1 - s.Integer(1000) ** 2 / nu**2}
    frequency, squeeze = _base_polynomials(kind)
    omega = [nu * x.subs(point) for x in frequency]
    S = [s.Matrix([[x.subs(point)]]) for x in squeeze]
    ref = riccati.reference_jets(omega, [s.zeros(1)] * 11, S, 10)
    N = 10
    lam = wkb.background()["lambda"]
    H = wkb.background()["H"]
    d = H / 2 if kind == "transverse" else H * (s.Rational(1, 2) + wkb.z)
    w = [s.Integer(0)] * (N + 1)
    wd = [s.Integer(0)] * (N + 1)
    w[0] = s.Integer(1)
    for n in range(1, 6):
        p = series.coefficient(kind, n)
        w[2 * n] = p.subs(point)
        wd[2 * n] = (series.D0(p) - 2 * n * lam * p).subs(point)
    ratio = _truncated_product(wd, _truncated_inverse(w, N), N)
    beta = [s.Integer(0)] + [ratio[n - 1] / 2 for n in range(1, N + 1)]
    beta[1] += (d + lam / 2).subs(point)
    denominator = [w[n] - s.I * beta[n] for n in range(N + 1)]
    denominator[0] += 1
    inv = _truncated_inverse(denominator, N)
    graph = [2 * inv[n] for n in range(N + 1)]
    graph[0] -= 1
    return ref, graph, point


@pytest.mark.parametrize("kind", ("transverse", "longitudinal"))
@pytest.mark.parametrize("nu", (s.Integer(1000), s.Integer(10) ** 16))
@pytest.mark.parametrize("order", range(1, 11))
def test_actual_W10_and_full_Riccati_coefficients(kind, nu, order):
    ref, graph, _ = _actual_initial_series(kind, nu)
    assert s.cancel(ref[order][0][0] - graph[order] / nu**order) == 0


@pytest.mark.parametrize("kind", ("transverse", "longitudinal"))
def test_actual_high_band_all_order_state_not_reset(kind):
    nu = s.Integer(10) ** 16
    actual = preparation.initial_data(kind, 1000, nu)
    assert actual["active_higher_derivative_orders"][:3] == [6, 8, 10]
    ref, _, point = _actual_initial_series(kind, nu)
    H = wkb.background()["H"]
    d = (H / 2 if kind == "transverse" else H * (s.Rational(1, 2) + wkb.z)).subs(point)
    W = actual["all_order_frequency"]
    v = actual["all_order_half_log_rate"] + d
    r = (nu - W + s.I * v) / (nu + W - s.I * v)
    rhat = sum(ref[n][0][0] for n in ref)
    difference = s.cancel(r - rhat)
    assert (
        s.cancel(s.re(difference) ** 2 + s.im(difference) ** 2)
        < (s.Integer(10) ** 33 / nu**11) ** 2
    )
    assert actual["all_order_frequency"] != actual["frozen_frequency"]


@pytest.mark.parametrize("kind", ("transverse", "longitudinal"))
def test_actual_partially_active_first_cutoff(kind):
    actual = preparation.initial_data(kind, 1000, 3000)
    assert actual["active_higher_derivative_orders"] == [6]
    assert actual["included_correction_terms"][6]["weight"] == s.Rational(1, 2)
    assert actual["all_order_frequency"] != actual["frozen_frequency"]


@pytest.mark.parametrize("ratio", (s.Rational(1, 2), 1, s.Rational(3, 2), 2, 3))
def test_higher_cutoff_weight_is_bounded_without_dropping_term(ratio):
    value = cutoffs.turn_on(s.sympify(ratio))
    assert 0 <= value <= 1
    if ratio == s.Rational(3, 2):
        assert value != 0 and value != 1


@pytest.mark.parametrize("lower", (10**16, 2 * 10**16, 10**20))
def test_independent_infinite_energy_tail_integral(lower):
    with mp.workdps(60):
        b = mp.mpf(99) / 100 / (mp.mpf(25) / 16) ** 2
        K = mp.mpf(lower)
        transformed = mp.quad(lambda x: x**5, [0, 1])
        exact_majorant = (
            3
            * mp.mpf(2)
            * 10**31
            / (2 * mp.pi**2 * b ** mp.mpf("1.5"))
            * K**-6
            * transformed
        )
        assert abs(transformed - mp.mpf(1) / 6) < mp.mpf("1e-55")
        assert exact_majorant < mp.mpf("1e-65")
        assert (
            initial.data()["integrable_energy_weighted_covariance_tail_upper"]
            < s.Rational(1, 10) ** 65
        )


@pytest.mark.parametrize("order", (0, -1, 11, True, False, 1.0, s.Integer(1)))
def test_reference_rejects_unsupported_finite_order(order):
    with pytest.raises(ValueError):
        riccati.reference_jets(
            [s.Integer(1000)] * 12, [s.zeros(3)] * 12, [s.eye(3)] * 12, order
        )


def test_reference_rejects_missing_derivative_jet():
    with pytest.raises(ValueError):
        riccati.reference_jets([s.Integer(1000)], [s.zeros(3)], [s.eye(3)], 10)


def test_omitted_time_dependent_frame_is_a_negative_control():
    t = s.Symbol("time", real=True)
    B = s.eye(3) / (1 + t * t)
    L = (B.inv() * s.diff(B, t)).subs(t, s.Rational(1, 4))
    assert L == -s.Rational(8, 17) * s.eye(3)
    actual_squeeze = -(L + L.T) / 2
    incorrectly_static_frame = s.zeros(3)
    assert actual_squeeze - incorrectly_static_frame == s.Rational(8, 17) * s.eye(3)


def test_covariance_tail_does_not_claim_response_or_cutoff():
    text = str(initial.data()) + str(audit.observable())
    assert (
        "not a physical cutoff" in text
        or "not_parameter" in text
        or "not a response" in text
        or "not a response-parameter" in text
    )
    assert "subtraction" in text and "derivatives" in text
    assert audit.matching()[-1]["status"] == audit.ITEM["status"]
