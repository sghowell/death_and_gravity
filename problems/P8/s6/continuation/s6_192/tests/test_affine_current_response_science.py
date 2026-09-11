"""Independent actual metric, covariance flow and holomorphic comparison tests."""

from functools import cache
from math import comb

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_current_response import adiabatic, audit, low, tail, vertices
from p8_vacuum_affine_matrix_adiabatic import riccati


@pytest.mark.parametrize(
    "name,value",
    audit.residuals().items(),
    ids=lambda x: x if isinstance(x, str) else None,
)
def test_exact_identity(name, value):
    assert (
        all(s.cancel(x) == 0 for x in value)
        if isinstance(value, s.MatrixBase)
        else s.cancel(value) == 0
    ), name


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
    assert len(audit.residuals()) == 62 and audit.scalar_entry_count() == 641
    assert len(audit.gates()) == 41 and len(audit.controls()) == 9
    assert audit.rejected_inputs() == 114 and len(audit.frontier()) == 9
    assert audit.matching()[-1] == audit.ITEM


def block(A, B):
    result = mp.zeros(A.rows + B.rows, A.cols + B.cols)
    result[: A.rows, : A.cols] = A
    result[A.rows :, A.cols :] = B
    return result


def opnorm(M):
    values = mp.eigsy((M.T * M + (M.T * M).T) / 2, eigvals_only=True)
    return mp.sqrt(max(abs(x) for x in values))


def positive_root(M):
    values, U = mp.eigsy((M + M.T) / 2)
    return U * mp.diag([mp.sqrt(x) for x in values]) * U.T


def actual_metric(e, momentum):
    m = mp.mpf(1000)
    a = mp.mpf(25) / 16
    eigen = [mp.mpf(1) / 2, -mp.mpf(1) / 6, -mp.mpf(1) / 3]
    D = mp.matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 10
    k = mp.matrix([1, 2, 3]) * mp.mpf(momentum) / mp.sqrt(14)
    C = mp.matrix([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    E = mp.diag([mp.exp(e * x) for x in eigen])
    Q = mp.diag([mp.exp(-e * x) for x in eigen])

    def frechet(sign):
        result = mp.zeros(3)
        for i in range(3):
            for j in range(3):
                x = sign * e * (eigen[i] - eigen[j])
                divided = mp.expm1(x) / x if x else mp.mpf(1)
                result[i, j] = sign * D[i, j] * mp.exp(sign * e * eigen[j]) * divided
        return result

    ED, QD = frechet(1), frechet(-1)
    K = E / a + k * k.T / (a**3 * m**2)
    V = a * m**2 * Q + C.T * E * C / a
    M = block(V, K)
    MD = block(a * m**2 * QD + C.T * ED * C / a, ED / a)
    omega = mp.sqrt(m * m + (k.T * Q * k)[0] / a**2)
    B = positive_root(K)
    T = block(B / mp.sqrt(omega), mp.sqrt(omega) * B**-1)
    G = T.T * MD * T / omega
    return M, MD, omega, T, G


def to_sympy(M):
    return s.Matrix(
        [
            [s.Float(mp.nstr(M[i, j], 85), 85) for j in range(M.cols)]
            for i in range(M.rows)
        ]
    )


def to_mp(M):
    return mp.matrix(
        [[mp.mpf(str(s.N(M[i, j], 80))) for j in range(M.cols)] for i in range(M.rows)]
    )


@cache
def actual_fixture(momentum):
    with mp.workdps(90):
        e = mp.mpf(1) / 200
        result = actual_metric(e, momentum)
        frames = [
            mp.diff(lambda x: actual_metric(x, momentum)[3], e, n) for n in range(3)
        ]
        metric = [
            mp.diff(lambda x: actual_metric(x, momentum)[1], e, n) for n in range(3)
        ]
        omega = [
            mp.diff(lambda x: actual_metric(x, momentum)[2], e, n) for n in range(3)
        ]
        direct = [
            mp.diff(lambda x: actual_metric(x, momentum)[4], e, n) for n in range(3)
        ]
        rebuilt = vertices.vertex_jet(
            list(map(to_sympy, frames)),
            [s.Float(mp.nstr(x, 85), 85) for x in omega],
            list(map(to_sympy, metric)),
        )
        return result, frames, direct, list(map(to_mp, rebuilt))


@pytest.mark.parametrize("momentum", (0, 1000, 10**12))
@pytest.mark.parametrize("order", range(3))
def test_actual_anisotropic_metric_vertex_derivatives(momentum, order):
    with mp.workdps(70):
        _, _, direct, rebuilt = actual_fixture(momentum)
        assert mp.norm(direct[order] - rebuilt[order]) < mp.mpf("1e-40") * max(
            1, mp.norm(direct[order])
        )
        assert opnorm(direct[order]) < int(vertices.G[order])
        assert mp.norm(direct[order][:3, 3:]) < mp.mpf("1e-50")


@pytest.mark.parametrize("momentum", (0, 1000, 10**12))
@pytest.mark.parametrize("order", (1, 2))
def test_actual_noncommuting_positive_root_frame_jets(momentum, order):
    with mp.workdps(70):
        result, frames, _, _ = actual_fixture(momentum)
        M, _, omega, T, _ = result
        assert mp.norm(T.T * M * T - omega * mp.eye(6)) < mp.mpf("1e-40") * omega
        assert opnorm(frames[0] ** -1 * frames[order]) < int(vertices.FRAME[order])
        if momentum:
            assert mp.norm(frames[0] * frames[1] - frames[1] * frames[0]) > mp.mpf(
                "1e-30"
            )


@pytest.mark.parametrize("e", (-s.Rational(1, 100), 0, s.Rational(1, 100)))
@pytest.mark.parametrize("order", range(3))
def test_actual_energy_relative_higher_metric_vertex(e, order):
    with mp.workdps(70):
        ep = mp.mpf(str(s.N(e, 70)))
        M = actual_metric(ep, 1000)[0]
        H = mp.diff(lambda x: actual_metric(x, 1000)[1], ep, order)
        inv = positive_root(M) ** -1
        assert opnorm(inv * H * inv) < 2


def numeric_generator(w, R, S):
    n = R.rows
    A = mp.zeros(2 * n)
    A[:n, :n] = R + S
    A[:n, n:] = w * mp.eye(n)
    A[n:, :n] = -w * mp.eye(n)
    A[n:, n:] = R - S
    return A


@cache
def flow_fixture(frequency, length):
    with mp.workdps(75):
        T = mp.mpf(length) / 5
        R = mp.matrix([[0, 1, -2], [-1, 0, 1], [2, -1, 0]]) / 13
        S = mp.matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 17
        D = mp.matrix([[2, -1, 1], [-1, 0, 2], [1, 2, -2]]) / 19
        A = numeric_generator(mp.mpf(frequency), R, S)
        A1 = numeric_generator(mp.mpf(2), R / 7, D)
        A2 = numeric_generator(mp.mpf(-3), R / 11, S / 5)
        C = mp.diag([1, 2, 3, 4, 5, 6]) / 12
        augmented = mp.zeros(18)
        for j in range(3):
            augmented[6 * j : 6 * (j + 1), 6 * j : 6 * (j + 1)] = A
        augmented[:6, 6:12] = A1
        augmented[6:12, 12:] = A1
        augmented[:6, 12:] = A2 / 2
        E = mp.expm(T * augmented)
        U = [E[:6, :6], E[:6, 6:12], 2 * E[:6, 12:]]
        Sigma = [
            sum((comb(n, j) * U[j] * C * U[n - j].T for j in range(n + 1)), mp.zeros(6))
            for n in range(3)
        ]

        def direct(e):
            exp = mp.expm(T * (A + e * A1 + e**2 * A2 / 2))
            return exp * C * exp.T

        expected = [mp.diff(direct, mp.mpf(0), n) for n in range(3)]
        rate = opnorm(S)
        initial = opnorm(C)
        b1, b2 = opnorm(A1), opnorm(A2)
        bounds = [
            initial * mp.exp(2 * rate * T),
            initial * mp.exp(2 * rate * T) * 2 * b1 * T,
            initial * mp.exp(2 * rate * T) * (2 * b2 * T + 4 * b1**2 * T**2),
        ]
        return U, Sigma, expected, bounds, rate, T


@pytest.mark.parametrize("frequency", (1, 1000, 10**6))
@pytest.mark.parametrize("length", (1, 5))
@pytest.mark.parametrize("order", range(3))
def test_independent_full_noncommuting_covariance_Duhamel_response(
    frequency, length, order
):
    with mp.workdps(60):
        U, actual, expected, bounds, rate, T = flow_fixture(frequency, length)
        assert mp.norm(actual[order] - expected[order]) < mp.mpf("1e-40")
        assert opnorm(actual[order]) < bounds[order]
        assert opnorm(U[0]) < mp.exp(rate * T)


@pytest.mark.parametrize("frequency", (1, 1000, 10**8))
def test_pure_fast_flow_is_orthogonal_without_momentum_growth(frequency):
    with mp.workdps(70):
        R = mp.matrix([[0, 1, -2], [-1, 0, 1], [2, -1, 0]]) / 13
        A = numeric_generator(mp.mpf(frequency), R, mp.zeros(3))
        U = mp.expm(A)
        assert mp.norm(U.T * U - mp.eye(6)) < mp.mpf("1e-50")


def test_omitting_iterated_first_response_changes_second_covariance():
    e = s.Symbol("epsilon", real=True)
    U = s.diag(s.exp(e), s.exp(-e))
    C = U * U.T / 2
    assert C.diff(e, 2).subs(e, 0) == 2 * s.eye(2)
    assert s.zeros(2) != C.diff(e, 2).subs(e, 0)


def test_omitting_physical_vertex_derivative_changes_current():
    e = s.Symbol("epsilon", real=True)
    Sigma = s.eye(6)
    G = (1 + e + e**2) * s.eye(6)
    J = -s.trace(G * Sigma) / 2
    assert s.diff(J, e).subs(e, 0) == -3
    assert s.diff(J, e, 2).subs(e, 0) == -6


def adjoint(M):
    return mp.matrix([[mp.conj(M[j, i]) for j in range(M.rows)] for i in range(M.cols)])


def conjugate(M):
    return mp.matrix([[mp.conj(M[i, j]) for j in range(M.cols)] for i in range(M.rows)])


def marker_graph_coefficients(e):
    G = mp.matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 5
    D = mp.matrix([[2, -1, 1], [-1, 0, 2], [1, 2, -2]]) / 7
    F = mp.matrix([[0, 1, 2], [1, 1, -1], [2, -1, -1]]) / 11
    return [mp.zeros(3)] + [
        (mp.j**n)
        * (G / (n + 1) + e * D / (n + 2) + e**2 * F / (n + 3))
        / mp.mpf(20) ** n
        for n in range(1, 7)
    ]


def ordered_holomorphic_covariance(z, e):
    graph = marker_graph_coefficients(e)

    def W(point):
        r = sum((point**n * graph[n] for n in range(1, 7)), mp.zeros(3))
        sharp = sum((point**n * adjoint(graph[n]) for n in range(1, 7)), mp.zeros(3))
        F = mp.zeros(6, 3)
        F[:3, :] = (mp.eye(3) + r) / mp.sqrt(2)
        F[3:, :] = -mp.j * (mp.eye(3) - r) / mp.sqrt(2)
        Fsharp = mp.zeros(3, 6)
        Fsharp[:, :3] = (mp.eye(3) + sharp) / mp.sqrt(2)
        Fsharp[:, 3:] = mp.j * (mp.eye(3) - sharp) / mp.sqrt(2)
        return F * (mp.eye(3) - sharp * r) ** -1 * Fsharp

    return (W(z) + conjugate(W(mp.conj(z)))) / 2


@cache
def marker_fixture():
    e = s.Symbol("epsilon", real=True)
    G = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 5
    D = s.Matrix([[2, -1, 1], [-1, 0, 2], [1, 2, -2]]) / 7
    F = s.Matrix([[0, 1, 2], [1, 1, -1], [2, -1, -1]]) / 11
    graph = [s.zeros(3)] + [
        s.I**n
        * (G / (n + 1) + e * D / (n + 2) + e**2 * F / (n + 3))
        / s.Integer(20) ** n
        for n in range(1, 7)
    ]
    coeff = adiabatic.covariance_coefficients(graph, 6)
    with mp.workdps(70):
        ep = mp.mpf(1) / 200
        literal = {
            (n, a): to_mp(coeff[n].diff(e, a).subs(e, s.Rational(1, 200)))
            for n in range(7)
            for a in range(3)
        }
        fourier = {(n, a): mp.zeros(6) for n in range(7) for a in range(3)}
        count = 32
        radius = mp.mpf(1) / 5
        for j in range(count):
            z = radius * mp.exp(2 * mp.pi * mp.j * j / count)
            samples = [
                mp.diff(lambda x, z=z: ordered_holomorphic_covariance(z, x), ep, a)
                for a in range(3)
            ]
            for n in range(7):
                for a in range(3):
                    fourier[n, a] += samples[a] / z**n / count
        return literal, fourier


@pytest.mark.parametrize("order", range(7))
@pytest.mark.parametrize("parameter", range(3))
def test_independent_complex_Cauchy_complete_covariance_coefficients(order, parameter):
    with mp.workdps(60):
        literal, fourier = marker_fixture()
        assert mp.norm(literal[order, parameter] - fourier[order, parameter]) < mp.mpf(
            "1e-48"
        )


@pytest.mark.parametrize("order", (1, 3, 5))
@pytest.mark.parametrize("parameter", range(3))
def test_physical_current_even_but_full_covariance_odd(order, parameter):
    with mp.workdps(60):
        _, fourier = marker_fixture()
        value = fourier[order, parameter]
        assert mp.norm(value[:3, :3]) + mp.norm(value[3:, 3:]) < mp.mpf("1e-48")
        if parameter == 0:
            assert mp.norm(value[:3, 3:]) > mp.mpf("1e-20")


def test_conjugating_complex_marker_changes_the_analytic_reference():
    with mp.workdps(60):
        z = (1 + mp.j) / 10
        e = mp.mpf(1) / 200
        graph = marker_graph_coefficients(e)
        r = sum((z**n * graph[n] for n in range(1, 7)), mp.zeros(3))
        F = mp.zeros(6, 3)
        F[:3, :] = (mp.eye(3) + r) / mp.sqrt(2)
        F[3:, :] = -mp.j * (mp.eye(3) - r) / mp.sqrt(2)
        wrong = F * (mp.eye(3) - adjoint(r) * r) ** -1 * adjoint(F)
        wrong = mp.matrix([[mp.re(wrong[i, j]) for j in range(6)] for i in range(6)])
        assert mp.norm(wrong - ordered_holomorphic_covariance(z, e)) > mp.mpf("1e-6")


@pytest.mark.parametrize("frequency", (1000, 10**6, 10**16))
def test_independent_finite_reference_even_current_remainder(frequency):
    order = 10
    R = s.Matrix([[0, 1, -2], [-1, 0, 1], [2, -1, 0]]) / 7
    S = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 11
    omega = [s.Integer(frequency)] + [
        s.Rational(frequency, j + 1) for j in range(1, order + 1)
    ]
    b = riccati.reference_jets(
        omega,
        [R / (j + 1) for j in range(order + 1)],
        [S / (j + 1) for j in range(order + 1)],
        order,
    )
    graph = [s.zeros(3)] + [b[n][0] for n in range(1, order + 1)]
    coeff = adiabatic.covariance_coefficients(graph, 4)
    with mp.workdps(160):

        def complex_mp(M):
            return mp.matrix(
                [
                    [
                        mp.mpc(str(s.N(s.re(x), 155)), str(s.N(s.im(x), 155)))
                        for x in row
                    ]
                    for row in M.tolist()
                ]
            )

        r = complex_mp(sum(graph, s.zeros(3)))
        F = mp.zeros(6, 3)
        F[:3, :] = (mp.eye(3) + r) / mp.sqrt(2)
        F[3:, :] = -mp.j * (mp.eye(3) - r) / mp.sqrt(2)
        sigma = F * (mp.eye(3) - adjoint(r) * r) ** -1 * adjoint(F)
        sigma = mp.matrix([[mp.re(sigma[i, j]) for j in range(6)] for i in range(6)])
        reference = sum((complex_mp(coeff[n]) for n in (0, 2, 4)), mp.zeros(6))
        G = mp.diag([1, -1, 0, 2, -2, 0]) / 2
        current = (
            -mp.mpf(frequency)
            * sum((G * (sigma - reference))[i, i] for i in range(6))
            / 2
        )
        bound = 24 * mp.mpf(1000) ** 6 / mp.mpf(frequency) ** 5
        assert 0 < abs(current) < bound


@pytest.mark.parametrize("parameter", range(3))
def test_independent_complete_actual_reference_current_radial_tail(parameter):
    from p8_vacuum_affine_matrix_response_tail import covariance

    with mp.workdps(70):
        K = mp.mpf(10) ** 16
        m = mp.mpf(1000)
        b = mp.mpf(99) / 100 / (mp.mpf(25) / 16) ** 2
        weighted = covariance.constants()[
            "actual_weighted_parameter_error_coefficients"
        ]
        value = mp.mpf(0)
        for j in range(parameter + 1):
            a = parameter - j
            power = 9 - a
            coefficient = mp.mpf(str(s.N(weighted[a], 70)))
            radial = mp.quad(
                lambda x, power=power: (
                    x * mp.sqrt(x * x - (m / K) ** 2) * x ** (-power)
                ),
                [1, mp.inf],
            )
            value += (
                3
                * comb(parameter, j)
                * int(vertices.G[j])
                * coefficient
                * K ** (3 - power)
                * radial
                / (2 * mp.pi**2 * b ** mp.mpf("1.5"))
            )
        assert value < mp.mpf(str(s.N(tail.STATE_TAIL[parameter], 70)))


@pytest.mark.parametrize("parameter", range(3))
def test_independent_reference_comparison_infinite_radial_tail(parameter):
    with mp.workdps(70):
        K = mp.mpf(10) ** 16
        m = mp.mpf(1000)
        b = mp.mpf(99) / 100 / (mp.mpf(25) / 16) ** 2
        radial = mp.quad(
            lambda x: x * mp.sqrt(x * x - (m / K) ** 2) * x**-5, [1, mp.inf]
        )
        value = (
            2
            * int(adiabatic.CONTOUR_CURRENT[parameter])
            * m**6
            / K**2
            * radial
            / (2 * mp.pi**2 * b ** mp.mpf("1.5"))
        )
        assert value < mp.mpf(str(s.N(tail.REFERENCE_TAIL[parameter], 70)))


@pytest.mark.parametrize("parameter", range(3))
def test_independent_complete_finite_band_comparison_integral(parameter):
    with mp.workdps(70):
        K = mp.mpf(10) ** 16
        m = mp.mpf(1000)
        b = mp.mpf(99) / 100 / (mp.mpf(25) / 16) ** 2
        ratio = m / K
        first = mp.quad(
            lambda x: x ** (2 + parameter) * mp.sqrt(x * x - ratio**2), [ratio, 1]
        )
        zeroth = mp.quad(lambda x: x**2 * mp.sqrt(x * x - ratio**2), [ratio, 1])
        value = (
            int(low.CURRENT[parameter]) * K ** (4 + parameter) * first
            + 9 * int(adiabatic.CONTOUR_CURRENT[parameter]) * K**4 * zeroth
        ) / (2 * mp.pi**2 * b ** mp.mpf("1.5"))
        assert value < int(tail.COMPARISON[parameter])


@pytest.mark.parametrize(
    "epsilon",
    (
        -s.Rational(1, 100),
        -s.Rational(1, 10**8),
        0,
        s.Rational(1, 10**8),
        s.Rational(1, 100),
    ),
)
def test_finite_amplitude_Taylor_including_zero_and_negative_direction(epsilon):
    with mp.workdps(70):
        e = mp.mpf(str(s.N(epsilon, 70)))
        remainder = mp.expm1(e) - e
        bound = e * e * mp.exp(mp.mpf(1) / 100) / 2
        assert abs(remainder) <= bound
        if epsilon:
            assert abs(remainder) < bound
        else:
            assert remainder == bound == 0


@pytest.mark.parametrize("order", (True, False, -1, 11, s.Integer(2), 2.0, "2", None))
def test_non_native_or_unsupported_covariance_marker_order(order):
    with pytest.raises(ValueError):
        adiabatic.covariance_coefficients([s.zeros(3)] * 11, order)


@pytest.mark.parametrize("graph", ([], [s.eye(3)], [s.zeros(3)]))
def test_missing_or_nonzero_constant_graph_rejected(graph):
    with pytest.raises(ValueError):
        adiabatic.covariance_coefficients(graph, 2)


@pytest.mark.parametrize("which", range(3))
def test_incomplete_physical_vertex_parameter_jets_rejected(which):
    args = [[s.eye(3)] * 3, [s.Integer(1)] * 3, [s.eye(3)] * 3]
    args[which] = args[which][:-1]
    with pytest.raises(ValueError):
        vertices.vertex_jet(*args)


def test_no_fixed_covariant_matching_or_full_feedback_claim():
    observable = audit.observable()
    assert "not yet matched" in observable["boundary"]
    assert "mode comparison" in adiabatic.data()["comparison_definition"]
    assert "not yet" in tail.data()["complete_limit"]
    assert all(row["status"] != "COMPLETE" for row in audit.frontier())
