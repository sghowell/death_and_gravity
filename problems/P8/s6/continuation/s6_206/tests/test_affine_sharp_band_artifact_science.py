"""Independent all-polarization ultraviolet and original-shell checks."""

from functools import cache
from itertools import product

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_cd_metric_noise import stress
from p8_vacuum_affine_sharp_band_artifact import (
    angular,
    asymptotics,
    audit,
    geometry,
    polarizations,
)
from p8_vector_hadamard import series
from p8_vector_state import wkb


def cross(a, b):
    return mp.matrix(
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    )


def conjugate(v):
    return mp.matrix([mp.conj(x) for x in v])


def real_frame(n):
    index = min(range(3), key=lambda j: abs(n[j]))
    axis = mp.matrix([1 if j == index else 0 for j in range(3)])
    first = axis - n * (n.T * axis)[0]
    first /= mp.norm(first)
    return first, cross(n, first)


def joint_frame(k, p):
    magnitude = mp.norm(k)
    n0 = -k / magnitude
    ell = -k + p
    rho = mp.sqrt((ell.T * ell)[0])
    n = ell / rho
    skew = n * n0.T - n0 * n.T
    R = mp.eye(3) + skew + skew * skew / (1 + (n.T * n0)[0])
    first, second = real_frame(n0)
    return ell, rho, n, R, R * first, R * second


@cache
def coefficient_function(kind):
    expressions = []
    for order in range(1, 5):
        P = series.coefficient(kind, order)
        expressions.extend((P, s.diff(P, wkb.u), s.diff(P, wkb.z)))
    return s.lambdify((wkb.u, wkb.z), expressions, "mpmath", cse=True)


def W_data(kind, u, ell):
    a = (1 + u * u) ** 2
    H = 4 * u / (1 + u * u)
    omega = mp.sqrt(1000**2 + (ell.T * ell)[0] / a**2)
    z = 1 - 1000**2 / omega**2
    omega_prime = -H * z * omega
    z_prime = -2 * H * z * (1 - z)
    values = coefficient_function(kind)(u, z)
    W = omega
    derivative = omega_prime
    for order in range(1, 5):
        P, du, dz = values[3 * (order - 1) : 3 * order]
        W += P * omega ** (1 - 2 * order)
        derivative += (du + dz * z_prime) * omega ** (1 - 2 * order) + (
            1 - 2 * order
        ) * P * omega ** (-2 * order) * omega_prime
    return a, H, omega, z, W, derivative


def physical_modes(u, ell, n, e1, e2, phase_sign=-1):
    result = []
    for kind, e in (("transverse", e1), ("transverse", e2), ("longitudinal", n)):
        a, H, omega, z, W, derivative = W_data(kind, u, ell)
        f = 1 / mp.sqrt(2 * W)
        rate = H / 2 if kind == "transverse" else (mp.mpf("0.5") + z) * H
        momentum = (phase_sign * 1j * W - derivative / (2 * W) - rate) * f
        if kind == "transverse":
            E = momentum * e
            B = -phase_sign * 1j * cross(ell, e) * f / a
            V0 = mp.mpf(0)
            V = 1000 * f * e
        else:
            E = 1000 / omega * momentum * n
            B = mp.zeros(3, 1)
            V0 = phase_sign * 1j * mp.sqrt((ell.T * ell)[0]) / (a * omega) * momentum
            V = omega * f * n
        result.append((mp.matrix([*E, *B, V0, *V]), W))
    return result


@cache
def complete_stress_matrices():
    def real_number(x):
        n, d = s.fraction(x)
        return mp.mpf(int(n)) / int(d)

    return [
        mp.matrix([[real_number(M[i, j]) for j in range(10)] for i in range(10)])
        for M in stress.data()["quadratic_matrices"].values()
    ]


def contracted_spatial_matrix(tensor):
    matrices = complete_stress_matrices()
    return sum(
        (
            tensor[i, j] * matrices[4 * (i + 1) + j + 1]
            for i in range(3)
            for j in range(3)
        ),
        mp.zeros(10),
    )


def trace(M):
    return sum(M[j, j] for j in range(M.rows))


def tensor_fixtures():
    D = mp.matrix([[2, 1, -1], [1, -3, 2], [-1, 2, 1]]) / 7
    G = mp.matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 17
    return D, G


def leading_formula(D, G, n, a):
    return (
        4 * trace(D * G)
        - 4 * (n.T * (D * G + G * D) * n)[0]
        + 3 * (n.T * D * n)[0] * (n.T * G * n)[0]
    ) / (8 * a)


def physical_endpoint(K, u, n, P, D, G):
    k = K * n
    e1, e2 = real_frame(n)
    ell, _, nl, _, f1, f2 = joint_frame(k, P)
    left = physical_modes(u, k, n, e1, e2)
    right = physical_modes(u, ell, nl, f1, f2)
    MD, MG = contracted_spatial_matrix(D), contracted_spatial_matrix(G)
    terms = mp.matrix(3, 3)
    for i, (v, W) in enumerate(left):
        for j, (w, V) in enumerate(right):
            AD = (v.T * MD * w)[0]
            AG = (v.T * MG * w)[0]
            terms[i, j] = -mp.im(-1j * mp.conj(AD) * AG / (W + V))
    return terms


@pytest.mark.parametrize("K", ("1e6", "1e10", "1e20"))
@pytest.mark.parametrize("time", ("-0.5", "0", "0.41"))
@pytest.mark.parametrize("direction", ((3, 0, 4), (1, 2, 2), (0, 0, 1)))
@pytest.mark.parametrize("equal", (False, True))
def test_actual_full_W8_all_nine_pair_endpoint_symbol(K, time, direction, equal):
    with mp.workdps(110):
        K, u = mp.mpf(K), mp.mpf(time)
        n = mp.matrix(direction)
        n /= mp.norm(n)
        P = mp.matrix([11, -7, 5])
        D, G = tensor_fixtures()
        if equal:
            G = D
        terms = physical_endpoint(K, u, n, P, D, G)
        actual = sum(terms) / K
        expected = leading_formula(D, G, n, (1 + u * u) ** 2)
        scale = mp.norm(D) * mp.norm(G)
        assert (
            abs(actual - expected) < (mp.mpf("1e4") / K + mp.mpf("1e7") / K**2) * scale
        )
        mixed = sum(terms[i, j] for i, j in ((0, 2), (1, 2), (2, 0), (2, 1))) / K
        assert abs(mixed) < (mp.mpf("1e4") / K + mp.mpf("1e7") / K**2) * scale


@pytest.mark.parametrize("time", ("-0.5", "0", "0.5"))
def test_actual_longitudinal_pair_not_Maxwell_discarded(time):
    with mp.workdps(100):
        u = mp.mpf(time)
        K = mp.mpf("1e20")
        n = mp.matrix([0, 0, 1])
        D = mp.diag([-1, -1, 2]) / mp.sqrt(6)
        terms = physical_endpoint(K, u, n, mp.zeros(3, 1), D, D)
        expected = (n.T * D * n)[0] ** 2 / (8 * (1 + u * u) ** 2)
        assert expected > 0
        assert abs(terms[2, 2] / K - expected) < mp.mpf("1e-30")
        assert (
            abs(
                sum(terms) / K
                - sum(terms[i, j] for i in range(2) for j in range(2)) / K
            )
            > expected / 2
        )


@pytest.mark.parametrize("module", (polarizations, angular, geometry, asymptotics))
def test_complete_exact_core(module):
    packet = module.data()
    for value in packet["checks"].values():
        entries = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.simplify(x) == 0 for x in entries)
    assert all(v is True for v in packet["gates"].values())


MOMENT_ORDERS = tuple(
    (i, j, k) for i, j, k in product(range(5), repeat=3) if i + j + k <= 4
)


@pytest.mark.parametrize("orders", MOMENT_ORDERS)
def test_independent_full_weighted_angular_monomials(orders):
    with mp.workdps(80):
        i, j, k = orders
        azimuth = mp.quad(
            lambda phi: mp.cos(phi) ** i * mp.sin(phi) ** j,
            [0, mp.pi / 2, mp.pi, 3 * mp.pi / 2, 2 * mp.pi],
        )
        polar = mp.quad(
            lambda u: u * (-u) ** k * (1 - u * u) ** (mp.mpf(i + j) / 2), [0, 1]
        )
        exact = mp.mpf(str(s.N(angular.hemisphere_moment(i, j, k), 85)))
        assert abs(azimuth * polar - exact) < mp.mpf("1e-70")


def five_tensors():
    return [
        mp.diag([1, -1, 0]) / mp.sqrt(2),
        mp.matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) / mp.sqrt(2),
        mp.matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) / mp.sqrt(2),
        mp.matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) / mp.sqrt(2),
        mp.diag([-1, -1, 2]) / mp.sqrt(6),
    ]


@cache
def azimuth_average_coefficients(channel):
    with mp.workdps(110):
        D = five_tensors()[channel]

        def average(u):
            out = mp.mpf(0)
            for j in range(12):
                phi = 2 * mp.pi * j / 12
                n = mp.matrix(
                    [
                        mp.sqrt(1 - u * u) * mp.cos(phi),
                        mp.sqrt(1 - u * u) * mp.sin(phi),
                        u,
                    ]
                )
                out += leading_formula(D, D, n, mp.mpf(1))
            return 2 * mp.pi * out / 12

        values = [average(mp.mpf(j) / 8) for j in range(5)]
        matrix = mp.matrix(
            [[mp.mpf(j) ** r / 8**r for r in range(5)] for j in range(5)]
        )
        return tuple(mp.lu_solve(matrix, mp.matrix(values)))


def angular_average(channel, u):
    return sum(c * u**j for j, c in enumerate(azimuth_average_coefficients(channel)))


@cache
def normalized_lost_shell(channel, epsilon):
    with mp.workdps(100):
        e = mp.mpf(epsilon)

        def integrand(u):
            r = e * u + mp.sqrt(1 - e * e * (1 - u * u))
            return angular_average(channel, u) * (1 - r**4) / (4 * e * (2 * mp.pi) ** 3)

        return mp.quad(integrand, [-1, 0, e / 2])


@pytest.mark.parametrize("channel", range(5))
@pytest.mark.parametrize("epsilon", ("0.5", "0.1", "0.01", "0.0001", "1e-12"))
def test_original_full_shell_converges_to_actual_five_channel_coefficient(
    channel, epsilon
):
    with mp.workdps(90):
        e = mp.mpf(epsilon)
        actual = normalized_lost_shell(channel, epsilon)
        eigen = (18, 18, 12, 12, mp.mpf(28) / 3)[channel]
        expected = eigen / (512 * mp.pi**2)
        assert actual > 0
        assert abs(actual - expected) < 10 * e * expected
        if e < mp.mpf("1e-6"):
            assert abs(actual / expected - 1) < 10 * e


@pytest.mark.parametrize("epsilon", ("0.9", "0.5", "0.01", "1e-12"))
@pytest.mark.parametrize("u", ("-1", "-0.2", "0", "0.25", "1"))
def test_exact_original_boundary_and_grazing_strip(epsilon, u):
    with mp.workdps(90):
        K = mp.mpf(1)
        P = mp.mpf(epsilon)
        u = mp.mpf(u)
        r = P * u + mp.sqrt(K * K - P * P * (1 - u * u))
        assert abs(r * r - 2 * P * u * r + P * P - K * K) < mp.mpf("1e-80")
        assert (
            (r < K) == (u < P / (2 * K))
            if u != P / (2 * K)
            else abs(r - K) < mp.mpf("1e-80")
        )
        assert K - P <= r <= K + P
        grazing = P * (P / (2 * K)) + mp.sqrt(K * K - P * P * (1 - (P / (2 * K)) ** 2))
        assert abs(grazing - K) < mp.mpf("1e-80")


@pytest.mark.parametrize("axis", ((0, 0, 1), (3, 0, 4), (1, 2, 2)))
def test_full_tensor_coefficient_rotation_and_positive_channels(axis):
    with mp.workdps(80):
        n = mp.matrix(axis)
        n /= mp.norm(n)
        e1, e2 = real_frame(n)
        R = mp.matrix([[e1[i], e2[i], n[i]] for i in range(3)])
        tensors = [R * D * R.T for D in five_tensors()]
        for i, D in enumerate(tensors):
            for j, G in enumerate(tensors):
                actual = (
                    18 * trace(D * G)
                    - 12 * ((D * n).T * (G * n))[0]
                    - (n.T * D * n)[0] * (n.T * G * n)[0]
                )
                expected = (
                    (mp.mpf(18), mp.mpf(18), mp.mpf(12), mp.mpf(12), mp.mpf(28) / 3)[i]
                    if i == j
                    else 0
                )
                assert abs(actual - expected) < mp.mpf("1e-70")


def test_nonpolynomial_cusp_and_missing_second_leg_control():
    with mp.workdps(80):
        A = mp.mpf(18) / (512 * mp.pi**2)

        def conversion(p):
            return -abs(p) * A

        e = mp.mpf("1e-20")
        right = (conversion(e) - conversion(0)) / e
        left = (conversion(-e) - conversion(0)) / (-e)
        assert abs((right - left) + 2 * A) < mp.mpf("1e-70")
        assert right != left and conversion(0) == 0
        one_leg_minus_itself = mp.mpf(0)
        assert one_leg_minus_itself != conversion(e)


BAD_MOMENT = (True, False, 1.0, s.Float(1), "1", None, s.oo, s.I, s.Symbol("j"), -1, 5)


@pytest.mark.parametrize("value", BAD_MOMENT)
@pytest.mark.parametrize("position", range(3))
def test_exact_moment_guards_after_cache_warmup(value, position):
    angular.hemisphere_moment(0, 0, 0)
    angular.hemisphere_moment(1, 0, 0)
    orders = [0, 0, 0]
    orders[position] = value
    with pytest.raises((TypeError, ValueError)):
        angular.hemisphere_moment(*orders)


@pytest.mark.parametrize("orders", ((4, 1, 0), (2, 2, 1), (0, 4, 1), (3, 3, 3)))
def test_total_degree_guard(orders):
    with pytest.raises(ValueError):
        angular.hemisphere_moment(*orders)


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_all_named_scoped_residuals(name):
    value = audit.residuals()[name]
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(x == 0 for x in entries)


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_all_unsupported_scope_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_exact_scientific_counts():
    assert len(audit.residuals()) == 34 and audit.scalar_entry_count() == 272
    assert len(audit.gates()) == 38 and len(audit.controls()) == 9
    assert audit.rejected_inputs() == 128


def test_unchanged_open_frontier():
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert all(v is True for v in audit.gates().values())
