"""Independent mixed-frame, ordered-reference and C2 covariance-tail fixtures."""

from functools import cache
from math import comb

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_matrix_adiabatic import riccati as old_reference
from p8_vacuum_affine_matrix_response_tail import (
    audit,
    covariance,
    mixed,
    reference,
    variation,
)


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
def test_all_packet_gates(name):
    assert all(audit.packets()[name]["gates"].values())


def test_exact_counts_and_unclosed_frontier():
    assert len(audit.residuals()) == 225
    assert audit.scalar_entry_count() == 370
    assert len(audit.gates()) == 39
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 113
    assert len(audit.frontier()) == 9
    assert audit.matching()[-1] == audit.ITEM


@pytest.mark.parametrize("j", range(7))
@pytest.mark.parametrize("a", range(3))
def test_independent_mixed_exponential_generating_function(j, a):
    t = s.Symbol("t", real=True)
    scalar = s.exp(s.Rational(1, 100) * (s.exp(t) - 1) + a * t)
    expected = s.diff(scalar, t, j).subs(t, 0)
    actual = mixed.constants()["exp"][j, a]
    if (j, a) == (0, 0):
        assert actual == expected
    else:
        assert actual == 2 * expected


@pytest.mark.parametrize(
    "j,a", ((0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2), (3, 1))
)
def test_actual_noncommuting_exponential_mixed_jets(j, a):
    with mp.workdps(60):
        G = mp.matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 20
        D = mp.matrix([[2, -1, 1], [-1, 0, 2], [1, 2, -2]]) / 20
        B = mp.matrix([[0, 1, 2], [1, 1, -1], [2, -1, -1]]) / 30

        def direction(t):
            return G + t * D + t * t * B / 2

        epsilon = mp.mpf(1) / 200

        def exponent(t, e):
            return mp.expm(e * direction(t))

        derivative = mp.diff(exponent, (mp.mpf(0), epsilon), (j, a))
        balance = mp.expm(-epsilon * G / 2)
        actual = balance * derivative * balance
        actual = (actual + actual.T) / 2
        norm = max(abs(x) for x in mp.eigsy(actual, eigvals_only=True))
        bound = mp.mpf(str(s.N(mixed.constants()["exp"][j, a], 60)))
        assert norm < bound
        assert mp.norm(G * D - D * G) > 0


@pytest.mark.parametrize(
    "j,a", ((0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2))
)
def test_independent_complete_mixed_square_root_product(j, a):
    t, e = s.symbols("time parameter", real=True)
    B0 = s.diag(1, 7, 13)
    G = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 31
    D = s.Matrix([[2, -1, 1], [-1, 0, 2], [1, 2, -2]]) / 37
    F = s.Matrix([[0, 1, 2], [1, 1, -1], [2, -1, -1]]) / 41
    B = B0 + t * G + e * D + t * e * F + t * t * e * e * G / 4
    M = {
        (n, b): B0.inv()
        * B.applyfunc(lambda x, n=n, b=b: s.diff(x, t, n, e, b).subs({t: 0, e: 0}))
        for n in range(j + 1)
        for b in range(a + 1)
    }
    K = (B * B).applyfunc(s.expand)
    Q = (
        B0.inv()
        * K.applyfunc(lambda x: s.diff(x, t, j, e, a).subs({t: 0, e: 0}))
        * B0.inv()
    )
    right = Q
    for n in range(j + 1):
        for b in range(a + 1):
            if (n, b) in ((0, 0), (j, a)):
                continue
            right -= comb(j, n) * comb(a, b) * M[n, b] * M[j - n, a - b].T
    assert M[j, a] + M[j, a].T == right
    direct = s.Matrix(3, 3, lambda n, b: B0[b, b] * right[n, b] / (B0[n, n] + B0[b, b]))
    assert direct == M[j, a]
    assert sum(abs(x) ** 2 for x in M[j, a]) <= sum(abs(x) ** 2 for x in right)
    if (j, a) == (1, 1):
        assert Q != right


@pytest.mark.parametrize("momentum", (0, 1000, 10**6, 10**20))
@pytest.mark.parametrize("a", (1, 2))
def test_actual_frequency_parameter_jets(momentum, a):
    with mp.workdps(80):
        time = mp.mpf(1) / 5
        scale = (1 + time * time) ** 2
        k = mp.mpf(momentum) * mp.matrix([mp.mpf(2) / 3, mp.mpf(1) / 3, mp.mpf(2) / 3])
        G = mp.matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 20

        def omega(e):
            return mp.sqrt(1000**2 + (k.T * mp.expm(-e * G) * k)[0] / scale**2)

        epsilon = mp.mpf(1) / 200
        derivative = mp.diff(omega, epsilon, a)
        bound = mixed.constants()["omega"][0, a]
        assert abs(derivative) / omega(epsilon) <= int(bound)


@cache
def _input_jets():
    R = s.Matrix([[0, 1, -2], [-1, 0, 3], [2, -3, 0]]) / 37
    G = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 41
    D = s.Matrix([[2, -1, 1], [-1, 0, 2], [1, 2, -2]]) / 43
    omega, rotation, squeeze = {}, {}, {}
    for j in range(11):
        for a in range(3):
            omega[j, a] = s.Rational((-1) ** (j + a), j + a + 1)
            rotation[j, a] = R * s.Rational((-1) ** (j + a), j + a + 1)
            squeeze[j, a] = (G if (j + a) % 2 == 0 else D) * s.Rational(
                (-1) ** a, j + a + 1
            )
    omega[0, 0] = s.Integer(1000)
    return omega, rotation, squeeze


@cache
def _parameter_reference():
    return reference.parameter_reference(*_input_jets(), 10)


@cache
def _old_reference_at(epsilon):
    omega, R, S = _input_jets()
    w = [
        sum(omega[j, a] * epsilon**a / s.factorial(a) for a in range(3))
        for j in range(11)
    ]
    rotation = [
        sum((R[j, a] * epsilon**a / s.factorial(a) for a in range(3)), s.zeros(3))
        for j in range(11)
    ]
    squeeze = [
        sum((S[j, a] * epsilon**a / s.factorial(a) for a in range(3)), s.zeros(3))
        for j in range(11)
    ]
    return old_reference.reference_jets(w, rotation, squeeze, 10)


@pytest.mark.parametrize("order", range(1, 11))
@pytest.mark.parametrize("parameter_order", range(3))
def test_independent_parameter_differentiation_of_full_reference(
    order, parameter_order
):
    h = s.Rational(1, 10) ** 12
    base = _old_reference_at(s.Integer(0))[order][0]
    expected = _parameter_reference()[order, 0, parameter_order]
    if parameter_order == 0:
        difference = (base - expected).applyfunc(s.cancel)
        assert all(x == 0 for x in difference)
        return
    plus = _old_reference_at(h)[order][0]
    minus = _old_reference_at(-h)[order][0]
    actual = (
        (plus - minus) / (2 * h)
        if parameter_order == 1
        else (plus - 2 * base + minus) / h**2
    )
    difference = (actual - expected).applyfunc(s.cancel)
    assert sum(abs(x) ** 2 for x in difference) < s.Rational(1, 10) ** 40


@pytest.mark.parametrize("parameter_order", range(3))
def test_complete_mixed_reference_norm_majorants(parameter_order):
    exact = _parameter_reference()
    bounds = reference.constants()["ordered_mixed_coefficient_bounds"]
    for n in range(1, 11):
        for d in range(11 - n + 1):
            norm_squared = sum(abs(x) ** 2 for x in exact[n, d, parameter_order])
            assert norm_squared <= bounds[n, d, parameter_order] ** 2 / s.Integer(
                1000
            ) ** (2 * n)


def _mpmatrix(value):
    return mp.matrix(
        [
            [
                mp.mpf(str(s.N(s.re(value[i, j]), mp.mp.dps)))
                + 1j * mp.mpf(str(s.N(s.im(value[i, j]), mp.mp.dps)))
                for j in range(value.cols)
            ]
            for i in range(value.rows)
        ]
    )


@cache
def _covariance_inputs(seed):
    G = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]])
    D = s.Matrix([[2, -1, 1], [-1, 0, 2], [1, 2, -2]])
    B = s.Matrix([[0, 1, 2], [1, 1, -1], [2, -1, -1]])
    return [
        (G + s.I * D) / (100 + seed),
        (D + s.I * B) / (101 + seed),
        (B - s.I * G) / (103 + seed),
        (G - s.I * B) / (107 + seed),
    ]


@pytest.mark.parametrize("seed", range(3))
@pytest.mark.parametrize("order", range(4))
def test_independent_full_covariance_parameter_derivatives(seed, order):
    graph = _covariance_inputs(seed)
    expected = covariance.covariance_jet(graph)[order]
    with mp.workdps(70):
        matrices = [_mpmatrix(x) for x in graph]
        unit = mp.eye(3)

        def exact(e):
            r = mp.matrix(3)
            for j, value in enumerate(matrices):
                r += value * e**j / mp.factorial(j)
            C = (unit - r.H * r) ** -1
            F = mp.matrix(6, 3)
            for i in range(3):
                for j in range(3):
                    F[i, j] = (unit[i, j] + r[i, j]) / mp.sqrt(2)
                    F[i + 3, j] = -1j * (unit[i, j] - r[i, j]) / mp.sqrt(2)
            return (F * C * F.H).apply(lambda x: mp.re(x))

        actual = mp.diff(exact, mp.mpf(0), order)
        assert mp.norm(actual - _mpmatrix(expected)) < mp.mpf("1e-55")


@pytest.mark.parametrize("seed", range(3))
@pytest.mark.parametrize("order", (1, 2, 3))
def test_independent_covariance_Frechet_norms(seed, order):
    r, d, *_ = _covariance_inputs(seed)
    exact = covariance.covariance_jet([r, d, s.zeros(3), s.zeros(3)])[order]
    with mp.workdps(65):
        R, D = _mpmatrix(r), _mpmatrix(d)
        normr = mp.sqrt(max(mp.eighe(R.H * R, eigvals_only=True)))
        normd = mp.sqrt(max(mp.eighe(D.H * D, eigvals_only=True)))
        matrix = _mpmatrix(exact).apply(lambda x: mp.re(x))
        norm = max(abs(x) for x in mp.eigsy(matrix, eigvals_only=True))
        assert normr < mp.mpf(1) / 10
        assert norm < {1: 16, 2: 128, 3: 512}[order] * normd**order


@pytest.mark.parametrize("order", range(3))
@pytest.mark.parametrize("lower", (10**16, 10**17, 10**20))
def test_independent_infinite_parameter_tail_integral(order, lower):
    with mp.workdps(65):
        b = mp.mpf(99) / 100 / (mp.mpf(25) / 16) ** 2
        coefficient = mp.mpf(
            str(
                s.N(
                    covariance.constants()[
                        "actual_weighted_parameter_error_coefficients"
                    ][order],
                    65,
                )
            )
        )
        value = mp.quad(lambda x: x ** (5 - order), [0, 1])
        integral = (
            coefficient
            / (2 * mp.pi**2 * b ** mp.mpf("1.5"))
            * mp.mpf(lower) ** (-6 + order)
            * value
        )
        assert abs(value - mp.mpf(1) / (6 - order)) < mp.mpf("1e-60")
        assert integral < mp.mpf(str(s.N(covariance.TAIL[order], 65)))


@pytest.mark.parametrize(
    "epsilon",
    (
        s.Rational(-1, 100),
        s.Rational(-1, 10000),
        s.Integer(0),
        s.Rational(1, 10000),
        s.Rational(1, 100),
    ),
)
def test_independent_Taylor_integral_remainder_including_zero(epsilon):
    with mp.workdps(70):
        e = mp.mpf(str(s.N(epsilon, 70)))
        amplitude = mp.mpf("1e-28") / 4
        actual = amplitude * abs(mp.expm1(e) - e)
        upper = e * e * mp.mpf("1e-28") / 2
        assert actual <= upper
        if epsilon == 0:
            assert actual == upper == 0
        else:
            integral = (
                e * e * mp.quad(lambda x: (1 - x) * mp.exp(x * e), [0, 1]) * amplitude
            )
            assert abs(actual - integral) < mp.mpf("1e-90")


def test_quadratic_first_error_cannot_be_omitted():
    epsilon = s.Symbol("epsilon", real=True)
    S = s.Matrix([[1, 2, 0], [2, 3, 1], [0, 1, -1]]) / 100
    first = s.Matrix([[1, 0, 1], [0, -1, 2], [1, 2, 0]]) / 100
    hat = s.eye(3) / 100
    value = variation.generator(
        1000, s.zeros(3), S, hat + epsilon * first
    ) - variation.generator(1000, s.zeros(3), S, hat)
    actual = value.applyfunc(lambda x: s.diff(s.expand(x), epsilon, 2).subs(epsilon, 0))
    assert actual == -2 * first * S * first
    assert actual != s.zeros(3)


def test_varying_frequency_weight_cannot_be_omitted():
    epsilon = s.Symbol("epsilon", real=True)
    omega = 1000 * (1 + epsilon + epsilon**2)
    delta = 1 + 2 * epsilon + 3 * epsilon**2
    actual = s.diff(omega * delta, epsilon, 2).subs(epsilon, 0)
    wrong = (omega * s.diff(delta, epsilon, 2)).subs(epsilon, 0)
    assert actual == 12000 and wrong == 6000


@pytest.mark.parametrize("order", (0, -1, 11, True, False, 1.0, s.Integer(1)))
def test_mixed_reference_rejects_unsupported_order(order):
    with pytest.raises(ValueError):
        reference.parameter_reference(*_input_jets(), order)


def test_mixed_reference_rejects_incomplete_parameter_jets():
    omega, R, S = _input_jets()
    bad = dict(omega)
    del bad[10, 2]
    with pytest.raises(ValueError):
        reference.parameter_reference(bad, R, S, 10)


def test_tail_scope_remains_distinct_from_complete_response():
    text = str(audit.observable())
    assert "low-band" in text and "contacts" in text and "feedback" in text
    assert audit.matching()[-1]["status"] == audit.ITEM["status"]
