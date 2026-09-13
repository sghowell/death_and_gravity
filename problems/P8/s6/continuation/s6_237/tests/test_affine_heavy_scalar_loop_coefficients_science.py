"""Independent moment integrals, analytic jets, full tails and coefficient checks."""

import itertools

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_heavy_scalar_loop_coefficients import (
    audit,
    coefficients,
    jets,
)
from p8_vacuum_affine_heavy_scalar_loop_remainder import audit as previous
from p8_vacuum_affine_heavy_scalar_tree_matching import model


def mpq(value):
    value = s.Rational(value)
    return mp.mpf(int(value.p)) / int(value.q)


def analytic_PQ(w, L, b):
    d = mp.sqrt(1 - 4 * L * w + 4 * L * b * w * w)
    a = (1 - 2 * L * w + d) / 2
    v = L * w / a
    beta = (L - b) * w / a
    h = (1 + beta) * mp.log1p(beta) / beta if beta else mp.mpf(1)
    return (1 + v) / d, ((1 + v) * mp.log((1 + d) / 2) - h) / d


def full_T_D(n, L, b):
    d = mp.sqrt(n * n - 4 * n * L + 4 * L * b)
    a = (n - 2 * L + d) / 2
    beta = (L - b) / a
    v = L / a
    h = (1 + beta) * mp.log1p(beta) / beta if beta else mp.mpf(1)
    hp = (beta - mp.log1p(beta)) / (beta * beta) if beta else mp.mpf("0.5")
    logarithm = mp.log(a + L) - mp.log(L)
    numerator = (1 + v) * logarithm - h
    dn = (n - 2 * L) / d
    an = (1 + dn) / 2
    T = numerator / d
    D = numerator * dn / (d * d) - (an / a) * (1 - v * logarithm + beta * hp) / d
    return T, D


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_every_exact_entry(name, value):
    assert s.cancel(value) == 0, name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_every_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_all_unsupported_inputs(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_log1p_helper_origin_and_shift_regression():
    w = jets.w
    assert jets.logarithm_one(0) == 0
    direct = jets.logarithm_one(w)
    assert s.expand(direct - w + w * w / 2 - w**3 / 3) == 0
    # A constant1 must never be passed where a small log1p increment is meant.
    assert jets.logarithm_one(s.Integer(1)) == s.Rational(5, 6)
    assert jets.derived()[1][0] == -1
    assert s.cancel(jets.derived()[1][0] + jets.logarithm_one(s.Integer(1)) + 1) != 0


@pytest.mark.parametrize(
    "L,b",
    [(s.Rational(2, 3), 0), (1, s.Rational(1, 3)), (1, 1), (s.Rational(3, 2), -1)],
)
def test_PQ_coefficients_from_independent_full_analytic_derivatives(L, b):
    with mp.workdps(70):
        ll, bb = mpq(L), mpq(b)
        P, Q = jets.derived()[:2]
        for which, expected in enumerate((P, Q)):
            actual = mp.taylor(
                lambda w, index=which: analytic_PQ(w, ll, bb)[index], 0, 3
            )
            for j in range(4):
                wanted = (
                    mpq(expected[j].subs({jets.L: L, jets.b: b}))
                    if isinstance(expected[j], s.Basic)
                    else mp.mpf(expected[j])
                )
                assert abs(actual[j] - wanted) < mp.mpf("1e-55")


@pytest.mark.parametrize("parameter", [s.Rational(1, 4), 2, s.Rational(19, 8)])
@pytest.mark.parametrize("j", [1, 2, 3])
def test_full_log_moment_recurrence_against_direct_quadrature(parameter, j):
    with mp.workdps(70):
        S = mpq(parameter) + mp.j / 16
        L = lambda x: 1 - S * x * (1 - x)
        M = lambda k: mp.quad(lambda x: L(x) ** k, [0, mp.mpf("0.5"), 1])
        J = lambda k: mp.quad(lambda x: L(x) ** k * mp.log(L(x)), [0, mp.mpf("0.5"), 1])
        actual = J(j)
        reduced = (
            2 * j * (1 - S / 4) * J(j - 1) - 2 * M(j) + 2 * (1 - S / 4) * M(j - 1)
        ) / (2 * j + 1)
        assert abs(actual - reduced) < mp.mpf("1e-60")
        assert abs(actual) > mp.mpf("1e-4")


@pytest.mark.parametrize("degree", range(4))
def test_independent_heavy_and_light_polynomial_moments(degree):
    x, S = s.symbols("x S")
    moment = s.integrate((x * (1 - x)) ** degree, (x, 0, 1))
    assert moment == s.factorial(degree) ** 2 / s.factorial(2 * degree + 1)
    light = sum(
        s.binomial(degree, k) * (-S) ** k * s.factorial(k) ** 2 / s.factorial(2 * k + 1)
        for k in range(degree + 1)
    )
    assert (
        s.expand(s.integrate((1 - S * x * (1 - x)) ** degree, (x, 0, 1)) - light) == 0
    )


@pytest.mark.parametrize("phase_L,phase_b", [(0, 0), (0, 3), (2, 5), (5, 1), (7, 7)])
def test_complex_parameter_Cauchy_circle_not_only_real_parameters(phase_L, phase_b):
    with mp.workdps(70):
        L = mp.mpf("1.8") * mp.exp(mp.j * mp.pi * phase_L / 4)
        b = mp.mpf("1.7") * mp.exp(mp.j * mp.pi * phase_b / 4)
        for j in range(16):
            w = mp.exp(2 * mp.j * mp.pi * j / 16) / 32
            P, Q = analytic_PQ(w, L, b)
            assert abs(P) < 2 and abs(Q) < 4


@pytest.mark.parametrize("n", [128, 10**6, model.MASS2])
@pytest.mark.parametrize("choice", [0, 1, 2])
def test_complete_four_order_primitive_tails_for_complex_parameters(n, choice):
    with mp.workdps(1000 if n == model.MASS2 else 90):
        nn = mpq(n)
        fixtures = [
            (1 + mp.j / 3, mp.mpf("0.5") + mp.j / 4),
            (mp.mpf("0.5") - mp.j / 5, -mp.mpf("0.7")),
            (mp.mpf("1.5"), mp.mpf("1.5")),
        ]
        L, b = fixtures[choice]
        P = [1, 3 * L, 10 * L * L - 2 * L * b, 35 * L**3 - 15 * L * L * b]
        Q = [
            -1,
            (b - 7 * L) / 2,
            (-74 * L * L + 28 * L * b + b * b) / 6,
            (-533 * L**3 + 327 * L * L * b - 9 * L * b * b + b**3) / 12,
        ]
        H = mp.log(nn) - mp.log(L)
        tp = sum((P[j] * H + Q[j]) / nn ** (j + 1) for j in range(4))
        dp = sum(
            ((j + 1) * P[j] * H + (j + 1) * Q[j] - P[j]) / nn ** (j + 2)
            for j in range(4)
        )
        T, D = full_T_D(nn, L, b)
        assert 0 < abs(T - tp) < 2**22 * 500 / nn**5
        assert 0 < abs(D - dp) < 3 * 2**23 * 500 / nn**6


@pytest.mark.parametrize(
    "v,t",
    [
        (0, 0),
        (s.Rational(1, 4), 0),
        (-s.Rational(1, 4), s.I / 4),
        (s.Rational(3, 20) + s.I / 5, s.Rational(1, 4)),
        (-s.I / 4, -s.Rational(1, 4)),
        (s.I / 4, s.I / 4),
    ],
)
def test_joint_complex_domain_and_every_parameter_denominator(v, t):
    vr, vi = s.re(v), s.im(v)
    tr, ti = s.re(t), s.im(t)
    vv, tt = audit.require_bidisk(vr, vi, tr, ti)
    assert vv == v and tt == t
    channels = (2 + v - t / 2, t, 2 - v - t / 2)
    n = model.MASS2
    for channel in channels:
        assert s.simplify(channel * s.conjugate(channel)) <= s.Rational(19, 8) ** 2
        for x in (0, s.Rational(1, 3), s.Rational(1, 2), 1):
            L = 1 - channel * x * (1 - x)
            assert s.re(L) >= s.Rational(13, 32)
            for other in channels:
                b = other * s.Rational(1, 4)
                for z in (0, s.Rational(1, 10), 1):
                    denominator = n * z + (1 - z) ** 2 * L - b * z * z
                    assert s.re(denominator) > 0


def test_all_complete_product_norms_independently_enumerated():
    a = [0, 2, 7, 13, 27]
    w = s.Rational(1, 64)
    bubble = sum(
        s.Rational(a[i] * a[j], 2) * w ** (i + j - 6)
        for i, j in itertools.product(range(1, 5), repeat=2)
        if i + j >= 6
    )
    triangle = 2 * sum(
        a[i] * 2 * 32 ** (j - 1) * w ** (i + j - 6)
        for i, j in itertools.product(range(1, 5), repeat=2)
        if i + j >= 6
    )
    assert bubble == s.Rational(2286169, 8192) < 512
    assert triangle == 1003424 < 2**20
    assert 6 * 2**22 + 1024 + 2**20 + 2 * 3 * 2**23 + 2048 < 2**27
    assert 3 * 2**27 * 500 < 10**12


def test_Cauchy_coefficient_factors_and_no_factorial_error():
    v, t = s.symbols("v t")
    a, b, c = s.symbols("a b c")
    expression = a * v * v + b * v * v * t + c * v**4
    assert s.diff(expression, v, 2).subs({v: 0, t: 0}) / 2 == a
    assert s.diff(expression, v, 2, t).subs({v: 0, t: 0}) / 2 == b
    assert s.diff(expression, v, 4).subs({v: 0, t: 0}) / 24 == c
    assert [s.Rational(1, 4) ** (-k) for k in (2, 3, 4)] == [16, 64, 256]


def test_actual_coefficient_sign_margins_and_retained_scope():
    n, g2 = model.MASS2, model.G2
    assert s.Rational(1500, 288) * g2 / n < s.Rational(1, 10**203)
    assert s.Rational(1500, 432) * g2 / n < s.Rational(1, 10**203)
    assert s.Rational(8 * 10**12, 9) * g2 / n < s.Rational(1, 10**192)
    assert audit.frontier() == previous.frontier()
    assert audit.matching()[:-1] == previous.matching()
    assert len(audit.matching()) == 93
    assert "no sign" in audit.observable()["first_loop_matching"]
    assert "truncated" in audit.observable()["specified_truncation"]
    assert len({name for name, _, _ in audit.bad_cases()}) == len(audit.bad_cases())
    for name in ("b20", "b21", "b40"):
        assert audit.require_jet(name, 1) == (name, 1)
    constant = s.Symbol("fixed_OS4_contact")
    assert s.diff(constant, jets.v, 2) == 0
    assert "delta_b40" in coefficients.data()["complete_coefficient_remainder_bounds"]


def complete_loop_quadrature(n, v, t, grid):
    """Complete unexpanded S235 formula; quadrature is a diagnostic only."""
    channels = (2 + v - t / 2, t, 2 - v - t / 2)
    contact = -3 / (n - 2) + 2 / (n - 2) ** 2
    fixed = lambda a: contact + 1 / (n - a)
    result = 0
    for i, a in enumerate(channels):
        B = 0
        C = 0
        boxes = 0
        for x, wx in grid:
            L = 1 - a * x * (1 - x)
            B -= wx * mp.log(L)
            C += wx * full_T_D(n, L, 0)[0]
            for j, b in enumerate(channels):
                if i != j:
                    if b == 0:
                        boxes += wx * full_T_D(n, L, 0)[1]
                    else:
                        boxes += wx * sum(
                            wy * full_T_D(n, L, b * y * (1 - y))[1] for y, wy in grid
                        )
        result += fixed(a) ** 2 * B / 2 + 2 * fixed(a) * C + boxes
    return result / (16 * mp.pi**2)


def test_three_coefficients_from_complete_unexpanded_loop_quadrature():
    with mp.workdps(85):
        n = mp.mpf(10) ** 16
        gx, gw = mp.gauss_quadrature(40, "legendre")
        grid = [((gx[j] + 1) / 2, gw[j] / 2) for j in range(40)]
        fun = lambda v, t: complete_loop_quadrature(n, v, t, grid)
        actual20 = mp.diff(lambda v: fun(v, mp.mpf(0)), 0, 2) / 2
        actual21 = mp.diff(fun, (mp.mpf(0), mp.mpf(0)), (2, 1)) / 2
        actual40 = mp.diff(lambda v: fun(v, mp.mpf(0)), 0, 4) / 24
        H = mp.log(n)
        approximate20 = (
            (3 * H - mp.mpf(157) / 45) / n**4 + (24 * H - mp.mpf(7969) / 315) / n**5
        ) / (16 * mp.pi**2)
        approximate21 = (-3 * H + mp.mpf(13) / 21) / (16 * mp.pi**2 * n**5)
        error = mp.mpf(10) ** 12 / (16 * mp.pi**2 * n**6)
        assert abs(actual20 - approximate20) < 16 * error
        assert abs(actual21 - approximate21) < 64 * error
        assert abs(actual40) < 256 * error
        assert actual20 > 0 and actual21 < 0
