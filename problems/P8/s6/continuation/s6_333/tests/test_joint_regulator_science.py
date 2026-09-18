"""Exact domains and independent quantitative joint-limit calibrations."""

from itertools import product

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_joint_soft_regulator import (
    audit,
    kernels,
    moments,
    regulator,
    source,
    uniform,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_residual(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_written_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_rejected_inputs_and_scope(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "bad",
    (
        None,
        True,
        False,
        "1",
        1.0,
        s.Float(1),
        [],
        {},
        s.Matrix([1]),
        s.I,
        s.oo,
        s.nan,
        s.Symbol("a"),
    ),
)
def test_exact_domains_reject_inexact_or_unresolved_values(bad):
    x, eta, e = s.Rational(1, 8), s.Rational(1, 64), s.Rational(1, 32)
    calls = [
        (kernels.regulator, (bad,)),
        (kernels.remaining_kernels, (bad, x)),
        (kernels.remaining_kernels, (e, bad)),
        (moments.conditional_log_square_moments, (bad, x)),
        (moments.conditional_log_square_moments, (1, bad)),
        (regulator.quotient_remainder_upper, (bad, e)),
        (regulator.quotient_remainder_upper, (x, bad)),
        (regulator.conditional_regulator_tv_upper, (1, x, bad)),
        (regulator.original_regulator_tv_upper, (bad, e)),
        (regulator.original_regulator_tv_upper, (x, bad)),
    ]
    for args in ((bad, x, eta, e), (1, bad, eta, e), (1, x, bad, e), (1, x, eta, bad)):
        calls.extend(
            ((uniform.regulated_cutoff_tv_upper, args), (uniform.joint_tv_upper, args))
        )
    for call, args in calls:
        with pytest.raises((TypeError, ValueError)):
            call(*args)


@pytest.mark.parametrize("e", (-1, -s.Rational(1, 100), s.Rational(1, 7), 1))
def test_regulator_closed_interval(e):
    with pytest.raises(ValueError):
        kernels.regulator(e)


@pytest.mark.parametrize(
    "a,x,eta",
    (
        (-1, s.Rational(1, 8), s.Rational(1, 64)),
        (s.Rational(11, 10), s.Rational(1, 8), s.Rational(1, 64)),
        (1, 0, s.Rational(1, 64)),
        (1, s.Rational(1, 7), s.Rational(1, 64)),
        (1, s.Rational(1, 8), 0),
        (1, s.Rational(1, 8), -1),
        (1, s.Rational(1, 8), s.Rational(1, 7)),
    ),
)
def test_ordered_cutoff_domain(a, x, eta):
    for call in (uniform.regulated_cutoff_tv_upper, uniform.joint_tv_upper):
        with pytest.raises(ValueError):
            call(a, x, eta, s.Rational(1, 32))


@pytest.mark.parametrize(
    "a,x,h,e",
    tuple(
        product(
            (s.S.Zero, s.Rational(1, 10**800), s.Rational(1, 5), s.S.One),
            (s.Rational(1, 10**60), s.Rational(1, 128), s.Rational(1, 8)),
            (s.Rational(1, 100), s.Rational(1, 4), s.S.One),
            (s.S.Zero, s.Rational(1, 32), s.Rational(1, 8)),
        )
    ),
)
def test_exact_joint_budgets_and_zero_cases(a, x, h, e):
    eta = h * x
    cut = 132000 * a * eta * (1 - s.log(eta)) ** 2 / source.KAPPA
    reg = 400000 * e * a * x * (1 - s.log(x)) ** 2 / source.KAPPA
    assert uniform.regulated_cutoff_tv_upper(a, x, eta, e) == cut
    assert regulator.conditional_regulator_tv_upper(a, x, e) == reg
    assert uniform.joint_tv_upper(a, x, eta, e) == cut + reg
    assert (
        regulator.quotient_remainder_upper(x, e)
        == 66000 * e * x * (1 - s.log(x)) ** 2 / source.KAPPA
    )
    values = moments.conditional_log_square_moments(a, x)
    assert all(not value.has(s.Float) for value in (*values, cut, reg))
    if a == 0:
        assert values == (0, 0, 0) and cut == reg == 0
    if e == 0:
        assert reg == 0 and kernels.remaining_kernels(e, x) == (1, s.log(x))
    assert regulator.quotient_remainder_upper(0, e) == 0


@pytest.mark.parametrize(
    "x,e",
    tuple(
        product(
            (s.Rational(1, 10**60), s.Rational(1, 128), s.Rational(1, 8)),
            (s.S.Zero, s.Rational(1, 32), s.Rational(1, 8)),
        )
    ),
)
def test_original_physical_kappa_power(x, e):
    eta = x / 4
    actual = uniform.original_joint_tv_upper(x, eta, e)
    expected = (
        106 * eta * (1 - s.log(eta)) ** 2 + 320 * e * x * (1 - s.log(x)) ** 2
    ) / s.Integer(10) ** 1597
    assert actual == expected
    assert (
        regulator.original_regulator_tv_upper(x, e)
        == 320 * e * x * (1 - s.log(x)) ** 2 / s.Integer(10) ** 1597
    )
    assert source.KAPPA == s.Integer(10) ** 800
    assert source.HEAVY_MASS2 == s.Integer(10) ** 200 / 512 + 2


@pytest.mark.parametrize(
    "estr,ystr,tstr",
    tuple(
        product(
            ("0.125", "0.0625", "1e-80"),
            ("0.01", "1e-20", "1e-1000"),
            ("0.001", "0.01"),
        )
    ),
)
def test_independent_regulated_kernel_inequalities(estr, ystr, tstr):
    with mp.workdps(220):
        e, y, t = mp.mpf(estr), mp.mpf(ystr), mp.mpf(tstr)
        g0, g1 = mp.power(y, 2 * e), mp.power(y + t, 2 * e)
        h0, h1 = (
            mp.expm1(2 * e * mp.log(y)) / (2 * e),
            mp.expm1(2 * e * mp.log(y + t)) / (2 * e),
        )
        lr = mp.log1p(t / y)
        assert 0 < g0 <= g1 < 1 and abs(h0) <= -mp.log(y)
        assert 0 <= h1 - h0 <= lr and 0 <= g1 - g0 <= 2 * e * lr
        assert abs(h0 - mp.log(y)) <= e * mp.log(y) ** 2
        exactg, exacth = kernels.remaining_kernels(s.Rational(estr), s.Rational(ystr))
        assert abs(mp.mpf(str(exactg.evalf(230))) - g0) < mp.mpf("1e-190")
        assert abs(mp.mpf(str(exacth.evalf(230))) - h0) < mp.mpf("1e-180")


@pytest.mark.parametrize(
    "astr,xstr",
    tuple(product(("0.1", "0.5", "1", "1e-800"), ("0.125", "0.001", "1e-10000"))),
)
def test_independent_all_three_conditional_log_moments(astr, xstr):
    with mp.workdps(120):
        a, x = mp.mpf(astr), mp.mpf(xstr)
        L, l = 1 - mp.log(x), -mp.log(x)
        actual = [
            mp.quad(lambda u: u**a * (L - mp.log(u)) ** 2 if u else 0, [0, 1]),
            mp.quad(
                lambda u: (
                    u**a * (L - mp.log(u)) * (l - mp.log(1 - u)) if 0 < u < 1 else 0
                ),
                [0, mp.mpf("0.5"), 1],
            ),
            mp.quad(
                lambda u: u**a * (l - mp.log(1 - u)) ** 2 if u < 1 else 0,
                [0, mp.mpf("0.5"), 1],
            ),
        ]
        aa, xx = s.Rational(astr), s.Rational(xstr)
        exact = moments.conditional_log_square_moments(aa, xx)
        for value, expression, bound in zip(actual, exact, (5, 3, 5)):
            scaled = mp.mpf(str((expression / (aa * xx)).evalf(130)))
            assert abs(value / scaled - 1) < mp.mpf("1e-95")
            assert 0 < value <= bound * L * L


@pytest.mark.parametrize("rstr", ("0.125", "0.001", "1e-60", "1e-10000"))
def test_independent_log_weighted_radial_majorant(rstr):
    with mp.workdps(150):
        R = mp.mpf(rstr)
        z = mp.log(mp.mpf(64000) / (100000 * R))
        length = 4 * z
        actual = mp.quad(lambda u: u, [0, length]) + mp.quad(
            lambda u: (length + u) * mp.exp(-u / 4), [0, mp.inf]
        )
        exact = 8 * z * z + 16 * z + 16
        assert abs(actual / exact - 1) < mp.mpf("1e-120")
        assert actual <= 16 * (1 - mp.log(R)) ** 2


@pytest.mark.parametrize(
    "rstr,estr",
    tuple(product(("0.125", "0.001", "1e-10000"), ("0.125", "0.01", "1e-80"))),
)
def test_saturating_majorant_regulator_model(rstr, estr):
    # An envelope calibration, NOT an actual physical current or uniform proof.
    with mp.workdps(240):
        R, e = mp.mpf(rstr), mp.mpf(estr)
        p = (4 * mp.pi) ** (-e) * mp.gamma(mp.mpf("1.5")) / mp.gamma(mp.mpf("1.5") + e)
        A = mp.gamma(mp.mpf("1.5") + e) / (mp.gamma(mp.mpf("1.5")) * mp.gamma(1 + e))
        c = e / (2 * (1 + e))

        def integral(aa, bb):
            length = 4 * mp.log(mp.mpf(bb) / (aa * R))
            finite = -mp.expm1(-e * length) / e + mp.exp(-e * length) / (
                e + mp.mpf("0.25")
            )
            return aa * R * finite, aa * R * (length + 4)

        ih, ih0 = integral(100000, 64000)
        iu, _ = integral(60000, 42000)
        k0, u0 = 50000 * R, 30000 * R
        ke = k0 + c * u0 + e * A * (ih + c * iu)
        delta = u0 / 2 + ih0 + (mp.euler - 2 - mp.log(mp.pi)) * k0
        error = abs((p * ke - k0) / e - delta) / (8 * mp.pi**2)
        assert error / (e * R) <= 66000 * (1 - mp.log(R)) ** 2


def test_not_a_constant_energy_independent_remainder():
    e = s.Rational(1, 32)
    assert regulator.quotient_remainder_upper(0, e) == 0
    with mp.workdps(80):
        large = mp.mpf(
            str(regulator.quotient_remainder_upper(s.Rational(1, 8), e).evalf(90))
        )
        small = mp.mpf(
            str(regulator.quotient_remainder_upper(s.Rational(1, 10**30), e).evalf(90))
        )
        assert 0 < small < large


def test_both_regulator_rates_and_all_historical_frontiers():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 189
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert "outer hard dimensional/evanescent" in audit.observable()["not_established"]
    assert "TV convergence" in audit.observable()["not_established"]
