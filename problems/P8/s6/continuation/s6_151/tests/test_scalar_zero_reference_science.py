"""Independent sector/hypergeometric integrals and regulator-first finite parts."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_scalar_zero_reference import (
    anchors,
    audit,
    bounds,
    calibration,
    conversion,
    sector,
)


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_exact_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[v[0] for v in audit.bad_cases()]
)
def test_invalid_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def sector_R(e):
    return mp.quad(
        lambda t: mp.quad(
            lambda r: (
                r ** (e - 1)
                * (1 + t) ** (-2 + e)
                * mp.expm1(
                    (1 - 2 * e) * mp.log1p(r * (1 + t))
                    - (2 - e) * mp.log1p(r * t / (1 + t))
                )
            ),
            [0, 1],
        ),
        [0, 1],
    )


def sector_A(e):
    P = mp.log(2) if e == 1 else (1 - 2 ** (e - 1)) / (1 - e)
    return 6 * (P / e + sector_R(e))


def stable_hyper(e, x):
    a = x * (1 - x)
    first = (1 - a) ** (-e) * (1 + e * a / ((1 - e) * (1 - a)))
    return 2 * (1 + e) * mp.gamma(1 + e) ** 2 / mp.gamma(1 + 2 * e) * first - (
        1 + e
    ) / (1 - e) * a**e * mp.hyp2f1(2, 2 * e, 1 + e, a)


@pytest.mark.parametrize("e", (".2", ".4", ".7", "1"))
def test_six_sector_integral_against_independent_bubble_radial_master(e):
    with mp.workdps(32):
        e = mp.mpf(e)
        direct = (
            3
            / (e * (1 + e))
            * mp.quad(
                lambda x: mp.hyp2f1(e, 2 - e, 2 + e, 1 - x * (1 - x)), [0, 0.5, 1]
            )
        )
        assert abs(sector_A(e) - direct) < mp.mpf("1e-24")


@pytest.mark.parametrize("e", (".1", "-.07+.08j", ".06-.04j"))
@pytest.mark.parametrize("x", (".2", ".01", ".000001"))
def test_stable_hypergeometric_connection(e, x):
    with mp.workdps(40):
        e = complex(e)
        e = mp.mpc(str(e.real), str(e.imag))
        x = mp.mpf(x)
        direct = mp.hyp2f1(e, 2 - e, 2 + e, 1 - x * (1 - x))
        assert abs(stable_hyper(e, x) - direct) < mp.mpf("1e-30")


@cache
def finite_sector_constant():
    with mp.workdps(38):

        def integrand(r, t):
            v, w = 1 + t + r * t, 1 + r * (1 + t)
            B = w / v**2
            D = (1 + t + t * t + t**3 - r * t * t) / ((1 + t) ** 2 * v * v)
            Udiff = mp.log1p(r * t / (1 + t)) - 2 * mp.log1p(r * (1 + t))
            return mp.log(r) * D + D * mp.log1p(t) + B * Udiff / r

        R1 = mp.quad(lambda t: mp.quad(lambda r: integrand(r, t), [0, 1]), [0, 1])
        p2 = (1 - mp.log(2) - mp.log(2) ** 2 / 2) / 2
        return p2 + R1


@cache
def independent_circle_samples():
    with mp.workdps(40):
        out = []
        for k in range(32):
            e = mp.mpf(".06") * mp.exp(2j * mp.pi * k / 32)
            # x=u^4/2 resolves both symmetric hypergeometric endpoints.
            H = 2 * mp.quad(
                lambda u, e=e: 2 * u**3 * stable_hyper(e, mp.mpf(".5") * u**4), [0, 1]
            )
            bare_e2 = mp.exp(2 * mp.euler * e) * mp.gamma(1 + 2 * e) * H / (2 * (1 + e))
            I0_times_e = mp.exp(mp.euler * e) * mp.gamma(1 + e)
            out.append((e, bare_e2, I0_times_e))
        return out


@pytest.mark.parametrize("ell", (0, 2, -1))
def test_independent_complex_finite_parts_keep_pole_products(ell):
    with mp.workdps(38):
        j = finite_sector_constant()
        expected = (
            (mp.mpf(".5"), ell + mp.mpf(".5"), j + ell + ell * ell + mp.pi**2 / 6),
            (mp.mpf("-.5"), mp.mpf(".5") - ell, j + ell - ell * ell),
            (mp.mpf("-.5"), mp.mpf(".5"), j + ell + ell * ell / 2 + mp.pi**2 / 12),
        )
        values = [[], [], []]
        for e, bare0, I00 in independent_circle_samples():
            bare = mp.exp(2 * ell * e) * bare0
            I0 = mp.exp(ell * e) * I00
            for i, value in enumerate((bare, bare - I0 * I0, bare - I0)):
                values[i].append((e, value))
        for i in range(3):
            for n in range(3):
                coeff = sum(v / e**n for e, v in values[i]) / 32
                assert abs(coeff - expected[i][n]) < mp.mpf("1e-23")


@pytest.mark.parametrize("r", (".0001", ".1", ".5", "1"))
@pytest.mark.parametrize("t", ("0", ".3", "1"))
def test_pointwise_finite_sector_majorant(r, t):
    with mp.workdps(35):
        r, t = mp.mpf(r), mp.mpf(t)
        v, w = 1 + t + r * t, 1 + r * (1 + t)
        B = w / v**2
        D = (1 + t + t * t + t**3 - r * t * t) / ((1 + t) ** 2 * v * v)
        ud = mp.log1p(r * t / (1 + t)) - 2 * mp.log1p(r * (1 + t))
        integrand = mp.log(r) * D + D * mp.log1p(t) + B * ud / r
        assert 0 < D <= 4 and 0 < B <= 3
        assert abs(ud) / r <= mp.mpf("4.5")
        assert abs(integrand) <= 4 * abs(mp.log(r)) + mp.mpf("17.5")


def test_finite_sector_constant_and_zero_endpoint():
    with mp.workdps(30):
        assert abs(sector_R(0) - mp.log(2) / 2) < mp.mpf("1e-25")
        j = finite_sector_constant()
        assert abs(j) < 22
        assert abs(j - mp.mpf("-1.49442065276884266353698872776")) < mp.mpf("1e-27")


@pytest.mark.parametrize("e", (".4", ".7"))
def test_equal_mass_vacuum_derivative_normalization(e):
    with mp.workdps(30):
        e = mp.mpf(e)
        A = sector_A(e)
        S = lambda b: (
            mp.exp(2 * mp.euler * e) * mp.gamma(-1 + 2 * e) * b ** (1 - 2 * e) * A
        )
        derivative = -mp.diff(S, 1) / 3
        master = mp.exp(2 * mp.euler * e) * mp.gamma(2 * e) * A / 3
        assert abs(derivative - master) < mp.mpf("1e-24")


@pytest.mark.parametrize("ell", (1, 3, -2))
def test_regulated_nonzero_momentum_conversion_independently(ell):
    e, z = s.symbols("e z")
    I0 = 1 / e + ell + e * (s.Rational(3, 2) + ell * ell / 2)
    IR = 1 / (7 - z) - s.Rational(1, 7)
    I = I0 + IR + e * z / (13 - z)
    J0 = s.Rational(1, 2) / e**2 + (ell + s.Rational(1, 2)) / e + 5 + ell * ell
    J = J0 + z / (5 - z) + (I - I0) / e
    K = -s.Rational(1, 2) / e**2 + s.Rational(1, 2) / e
    old = J - I0 * I - J0 + I0**2
    new = J - I / e - K
    new0 = J0 - I0 / e - K
    assert s.factor(new - old - new0 - (I0 - 1 / e) * (I - I0)) == 0
    assert s.simplify(s.limit(new - old - new0, e, 0) - ell * IR) == 0
    assert s.diff(ell * IR, z, 2) != 0


@pytest.mark.parametrize("ell", (0, 2, 5))
def test_early_one_loop_epsilon_truncation_is_detected(ell):
    e = s.symbols("e")
    I0 = 1 / e + ell + e * (ell**2 / 2 + s.pi**2 / 12)
    wrong = (1 / e + ell) / e
    assert s.simplify(s.limit(I0 / e - wrong, e, 0) - (ell**2 / 2 + s.pi**2 / 12)) == 0
    assert ell**2 / 2 + s.pi**2 / 12 != 0


@pytest.mark.parametrize("M", (40, 100, 1000))
def test_complete_triangle_integral_and_its_bound(M):
    with mp.workdps(35):
        M = mp.mpf(M)
        integral = mp.quad(
            lambda y: y / ((y + mp.mpf(".25")) ** 2 * (y + M)), [0, 1, M, mp.inf]
        )
        exact = M * mp.log(4 * M) / (M - mp.mpf(".25")) ** 2 - 1 / (M - mp.mpf(".25"))
        assert abs(integral - exact) < mp.mpf("1e-28")
        assert 0 < integral < 2 * mp.log(4 * M) / M


@pytest.mark.parametrize("theta", ("0", ".3", "1", "2", "3"))
def test_outer_light_bubble_on_complex_disc(theta):
    with mp.workdps(30):
        z = 2 + mp.exp(1j * mp.mpf(theta))
        IR = mp.quad(lambda x: -mp.log(1 - x * (1 - x) * z), [0, 1])
        assert abs(IR) < mp.log(4) < 2


def test_full_heavy_numerator_and_local_tree_shift():
    L, g, M, z, F = s.symbols("L g M z F")
    C = -L + g / (M - z)
    TA, TB = g / (M + 2 - z), g / (M + 3 + z)
    CQ = -L + g / (M + 5)
    complete = (C + TA + TB) * CQ**2
    remainder = L**2 * (TA + TB) + (C + TA + TB) * (CQ**2 - L**2)
    assert s.expand(complete - C * L**2 - remainder) == 0
    for missing in (TA, TB):
        assert s.factor(s.diff(missing * CQ**2, z, 2)) != 0
    b2 = 2 * g / (M - 2) ** 3
    shift = L * L * g * F * s.diff(b2, g)
    assert s.factor(shift / b2 - L * L * F) == 0


def test_actual_enclosure_and_scope():
    d = calibration.data()
    assert all(d["bounds"].values())
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 9
    assert sum(v["status"] != "OPEN" for v in audit.matching()) == 2
    assert audit.controls()["rejected_inputs"] == 88
    assert len(audit.gates()) == 23
    assert (
        bounds.enclosure(s.Rational(1, 100), 1, 0, 0)["added_complete_scale_term_upper"]
        == 0
    )
    assert anchors.data()["checks"]["proper_MS_poles_independent_of_scale"] == 0
    assert (
        conversion.data()["checks"]["regulated_core_difference_before_finite_parts"]
        == 0
    )
    assert sector.data()["finite_constant_absolute_upper"] == 22
