"""Independent regulated slopes, coefficient maps and full field covariance."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_finite_field_covariance import (
    audit,
    bounds,
    calibration,
    coefficients,
    covariance,
    kernel,
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


def raw_trace(sval, e, m, mu):
    bubble = (
        mp.gamma(e)
        * mp.exp(mp.euler * e)
        * mp.quad(lambda x: (mu**2 / (m**2 - x * (1 - x) * sval)) ** e, [0, 1])
    )
    tad = mp.exp(mp.euler * e) * mp.gamma(e - 1) * m**2 * (mu**2 / m**2) ** e
    return (4 * m**2 - sval) * bubble - 2 * tad


@pytest.mark.parametrize("m", (2, 40, 1000))
@pytest.mark.parametrize("e", (".2", ".4"))
@pytest.mark.parametrize("sval", (0, 1))
def test_dimensional_trace_derivative_before_setting_scale(m, e, sval):
    with mp.workdps(35):
        m, e, mu = mp.mpf(m), mp.mpf(e), mp.mpf(7)
        actual = mp.diff(lambda z: raw_trace(z, e, m, mu), sval)
        reduced = (
            -2
            * mp.exp(mp.euler * e)
            * mp.gamma(e)
            * (3 - 2 * e)
            * mp.quad(
                lambda x: x * (1 - x) * (mu**2 / (m**2 - x * (1 - x) * sval)) ** e,
                [0, 1],
            )
        )
        assert abs(actual - reduced) < mp.mpf("1e-27")


@pytest.mark.parametrize("m", (2, 40, 1000))
def test_complex_regulator_fixes_fermion_first_epsilon_coefficient(m):
    with mp.workdps(36):
        m = mp.mpf(m)
        a = lambda x: x * (1 - x)
        ellx = lambda x: -mp.log1p(-a(x) / m**2)
        J1 = mp.quad(lambda x: a(x) * ellx(x), [0, 1])
        J2 = mp.quad(lambda x: a(x) * ellx(x) ** 2, [0, 1])
        zero = mp.mpf(1) / 3 - 3 * J1
        first = 2 * J1 - mp.mpf("1.5") * J2 - mp.pi**2 / 24
        circle0, circle1 = 0, 0
        for j in range(24):
            e = mp.mpf(".03") * mp.exp(2j * mp.pi * j / 24)
            regular = (
                mp.exp(mp.euler * e)
                * mp.gamma(1 + e)
                * (3 - 2 * e)
                * mp.quad(lambda x, e=e: a(x) * mp.exp(e * ellx(x)), [0, 1])
            )
            fp = -(regular - mp.mpf(".5")) / e
            circle0 += fp / 24
            circle1 += fp / e / 24
        assert abs(circle0 - zero) < mp.mpf("1e-27")
        assert abs(circle1 - first) < mp.mpf("1e-26")
        direct_slope = (
            mp.diff(
                lambda z: (
                    (4 * m**2 - z)
                    * mp.quad(lambda x: -mp.log1p(-a(x) * z / m**2), [0, 1])
                ),
                1,
            )
            / 2
        )
        assert abs(direct_slope - zero) < mp.mpf("1e-28")


@pytest.mark.parametrize("M", (40, 1000, 1000000))
def test_scalar_epsilon_coefficient_uses_complete_dimensional_slope(M):
    with mp.workdps(35):
        M = mp.mpf(M)
        ell = mp.log(M) + 2
        D = lambda x: x * M + (1 - x) ** 2
        b = lambda x: x * (1 - x) / D(x)
        grid = [0, 1 / M, 1]
        direct = mp.diff(
            lambda e: (
                mp.exp(mp.euler * e + ell * e)
                * mp.gamma(1 + e)
                * mp.quad(lambda x: D(x) ** (-e) * b(x), grid)
            ),
            0,
        )
        moment = mp.quad(lambda x: b(x) * (ell - mp.log(D(x))), grid)
        assert abs(direct - moment) < mp.mpf("1e-27")
        assert 0 < moment < ell / (2 * M)


@pytest.mark.parametrize("m", (36, 100, 1000))
def test_positive_log_moments_enclose_fermion_epsilon_coefficient(m):
    with mp.workdps(32):
        m = mp.mpf(m)
        a = lambda x: x * (1 - x)
        log = lambda x: -mp.log1p(-a(x) / m**2)
        b = 1 / (4 * m**2 - 1)
        J1 = mp.quad(lambda x: a(x) * log(x), [0, 1])
        J2 = mp.quad(lambda x: a(x) * log(x) ** 2, [0, 1])
        fp1 = 2 * J1 - mp.mpf("1.5") * J2 - mp.pi**2 / 24
        assert 0 < J1 < b / 6
        assert 0 < J2 < b**2 / 6
        assert abs(fp1) < b / 3 + b**2 / 4 + mp.mpf(2) / 3 < 1


@pytest.mark.parametrize("w", (0, mp.mpf(".5"), 1, 2))
def test_independent_double_circle_finite_total_coefficient_map(w):
    with mp.workdps(35):
        X, p1, p22, p21 = map(mp.mpf, (2, 3, -1, 4))
        k0, k1, k2, t0, t1 = map(mp.mpf, (".2", "-.3", ".4", "-.1", ".7"))
        wanted = (w * (w + 1) * k0**2 / 2 - w * t0) * X - w * k1 * p1
        finite = 0
        for j in range(16):
            e = mp.mpf(".04") * mp.exp(2j * mp.pi * j / 16)
            k = k0 + e * k1 + e**2 * k2
            t = t0 + e * t1

            def total(h, e=e, k=k, t=t):
                return (X + h * p1 / e + h**2 * (p22 / e**2 + p21 / e)) / (
                    1 + h * k + h**2 * t
                ) ** w

            # Differentiate the literal unexpanded expression in loop order,
            # then take its finite Laurent coefficient on a separate circle.
            second = mp.diff(total, 0, 2) / 2
            finite += second / 16
        assert abs(finite - wanted) < mp.mpf("1e-27")


@pytest.mark.parametrize(
    "valences,internal",
    (
        ((4,), 0),
        ((2, 2), 0),
        ((4, 4), 2),
        ((4, 2, 2), 2),
        ((2, 2, 2, 2), 2),
        ((1, 1, 1, 1), 0),
    ),
)
@pytest.mark.parametrize("K", (s.Rational(2, 3), s.Rational(7, 5)))
def test_literal_vertex_and_propagator_factors(valences, internal, K):
    assert sum(valences) == 2 * internal + 4
    vertex = s.prod(K ** (-s.Rational(n, 2)) for n in valences)
    assert s.simplify(vertex * K**internal - K**-2) == 0


@pytest.mark.parametrize("k", (s.Rational(1, 7), s.Rational(-2, 9)))
def test_correlated_mass_and_kinetic_line_insertion_preserves_mass_one(k):
    h, p2 = s.symbols("h p2")
    K = 1 + h * k
    transformed_inverse = (p2 + 1) / K
    transformed_propagator = 1 / transformed_inverse
    assert s.factor(transformed_propagator - K / (p2 + 1)) == 0
    assert s.diff(transformed_propagator, h).subs(h, 0) == k / (p2 + 1)
    assert transformed_inverse.subs(p2, -1) == 0
    wrong_inverse = p2 / (1 + h * k) + 1
    assert s.factor(wrong_inverse.subs(p2, -1)) != 0


@pytest.mark.parametrize("k1", (s.Rational(1, 3), s.Rational(-2, 7)))
def test_raw_and_counterterm_epsilon_products_cancel_only_together(k1):
    e, k0, p, F = s.symbols("e k0 pole finite")
    field = -2 * (k0 + e * k1)
    raw = p / e + F
    ct = -p / e
    rawfinite = s.expand(field * raw).coeff(e, 0)
    ctfinite = s.expand(field * ct).coeff(e, 0)
    assert s.expand(rawfinite + ctfinite + 2 * k0 * F) == 0
    assert s.expand(rawfinite - (-2 * k0 * F)) == -2 * k1 * p
    assert ctfinite == 2 * k1 * p
    assert k1 != 0


@pytest.mark.parametrize("weight", (s.Rational(1, 2), s.Integer(1)))
def test_fundamental_coupling_square_includes_first_order_square(weight):
    h, x, k0, k1, t0, p = s.symbols("h x k0 k1 t0 p")
    x1 = -weight * k0 * x
    x2 = (weight * (weight + 1) * k0**2 / 2 - weight * t0) * x - weight * k1 * p
    actual = s.diff((x + h * x1 + h**2 * x2) ** 2, h, 2).subs(h, 0) / 2
    assert s.expand(actual - 2 * x * x2 - x1**2) == 0
    assert s.expand(actual - 2 * x * x2) == weight**2 * k0**2 * x**2


def test_second_slope_not_falsely_completed_and_no_double_counting():
    d = calibration.data()
    assert all(d["bounds"].values())
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 13
    assert sum(r["status"] != "OPEN" for r in audit.matching()) == 6
    assert len(audit.residuals()) == 62
    assert len(audit.gates()) == 23
    assert audit.controls()["rejected_inputs"] == 104
    assert "new_second_slope" in str(
        coefficients.data()["second_finite_coefficient_map"]
    )
    assert covariance.data()["new_second_normalization_amplitude_kept_separate"] != 0
    assert kernel.data()["complete_Phi_k_first_epsilon"] != 0
    assert bounds.enclosure(1, 1, 40, 1, 0)["one_loop_amputated_relative_upper"] > 0
