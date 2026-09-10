"""Independent regulator integrals, finite-part extraction and reference tests."""

import itertools

import mpmath as mp
import pytest
import sympy as sp
from p8_vacuum_fermion_outer_ms import audit, conversion, laurent, remainder


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_all_exact_identities(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", tuple(audit.gates()))
def test_all_proof_gates(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=lambda v: v if isinstance(v, str) else None
)
def test_unsupported_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def leading_regulated(e):
    return (
        mp.exp(2 * mp.euler * e)
        * 4 ** (-e)
        * mp.sqrt(mp.pi)
        / 2
        * (mp.mpf(3) / 2 - e)
        / (1 - e)
        * mp.gamma(e)
        * mp.gamma(2 * e)
        / mp.gamma(mp.mpf(5) / 2 + e)
    )


@pytest.mark.parametrize("digits", (6, 10))
def test_leading_Laurent_coefficients_with_complex_regulator(digits):
    with mp.workdps(80):
        h = mp.mpf(10) ** (-digits)
        value = leading_regulated(1j * h)
        finite = mp.re(value + 1 / (2 * h * h))
        expected = mp.mpf(47) / 18 + mp.pi**2 / 12
        assert abs(finite - expected) < 10000 * h * h
        assert abs(-h * mp.im(value) + mp.mpf(7) / 6) < 10000 * h * h


def test_independent_radial_representation_at_epsilon_one_half():
    with mp.workdps(60):

        def integrand(u):
            if abs(u) < mp.mpf("1e-8"):
                return sum(
                    (-1) ** n
                    * (mp.mpf(1) / (2 * n + 1) - mp.mpf(1) / (2 * n - 1))
                    * u ** (2 * n - 4)
                    for n in range(2, 10)
                )
            return ((1 + u * u) * (mp.atan(u) / u - 1) + u * u / 3) / u**4

        radial = mp.quad(integrand, [0, 1, mp.inf])
        direct = -2 * mp.exp(mp.euler) * radial
        assert abs(direct - leading_regulated(mp.mpf(1) / 2)) < mp.mpf("1e-35")
        assert abs(radial + mp.pi / 8) < mp.mpf("1e-35")


def compact_difference(e, r):
    def f(z):
        if z == 0 or z == 1:
            return 0
        k = r * z
        difference = (3 * k - 3 * k * k + k**3 - k ** (1 - e)) / (
            (1 - e) * (1 - k) ** 3
        )
        return z ** (2 * e - 1) * (1 - z) ** (mp.mpf(3) / 2 - e) * difference

    return mp.quad(f, [0, mp.mpf(1) / 2, 1])


def finite_correction(r):
    def f(z):
        if z == 0 or z == 1:
            return 0
        k = r * z
        pole = k * (2 - k) / (1 - k) ** 2
        logterm = -k * mp.log(k) / (1 - k) ** 3
        return (
            (1 - z) ** (mp.mpf(3) / 2)
            / z
            * ((3 - 4 * mp.log(2) + 2 * mp.log(z) - mp.log(1 - z)) * pole - logterm)
        )

    return mp.quad(f, [0, mp.mpf(1) / 2, 1])


@pytest.mark.parametrize("r_string", ("0.0625", "0.01", "1e-20"))
def test_exact_finite_mass_ratio_piece_from_complex_regulator(r_string):
    with mp.workdps(70):
        r = mp.mpf(r_string)
        h = mp.mpf("1e-8")
        e = 1j * h
        H = (
            mp.exp(2 * mp.euler * e)
            * 4 ** (-e)
            * mp.sqrt(mp.pi)
            / (2 * mp.gamma(mp.mpf(3) / 2 - e))
        )
        value = mp.gamma(e) * H * compact_difference(e, r)
        expected = finite_correction(r)
        pole = compact_difference(0, r)
        assert abs((mp.re(value) - expected) / r) < mp.mpf("1e-8")
        assert abs((-h * mp.im(value) - pole) / r) < mp.mpf("1e-10")
        assert 0 < pole < 3 * r
        assert abs(expected) < r * (20 - 2 * mp.log(r))


@pytest.mark.parametrize(
    "r,z",
    tuple(
        itertools.product(
            (sp.Rational(1, 16), sp.Rational(1, 100), sp.Rational(1, 10**20)),
            (sp.Rational(1, 1000), sp.Rational(1, 2), sp.Integer(1)),
        )
    ),
)
def test_uniform_rational_endpoint_prefactors(r, z):
    k = r * z
    h = 1 / (1 - k) ** 2 - 1
    assert 0 < h <= 3 * k
    assert 1 / (1 - k) ** 3 < 2


@pytest.mark.parametrize("m", (2, sp.Rational(9, 4), 10, 720, 10**20, 10**200))
def test_dyadic_log_enclosure_independently(m):
    d = remainder.enclosure(m, 1, 144)
    T = d["spectral_threshold"]
    n = d["least_dyadic_exponent"]
    assert 2 ** (n - 1) < T <= 2**n
    assert d["dimensionless_finite_correction_absolute_upper"] == (20 + 2 * n) / T
    assert d["finite_F_absolute_upper"] > d["finite_F_correction_absolute_upper"] > 0


@pytest.mark.parametrize(
    "k_string,e_string",
    tuple(itertools.product(("0.01", "0.0625", "0.5"), ("0.1", "0.25", "0.5"))),
)
def test_outer_Feynman_parameter_formula_independently(k_string, e_string):
    with mp.workdps(45):
        k, e = mp.mpf(k_string), mp.mpf(e_string)
        direct = mp.quad(lambda t: (t + (1 - t) * k) ** (-e), [0, 1])
        closed = (1 - k ** (1 - e)) / ((1 - k) * (1 - e))
        assert mp.almosteq(direct, closed)


@pytest.mark.parametrize("i", range(9))
def test_no_primitive_row_can_be_removed(i):
    rows = audit.frontier()
    del rows[i]
    with pytest.raises(ValueError):
        audit.validate_frontier(rows)


def test_finite_reference_changes_only_the_declared_family():
    rows = audit.frontier()
    changed = [
        r["id"]
        for r in rows
        if r["status"] == "BOUNDED_PAIRED_IN_COMMON_MS_INTERACTION_SCHEME"
    ]
    assert changed == ["scalar_Phi4_W2_F0"]
    assert sum(r["status"] == "UNEVALUATED" for r in rows) == 5
    assert sum(r["status"] == "BOUNDED_PAIRED_PRIMITIVE_ROW" for r in rows) == 2


def test_exact_counts_and_scope():
    assert len(audit.residuals()) == 54
    assert audit.scalar_entry_count() == 54
    assert len(audit.gates()) == 28
    assert audit.rejected_inputs() == 53
    assert len(audit.controls()) == 9
    assert "correction is not omitted" in laurent.data()["scope"]
    assert "separately owned" in conversion.data()["scope"]


def test_zero_Yukawa_and_conservative_wide_mass_domain():
    assert remainder.enclosure(2, 0, 144)["finite_F_absolute_upper"] == 0
    d = remainder.enclosure(2, 1, 144)
    assert d["dimensionless_finite_correction_absolute_upper"] > 1
    assert d["finite_F_absolute_upper"] > 5 * d["C_over_Q_upper"]
