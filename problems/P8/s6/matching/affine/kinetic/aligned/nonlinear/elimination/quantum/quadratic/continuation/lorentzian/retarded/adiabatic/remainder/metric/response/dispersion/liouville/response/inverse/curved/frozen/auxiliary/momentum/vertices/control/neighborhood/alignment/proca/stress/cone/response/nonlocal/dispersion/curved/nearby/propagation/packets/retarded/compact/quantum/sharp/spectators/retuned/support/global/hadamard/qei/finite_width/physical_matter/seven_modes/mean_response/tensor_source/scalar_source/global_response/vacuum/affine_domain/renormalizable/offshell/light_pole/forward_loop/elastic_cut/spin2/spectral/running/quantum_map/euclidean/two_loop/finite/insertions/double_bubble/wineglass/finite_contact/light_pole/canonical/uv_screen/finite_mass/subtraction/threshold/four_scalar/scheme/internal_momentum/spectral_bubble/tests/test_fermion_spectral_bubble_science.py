"""Independent spectral fractions, crossing factors and local-scope controls."""

from fractions import Fraction as F

import pytest
import sympy as sp
from p8_vacuum_fermion_spectral_bubble import (
    audit,
    calibration,
    density,
    outer,
    ownership,
    spectral,
)


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_named_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_every_explicit_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[v[0] for v in audit.bad_cases()]
)
def test_all_invalid_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("m", (2, 36, 10**200))
@pytest.mark.parametrize(
    "Y,L,Q", ((0, 1, 144), (1, 0, 144), (F(1, 17), F(2, 3), 144), (F(1, 1000), -2, 256))
)
def test_actual_family_formula_with_independent_fraction_arithmetic(m, Y, L, Q):
    m, Y, L, Q = map(F, (m, Y, L, Q))
    d = outer.enclosure(m, Y, L, Q)
    T = 4 * m * m
    C = 12 * Y / Q
    expected = 8 * C * L * L / (3 * Q * T * T)
    assert d["forward_second_coefficient_upper"] == sp.Rational(expected)
    assert expected == 2 * Y * L * L / (Q * Q * m**4)
    assert d["channel_subtracted_linear_envelope_coefficient"] == sp.Rational(
        4 * C * L * L / (Q * T)
    )
    assert d["strictly_positive_family_coefficient"] == (Y > 0 and L != 0)


@pytest.mark.parametrize("beta", (F(0), F(1, 3), F(3, 5), F(4, 5)))
@pytest.mark.parametrize("T", (16, 144, 4 * 10**400))
def test_parameter_cut_interval_and_density_independently(beta, T):
    T = F(T)
    u = T / (1 - beta * beta)
    lo, hi = (1 - beta) / 2, (1 + beta) / 2
    for x in (lo, hi):
        assert x * (1 - x) == T / (4 * u)
    assert hi - lo == beta
    assert (u - T) * (hi - lo) == u * beta**3


@pytest.mark.parametrize("T", (16, 144, 4 * 10**400))
@pytest.mark.parametrize("multiple", (1, 2, 5, 10**100))
def test_spectral_weight_upper_and_large_mass_lower(T, multiple):
    T, u = F(T), F(T) * multiple
    beta2 = 1 - T / u
    # Square the nonnegative comparisons to avoid irrational rounding.
    assert beta2 >= 0
    weight2 = (u / (u - 1) ** 2) ** 2 * beta2**3
    assert weight2 <= (4 / u) ** 2
    if multiple >= 2:
        assert beta2**3 > F(1, 16)
        assert weight2 > (F(1, 4) / u) ** 2


@pytest.mark.parametrize("u", (16, 36, 10**400))
@pytest.mark.parametrize("x", (F(0), F(1, 4), F(1, 2), F(3, 4), F(1)))
@pytest.mark.parametrize(
    "real,imag", ((-(10**100), 10**200), (0, 0), (2, 1), (4, 10**100))
)
def test_complex_mixed_mass_gap_and_derivative_majorants(u, x, real, imag):
    u = F(u)
    A = x * (1 - x)
    reD = x * u + 1 - x - A * real
    imD = -A * imag
    modulus2 = reD * reD + imD * imD
    assert reD >= x * (u - 4) + (1 - x) > 0
    assert reD - x * (u - 4) - (1 - x) == A * (4 - real) + 4 * x * x
    assert A * A / modulus2 <= (1 - x) ** 2 / (u - 4) ** 2
    assert A**4 / modulus2**2 <= (1 - x) ** 4 / (u - 4) ** 4


@pytest.mark.parametrize("T,Rfactor", ((16, 2), (144, 3), (4 * 10**400, 10**20)))
def test_convergent_tail_moments_from_fractions(T, Rfactor):
    T = F(T)
    R = T * Rfactor
    first = 1 / T - 1 / R
    second = 1 / (2 * T * T) - 1 / (2 * R * R)
    assert 0 < first < 1 / T
    assert 0 < second < 1 / (2 * T * T)


def test_independent_partial_fraction_spectral_primitive():
    u, s, T = sp.symbols("u s T")
    a = (T - s) / (s - 1) ** 2
    b = (T - 1) / (s - 1)
    c = -a
    integrand = (u - T) / ((u - 1) ** 2 * (u - s))
    partial = a / (u - 1) + b / (u - 1) ** 2 + c / (u - s)
    assert sp.factor(partial - integrand) == 0
    primitive = a * sp.log(u - 1) - b / (u - 1) + c * sp.log(u - s)
    assert sp.factor(sp.diff(primitive, u) - integrand) == 0
    assert sp.factor(a + c) == 0


def test_spectral_closed_form_matches_inherited_spacelike_remainder():
    tau, T, d = sp.symbols("tau T d", positive=True)
    s = 1 - d
    # The real positive logarithm ratio at spacelike s fixes the branch.
    raw_argument = (tau - s) / (tau - 1)
    normalized_argument = sp.cancel(1 + d / (tau - 1))
    assert sp.cancel(raw_argument - normalized_argument) == 0
    I = (T - s) / (s - 1) ** 2 * sp.log(sp.cancel(raw_argument)) + (T - 1) / (
        (s - 1) * (tau - 1)
    )
    result = -((s - 1) ** 2) * I
    expected = -(T - 1 + d) * sp.log(normalized_argument) + d * (T - 1) / (tau - 1)
    assert sp.factor(result - expected) == 0


@pytest.mark.parametrize("channel", ("s", "t", "u"))
def test_each_outer_channel_contains_both_internal_positions(channel):
    records = [r for r in ownership.data()["placements"] if r["channel"] == channel]
    assert len(records) == 2
    assert {r["inserted_light_line"] for r in records} == {1, 2}
    assert sum(r["outer_symmetry_weight"] for r in records) == 1


def test_missing_second_internal_position_loses_half_the_answer():
    records = ownership.data()["placements"]
    one_position = [r for r in records if r["inserted_light_line"] == 1]
    assert len(one_position) == 3
    assert sum(r["outer_symmetry_weight"] for r in one_position) == sp.Rational(3, 2)
    assert sum(r["outer_symmetry_weight"] for r in records) == 3


def test_crossed_channel_coefficient_with_independent_formal_series():
    v, b0, b1, b2 = sp.symbols("v b0 b1 b2")
    one = b0 + b1 * v + b2 * v * v / 2
    crossed = one + one.subs(v, -v)
    assert sp.expand(crossed).coeff(v, 1) == 0
    assert sp.expand(crossed).coeff(v, 2) == b2


def test_overall_contact_conversion_cannot_tune_family_second_coefficient():
    v, c, a, b = sp.symbols("v local_contact a b")
    forward = c + a * (2 + v) + a * (2 - v) + b * ((2 + v) ** 2 + (2 - v) ** 2)
    assert sp.diff(forward, v, 2) / 2 == 2 * b
    assert sp.diff(sp.diff(forward, v, 2), c) == 0


def test_formal_correction_is_not_the_resummed_inverse():
    h, D, Ferm = sp.symbols("h D Ferm")
    exact = 1 / (D + h * Ferm)
    formal = 1 / D - h * Ferm / D**2
    assert sp.factor(sp.diff(exact - formal, h).subs(h, 0)) == 0
    assert sp.factor(sp.diff(exact - formal, h, 2).subs(h, 0)) == 2 * Ferm**2 / D**3
    assert "total mass is infinite" in density.data()["total_measure_warning"]


def test_actual_family_bounds_and_limits():
    d = calibration.data()
    e = d["literal_family_enclosure"]["forward_second_coefficient_upper"]
    rel = d["family_second_coefficient_upper_relative_to_positive_tree"]
    assert 0 < e < sp.Rational(1, 10**1418)
    assert 0 < rel < sp.Rational(1, 10**819)
    assert (
        outer.enclosure(2, 0, 1, 144)["strictly_positive_family_coefficient"] is False
    )
    assert (
        outer.enclosure(2, 1, 0, 144)["strictly_positive_family_coefficient"] is False
    )
    assert (
        outer.enclosure(2, 1, -1, 144)["strictly_positive_family_coefficient"] is True
    )


def test_scope_does_not_promote_one_family_to_P8():
    assert "not the exact normalized propagator" in spectral.data()["scope"]
    assert (
        "not a complete new-model two-loop amplitude"
        in outer.enclosure(2, 1, 1, 144)["scope"]
    )
    assert "not included" in ownership.data()["scope"]
    assert audit.controls()["original_P8_not_closed"] is True


def test_exact_counts():
    assert len(audit.residuals()) == 52
    assert len(audit.gates()) == 27
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 54
