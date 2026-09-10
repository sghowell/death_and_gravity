"""Independent interval, scheme, local-reference and omission regressions."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_vacuum_full_one_loop_matching import (
    amplitude,
    audit,
    calibration,
    local,
    logarithm,
    scheme,
)


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_named_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_every_explicit_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[x[0] for x in audit.bad_cases()]
)
def test_invalid_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("terms", range(1, 9))
def test_scale_log_with_independent_fraction_series_and_nesting(terms):
    lows, tails = [], []
    for r in (Fraction(1, 3), Fraction(1, 9)):
        lows.append(
            2 * sum(r ** (2 * j + 1) / Fraction(2 * j + 1) for j in range(terms))
        )
        tails.append(2 * r ** (2 * terms + 1) / (Fraction(2 * terms + 1) * (1 - r * r)))
    low = 400 * (3 * lows[0] + lows[1])
    high = low + 400 * (3 * tails[0] + tails[1])
    d = logarithm.enclosure(terms)
    assert d["scale_log_lower"] == sp.Rational(low)
    assert d["scale_log_upper"] == sp.Rational(high)
    if terms > 1:
        previous = logarithm.enclosure(terms - 1)
        assert previous["scale_log_lower"] < d["scale_log_lower"]
        assert d["scale_log_upper"] < previous["scale_log_upper"]


@pytest.mark.parametrize(
    "T,E",
    ((1, 0), (1, Fraction(1, 4)), (1, 2), (Fraction(1, 10**600), Fraction(1, 10**607))),
)
@pytest.mark.parametrize(
    "lo,hi", ((0, 0), (-2, -1), (Fraction(-1, 100), Fraction(1, 100)))
)
def test_full_coefficient_interval_from_all_fraction_corners(T, E, lo, hi):
    T, E, lo, hi = map(Fraction, (T, E, lo, hi))
    values = [T + old + T * extra for old in (-E, E) for extra in (lo, hi)]
    d = calibration.coefficient_interval(T, E, lo, hi)
    assert d["lower"] == sp.Rational(min(values))
    assert d["upper"] == sp.Rational(max(values))
    assert d["positive_complete_one_loop_coefficient"] == (min(values) > 0)


@pytest.mark.parametrize("M", (3, 36, 10**197))
@pytest.mark.parametrize(
    "x", (Fraction(0), Fraction(1, 4), Fraction(1, 2), Fraction(1))
)
def test_scalar_anchor_parameter_gap_with_independent_fractions(M, x):
    x = Fraction(x)
    Delta = x * M + (1 - x) ** 2
    assert 1 <= Delta <= M
    assert Delta - 1 == x * (M - 2 + x)
    assert M - Delta == (1 - x) * (M - 1 + x)


def test_finite_bare_field_factor_is_essential_not_an_extra_LSZ():
    d = scheme.data()
    h = d["reference_symbols"]["h"]
    r = d["reference_symbols"]["r"]
    fp = d["reference_symbols"]["fp"]
    bare = d["bare_coupling_identities"]
    assert sp.simplify(sp.diff(bare["star_G"] - bare["MS_G"], h).subs(h, 0)) == 0
    # Omitting the finite scalar-plus-fermion field factor leaves this mismatch.
    star_no_finite_field = bare["star_G"].subs({r: 0, fp: 0}, simultaneous=True)
    # Alter the denominator only, retaining the correct converted parameters.
    correct_den = d["canonical_star_field_factor"]
    ms_den = d["MS_field_factor"]
    wrong = bare["star_G"] * correct_den / ms_den
    mismatch = sp.factor(sp.diff(wrong - bare["MS_G"], h).subs(h, 0))
    G = d["parameters"]["G"]
    assert mismatch == (fp - r) * G
    assert star_no_finite_field != bare["star_G"]


def test_MS_total_quartic_pole_and_bare_parameter_pole_are_distinct():
    d = scheme.data()
    L = d["parameters"]["L"]
    total = d["MS_total_vertex_counterterms"]["L"]
    z = d["minimal_scalar_field_UV_increment"]
    assert sp.factor(d["MS_bare_quartic_pole"] - total + 2 * z * L) == 0
    assert sp.factor(d["MS_bare_quartic_pole"] - total) != 0


def test_forward_scale_conversion_from_an_independent_series():
    v = sp.symbols("center_displacement")
    L, g, D = sp.symbols("L g D", positive=True)
    expression = (-L + g / (D - v)) ** 2 + (-L + g / (D + v)) ** 2
    coefficient = sp.expand(sp.series(expression, v, 0, 3).removeO()).coeff(v, 2)
    assert sp.factor(coefficient + 4 * L * g / D**3 - 6 * g * g / D**4) == 0
    assert sp.factor(coefficient / (2 * g / D**3) + 2 * L - 3 * g / D) == 0


def test_tadpole_counterterm_cancels_only_the_stationary_heavy_piece():
    G, T, M, L = sp.symbols("G T M L")
    J = -G * T / 2
    reduced = (L - G * G / M) * T / 2
    counter = -G * J / M
    assert sp.factor(reduced + counter - L * T / 2) == 0
    assert sp.factor(reduced - L * T / 2) != 0


def test_two_sector_inverse_signs_and_canonical_mass_residue():
    v = sp.symbols("s_minus_one")
    r, fp, bs, bf = sp.symbols("r fp b_scalar b_fermion")
    kappa = 1 + r - fp
    raw = -kappa * v + (bf - bs) * v * v
    normalized = raw / kappa
    assert normalized.subs(v, 0) == 0
    assert sp.factor(sp.diff(normalized, v).subs(v, 0)) == -1
    assert sp.expand(raw).coeff(v, 2) == bf - bs


def test_actual_full_error_pole_and_local_reference_bounds():
    d = calibration.data()
    p = d["actual_MS_interaction_parameters"]
    T = 4 * p["lambda"]
    band = d["complete_one_loop_canonical_coefficient_interval"]
    assert (
        T * (1 - sp.Rational(1, 10**6))
        < band["lower"]
        < band["upper"]
        < T * (1 + sp.Rational(1, 10**6))
    )
    assert d["complete_one_loop_relative_error_upper"] < sp.Rational(1, 10**6)
    assert d["complete_Phi_pole_remainder_upper"] < sp.Rational(2, 10**405)
    assert d["positive_local_Phi_curvature_lower"] > 0
    assert d["complete_finite_mass_reference_interval"][1] < 0
    assert 10**799 < d["complete_one_loop_vacuum_energy_interval"][0]
    assert d["two_scalar_vacuum_constant_absolute_upper"] > 0


def test_valid_but_inconclusive_error_is_not_a_pass():
    d = calibration.coefficient_interval(1, 2, 0, 0)
    assert d["lower"] == -1 and d["upper"] == 3
    assert d["positive_complete_one_loop_coefficient"] is False


def test_scope_does_not_promote_one_loop_to_P8():
    assert "two-loop and later errors" in amplitude.data()["scope"]
    assert "No global quantum potential" in local.data()["scope"]
    assert (
        "not a higher-loop error bound"
        in calibration.coefficient_interval(1, 0, 0, 0)["scope"]
    )
    assert audit.controls()["original_P8_not_closed"] is True


def test_exact_counts():
    assert len(audit.residuals()) == 57
    assert len(audit.gates()) == 36
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 61
