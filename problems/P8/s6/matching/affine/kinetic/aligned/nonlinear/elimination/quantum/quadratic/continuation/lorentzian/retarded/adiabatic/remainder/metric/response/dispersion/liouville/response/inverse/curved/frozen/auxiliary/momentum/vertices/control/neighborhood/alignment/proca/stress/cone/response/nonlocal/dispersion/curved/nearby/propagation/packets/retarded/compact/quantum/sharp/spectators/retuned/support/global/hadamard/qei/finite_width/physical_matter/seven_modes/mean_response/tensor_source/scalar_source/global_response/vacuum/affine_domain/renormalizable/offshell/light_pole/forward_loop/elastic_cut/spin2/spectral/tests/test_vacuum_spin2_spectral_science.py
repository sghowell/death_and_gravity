"""Independent exact spectral, normalization, threshold and enclosure tests."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_vacuum_spin2_spectral import audit, calibration, cut, dispersion, moments


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_exact_spectral_cut_or_moment_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_continuous_bound_and_explicit_scope_gate(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "name,call,args",
    calibration.bad_cases(),
    ids=[r[0] for r in calibration.bad_cases()],
)
def test_reject_unsupported_spectral_calibration(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


@pytest.mark.parametrize(
    "name,value",
    [(k, v) for k, v in audit.controls().items() if k != "rejected_inputs"],
    ids=[k for k in audit.controls() if k != "rejected_inputs"],
)
def test_negative_control_and_unclosed_obligation(name, value):
    assert bool(value), name


@pytest.mark.parametrize("argument", (Fraction(3, 2), 2, 10, 10**100))
@pytest.mark.parametrize("terms", (1, 2, 4))
def test_positive_series_enclosure_is_nested_and_keeps_tail(argument, terms):
    d = calibration.q2_enclosure(argument, terms)
    more = calibration.q2_enclosure(argument, terms + 1)
    assert (
        0
        < d["strict_lower"]
        < more["strict_lower"]
        < more["strict_upper"]
        < d["strict_upper"]
    )
    assert d["strict_upper"] - d["strict_lower"] == d["tail_upper"] > 0


@pytest.mark.parametrize("transfer", (4, 5, 6, Fraction(9, 2), sp.Rational(11, 2)))
def test_actual_low_window_density_and_threshold(transfer):
    d = calibration.point(transfer)
    assert d["physical_transfer"] == sp.Rational(transfer)
    assert d["heavy_cut_exact"] == 0
    if transfer == 4:
        assert d["threshold"]
        assert d["light_cut_lower"] == d["light_cut_upper"] == 0
    else:
        assert not d["threshold"]
        assert 0 < d["light_cut_lower"] < d["light_cut_upper"]


def test_validation_not_bypassed_by_equal_cached_float():
    calibration.q2_enclosure(2, 1)
    calibration.point(5)
    for value in (True, 1.0, sp.Float(1)):
        with pytest.raises((TypeError, ValueError)):
            calibration.q2_enclosure(2, value)
    with pytest.raises((TypeError, ValueError)):
        calibration.point(5.0)


@pytest.mark.parametrize("M,T", ((2, 10), (10, 50), (100, 500)))
def test_heavy_unequal_mass_argument_and_threshold_gap(M, T):
    d = cut.data()
    Z = d["heavy_pair_Q2_argument"].subs({d["M"]: M, d["T"]: T})
    assert sp.factor(Z * Z - 1) > 0
    assert Z > 0
    assert T > 4 * M


def test_actual_low_window_is_negligible_but_not_zero():
    d = moments.data()
    assert (
        0
        < d["actual_low_window_slope_moment_lower"]
        < d["actual_low_window_slope_moment_upper"]
    )
    assert d["actual_low_window_slope_moment_upper"] < d["actual_total_slope_lower"]
    assert 0 < d["actual_low_window_fraction_upper"] < sp.Rational(1, 10**196)


def test_full_species_moments_not_replaced_by_low_cut():
    d = moments.data()
    assert d["both_species_fraction_strict_lower"] == sp.Rational(1, 20)
    L, H = d["light_full_slope_moment"], d["heavy_full_slope_moment"]
    assert L.limits[0][1:] == H.limits[0][1:] == (0, 1)
    assert L.function != H.function


def test_fixed_order_heavy_threshold_not_stable_particle_claim():
    assert "zeroth-order" in cut.data()["threshold_prescription"]
    assert "not an exact stable-heavy" in cut.data()["threshold_prescription"]


def test_positive_measure_does_not_assume_finite_total_mass():
    assert "need not have finite total mass" in dispersion.data()["measure_boundary"]
    assert "Not a scattering-amplitude contour" in dispersion.data()["scope"]


def test_Q2_positive_series_coefficients_independently():
    z = sp.Symbol("z", real=True)
    P2 = (3 * z * z - 1) / 2
    for n in range(1, 6):
        coefficient = sp.integrate(P2 * z ** (2 * n), (z, 0, 1))
        assert coefficient == sp.Rational(2 * n, (2 * n + 1) * (2 * n + 3)) > 0


def test_light_threshold_leading_coefficient_from_positive_series():
    e = sp.Symbol("positive_threshold_gap", positive=True)
    M, g = sp.symbols("positive_mass_squared positive_cubic_squared", positive=True)
    leading = (
        g
        / (8 * sp.pi * sp.sqrt((4 + e) * e))
        * sp.Rational(2, 15)
        / (1 + 2 * M / e) ** 3
    )
    assert (
        sp.simplify(
            sp.limit(leading / e ** sp.Rational(5, 2), e, 0, dir="+")
            - g / (960 * sp.pi * M**3)
        )
        == 0
    )


def test_heavy_threshold_leading_coefficient_from_positive_series():
    e = sp.Symbol("positive_threshold_gap", positive=True)
    M = sp.Integer(10)
    T = 4 * M + e
    Z = (T - 2 * M) / sp.sqrt((T - 4) * e)
    leading = (
        e
        / (8 * sp.pi * sp.sqrt(T) * (T - 4) ** sp.Rational(3, 2))
        * sp.Rational(2, 15)
        / Z**3
    )
    assert (
        sp.simplify(
            sp.limit(leading / e ** sp.Rational(5, 2), e, 0, dir="+")
            - 1 / (960 * sp.pi * M ** sp.Rational(7, 2))
        )
        == 0
    )


def test_exact_counts_and_original_scope():
    assert len(audit.residuals()) == 52
    assert len(audit.gates()) == 30
    assert len(audit.controls()) == 10
    assert audit.rejected_inputs() == 47
