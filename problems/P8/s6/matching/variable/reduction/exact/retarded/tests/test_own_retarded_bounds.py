"""Rational identity, inherited-domain, and physical-clock controls."""

from fractions import Fraction as F

import pytest
from p8_own_retarded import bounds


def test_gap_identity_is_an_all_degree_polynomial_check():
    data = bounds.gap_polynomials()
    assert data["d4_minus_one_minus_vP"] == (0, 0, 0, 0, 0)
    assert data["four_d4_minus_P"] == (0, 10, 20, 15, 4)
    # Independently multiply (1+v) four times using coefficient convolution.
    coefficients = [1]
    for _ in range(4):
        coefficients = [coefficients[0],
                        *(coefficients[i-1] + coefficients[i] for i in range(1, len(coefficients))),
                        coefficients[-1]]
    assert tuple(coefficients) == data["d_fourth"]
    assert tuple(coefficients[1:]) == data["P"]


@pytest.mark.parametrize(("delta", "x"), [
    (F(1, 625), F(-1, 4)), (F(1, 625), F(-1, 8)),
    (F(1, 10**8), F(-1, 4)), (F(1, 1600), F(-1, 17)),
    (F(1, 625), F(0)),
])
def test_exact_positive_delta_ratio_and_domain(delta, x):
    v = delta*x*x
    actual = bounds.denominator_ratio(delta, x)
    literal = (2 + delta - 2/(1 + v)**4)/delta
    assert actual == literal
    assert bounds.denominator_identity(delta, x) == 0
    assert 1 <= actual <= 1 + 8*x*x <= F(3, 2)
    assert 0 <= v <= bounds.V_MAX


def test_literal_zero_delta_is_not_smuggled_in_by_cancellation():
    with pytest.raises(ValueError):
        bounds.denominator_ratio(0, 0)
    with pytest.raises(ValueError):
        bounds.denominator_identity(0, F(-1, 4))
    assert bounds.denominator_ratio(0, F(-1, 4), extension=True) == F(3, 2)
    assert bounds.denominator_ratio(0, 0, extension=True) == 1
    assert bounds.denominator_identity(0, F(-1, 4), extension=True) == 0


@pytest.mark.parametrize("bad", [True, False, 0.0, float("nan"), float("inf"), "1/625", None])
def test_inexact_inputs_are_rejected(bad):
    with pytest.raises(TypeError):
        bounds.denominator_ratio(bad, F(-1, 8))
    with pytest.raises(TypeError):
        bounds.denominator_ratio(F(1, 625), bad)
    with pytest.raises(TypeError):
        bounds.rational(bad)


@pytest.mark.parametrize(("delta", "x"), [
    (F(-1, 1000), F(-1, 8)), (F(1, 500), F(-1, 8)),
    (F(1, 625), F(-1, 3)), (F(1, 625), F(1, 100)),
])
def test_domain_does_not_grow_with_a_convenient_probe(delta, x):
    with pytest.raises(ValueError):
        bounds.denominator_ratio(delta, x)
    with pytest.raises(ValueError):
        bounds.denominator_ratio(delta, x, extension=True)


def test_extension_flag_is_explicitly_boolean():
    for bad in (1, "yes", None):
        with pytest.raises(TypeError):
            bounds.denominator_ratio(0, 0, extension=bad)


def test_inherited_subbox_has_the_required_exact_maximum():
    assert bounds.DELTA_MAX*bounds.WIDTH**2 == F(1, 10_000)
    assert bounds.DELTA_MAX < F(1, 100)
    source = bounds.inherited_enclosures()
    assert source["v"].lo == 0 and source["v"].hi == F(1, 10_000)
    assert source["delta"].lo == 0 and source["delta"].hi == F(1, 625)
    assert source["zeta"].lo == 11 and source["zeta"].hi == 13


def test_coefficient_enclosures_and_independent_cross_multiplied_margins():
    data = bounds.coefficient_box()
    assert (data["k_lower"], data["k_upper"]) == (F(987657, 250000), F(1008683, 250000))
    assert (data["s_lower"], data["s_upper"]) == (F(42355819, 10**6), F(1608, 25))
    # These integer comparisons do not invoke either interval implementation.
    assert 5*987657 > 19*250000
    assert 5*1008683 < 21*250000
    assert 42355819 > 42*10**6
    assert 1608 < 65*25
    assert all(bounds.checks().values())


def test_monotone_ratio_uses_actual_b_cubed_and_full_lapse():
    source = bounds.inherited_enclosures()
    data = bounds.coefficient_box()
    low = source["b_cubed"].lo/source["N"].hi
    high = source["b_cubed"].hi/source["N"].lo
    assert data["k_lower"] <= low <= high <= data["k_upper"]
    assert data["s_lower"] <= 4*source["Q"].lo*source["b"].lo/3
    assert 2*source["Q"].hi*source["b"].hi <= data["s_upper"]
    # Omitting the lapse at the limiting center would give k=8, not k=4.
    assert F(2)**3/F(2) == 4
    assert F(2)**3 > data["k_max"]


def test_physical_clock_and_flux_normalization_at_exact_nonunit_scales():
    M2, tau, epsilon = F(9), F(7), F(1, 40)
    delta = epsilon**2
    k, kx, s = F(4), F(1, 9), F(50)
    Q, q, Qx, Qxx = F(2, 7), F(-1, 3), F(4, 11), F(-5, 13)
    KT = M2*kx/(tau*epsilon)
    QT, QTT = Qx/(tau*epsilon), Qxx/(tau**2*delta)
    nu = M2*s/(tau**2*delta)
    physical = M2*k*QTT + KT*QT + nu*(Q - q)
    normalized = k*Qxx + kx*Qx + s*(Q - q)
    assert physical*tau**2*delta/M2 == normalized
    # Suppressing the time-dependent kinetic derivative changes the equation.
    assert normalized != k*Qxx + s*(Q - q)
    assert kx*Qx != 0


def test_own_f_and_relative_inner_normalizations_remain_distinct():
    for x in (F(-1, 4), F(-1, 8), F(0)):
        gap = 1 + 8*x*x
        k, spring = F(4), 64/gap
        assert spring/k == 16/gap
        assert spring*(1 + 1/k) == 80/gap
        assert spring/k != spring*(1 + 1/k)


def test_scope_and_strict_margins_are_reportable_without_new_authority():
    data = bounds.calibration()
    assert all(value > 0 for value in data["strict_margins"].values())
    assert data["domain_v_identity"] == 0
    assert not data["literal_delta_zero_allowed"]
    assert data["x_zero_is_an_allowed_literal_observation"]
    assert not data["full_parent_background_claim"]
    assert not data["fixed_physical_low_frequency_band_claim"]


def test_polynomial_and_rounding_validation():
    assert bounds.polynomial((4, 6, 4, 1), F(1, 2)) == F(65, 8)
    with pytest.raises(ValueError):
        bounds.polynomial((), 0)
    with pytest.raises(TypeError):
        bounds.polynomial((4, 6.0), 0)
    assert bounds._outward(F(-17, 100), F(-13, 100), 10) == (F(-1, 5), F(-1, 10))
    with pytest.raises(ValueError):
        bounds._outward(2, 1)
    for bad in (True, 10.0, F(10)):
        with pytest.raises(TypeError):
            bounds._outward(0, 1, bad)
    with pytest.raises(ValueError):
        bounds._outward(0, 1, 0)
