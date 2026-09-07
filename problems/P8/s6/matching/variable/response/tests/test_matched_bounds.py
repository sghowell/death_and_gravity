"""Uniform rational/interval enclosures, strict domains and error semantics."""

from fractions import Fraction as Q

import pytest
import sympy as sp
from p8_variable_response import bounds, exact, independent


def test_exact_continuous_mass_and_coefficient_margins():
    d = bounds.coefficient_checks()
    assert all(sp.simplify(value) == 0 for value in d["residuals"].values())
    assert all(value > 0 for value in d["strict_continuous_margins"].values())


def test_independent_entire_box_derivatives_are_not_point_samples():
    d = independent.analytic_boxes()
    assert d["closed_u_box"] == (Q(-1, 10), Q(1, 10))
    assert d["closed_c_box"] == (2, Q(201, 100))
    assert d["u_subintervals"] == 40
    for key in ("A", "C", "d_cross", "f_cross", "E", "b_analytic"):
        assert d["analytic_c_derivative_abs_upper"][key] < 8
    for key in ("j_light", "j_heavy"):
        assert d["analytic_c_derivative_abs_upper"][key] < 1
        assert d["physical_source_projection_abs_upper"][key] < 1


def test_all_independent_norm_constants_bridge():
    primary = bounds.calibration()
    for key, value in independent.norm_constants().items():
        assert primary[key] == sp.Rational(value)
    bounds.checks()


def test_nontrivial_example_and_vacuous_upper_box_are_distinct():
    assert bounds.transfer_error(sp.Rational(1, 10**21)) == sp.Rational(1, 250)
    assert bounds.transfer_error(exact.DELTA_MAX) == 40
    assert bounds.prepared_error(sp.Rational(1, 10**21), 1, 0) == sp.Rational(1, 250)


def test_prepared_leakage_upper_bound_is_not_sharp_delta_control():
    small, smaller = sp.Rational(1, 10**12), sp.Rational(1, 10**21)
    assert bounds.transfer_error(smaller) < bounds.transfer_error(small)
    assert bounds.transfer_error(smaller)/smaller > bounds.transfer_error(small)/small


def test_outer_volterra_tail_is_quantitative_and_decreases():
    r = exact.SLICE
    assert bounds.outer_picard_tail(r, 0) == 2*64*r
    assert bounds.outer_picard_tail(r, 12) < bounds.outer_picard_tail(r, 6)
    assert bounds.outer_picard_tail(2*r, 0) == 4*128*r


def test_outward_interval_root_and_mixed_derivative():
    interval, jet = independent.Interval, independent.Jet
    root = interval(Q(2, 3), Q(7, 3)).sqrt()
    assert root.lo**2 <= Q(2, 3) <= root.hi**2
    assert root.lo**2 <= Q(7, 3) <= root.hi**2
    u, c = jet(2, variable=(1, 0)), jet(3, variable=(0, 1))
    mixed = ((u*c)**2).get((1, 1))
    assert mixed.pair() == (24, 24)


def test_elementary_exponential_not_a_float_comparison():
    assert bounds.exponential_majorant(sp.Rational(2, 3)) == sp.Rational(41, 21)
    assert bounds.exponential_majorant(0) == 1
    assert independent.checks()["exp9_lower_taylor"] > 2002


@pytest.mark.parametrize("bad", [True, 1.0, "1", sp.Float(1), sp.oo, sp.nan, sp.Symbol("x", positive=True)])
def test_exact_public_number_rejection(bad):
    with pytest.raises((TypeError, ValueError)):
        bounds.prepared_error(sp.Rational(1, 10**21), bad, 1)


@pytest.mark.parametrize("bad", [True, 0.1, sp.Rational(1, 2), "1/2"])
def test_independent_fraction_domain_rejection(bad):
    with pytest.raises((TypeError, ValueError)):
        independent.Interval(bad)


def test_box_and_order_negative_controls():
    bad_calls = [lambda: independent.Interval(2, 1),
                 lambda: independent.Interval(-1, 1).inverse(),
                 lambda: independent.Interval(-1).sqrt(),
                 lambda: independent.Interval(2)**Q(1, 3),
                 lambda: bounds.outer_picard_tail(sp.Rational(3, 100), 1),
                 lambda: bounds.outer_picard_tail(0, True),
                 lambda: bounds.outer_picard_tail(0, -1),
                 lambda: bounds.exponential_majorant(1),
                 lambda: bounds.source_error(sp.Rational(1, 10**12), -1),
                 lambda: exact.side(1.0)]
    for call in bad_calls:
        with pytest.raises((TypeError, ValueError)):
            call()
