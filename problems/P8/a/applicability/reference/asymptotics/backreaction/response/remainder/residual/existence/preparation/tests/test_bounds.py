"""Independent arithmetic and exclusion controls for the full forced map."""

from fractions import Fraction as Q

import pytest
import sympy as sp
from p8a_preparation import bounds


def test_full_gate_has_strict_exact_margins():
    data = bounds.calibration()
    assert all(value > 0 for value in data["strict_margins"].values())
    assert data["rhs_lipschitz"] < sp.Rational(567, 1000)
    assert data["actual_gate"]["contraction"] < sp.Rational(3, 25)
    assert data["rounded"]["self_map"] == sp.Rational(183, 10**9)
    assert data["rounded"]["ball_margin"] == sp.Rational(817, 10**9)


def test_pair_P_constants_reassembled_with_fraction():
    data = bounds.calibration()
    length, delta = Q(1, 10**10), Q(1, 10**14)
    b, c = Q(2, 10**12), Q(12, 10**9)
    # Integrate the pointwise pair bounds term by term.
    einstein = (9*length**2/2+3*b*length**4/4)/(60*delta)
    curvature = (2*b*length**2/2)/4
    curvature += (Q(1, 4)*length**2/2+(b+Q(1, 2))*length**3/6)/30
    forcing = 8*c*(Q(9, 2)*length**2+9*length**2/2+27*b*length**3/3)
    for key, value in (("Einstein_P", einstein), ("curvature_P", curvature),
                       ("source_P", forcing), ("P", einstein+curvature+forcing)):
        assert str(data["pairs"][key]) == str(value)


def test_full_local_product_rule_all_terms_accounted_for():
    length, radius = sp.symbols("L r", positive=True)
    b, v = sp.symbols("B V", nonnegative=True)
    d_times_x = sp.Rational(17, 30)+radius*length**3/12
    changed_d_times_background = v*length**3/12
    dprime_times_primitive = length/4+radius*length**3/4
    changed_dprime_times_background = b*length**2/4
    expected = (sp.Rational(17, 30)+length/4+b*length**2/4
                +(radius+v)*length**3/12+radius*length**3/4)
    assert sp.expand(d_times_x+changed_d_times_background+dprime_times_primitive
                     +changed_dprime_times_background-expected) == 0
    data = bounds.calibration()
    values = {length: data["length"], radius: data["radius"],
              b: data["geometry"]["u_cap"],
              v: data["geometry"]["background_uprime_cap"]}
    assert expected.subs(values) == data["pairs"]["local_Wick_terms"]


def test_inverse_pole_and_log_normalization_retained():
    length = Q(1, 10**10)
    expected = 2*length/(1-length)+Q(2, 9*10**5)+Q(1, 5)
    assert str(bounds.inverse_bound()) == str(expected)
    assert expected < Q(21, 100)
    # The omitted pole would produce a strictly smaller, invalid exact replay.
    assert expected-(Q(2, 9*10**5)+Q(1, 5)) == 2*length/(1-length)


def test_shared_history_response_coarsening():
    m, duration, length = Q(1, 10**5), Q(3), Q(1, 10**10)
    expected = m*duration*length**2*(Q(5, 4)+Q(32, 2))
    expected += 18*m*m*duration**5*length*2
    data = bounds.calibration()["mode_response"]
    assert str(data["total"]) == str(expected)
    assert expected < Q(1, 10**15)
    assert 2*m*duration**3 < Q(1, 2)


def test_actual_A8_weighted_density_not_an_unproved_y_weight():
    data = bounds.calibration()
    delta = Q(1, 10**14)
    expected = 81*delta*(14000+8+Q(3, 2880*1922))
    assert str(data["A8_inputs"]["source_derived"]) == str(expected)
    assert expected < Q(12, 10**9)
    assert data["A8_inputs"]["S0"] < sp.Rational(12, 10**11)
    assert data["A8_inputs"]["S1"] < sp.Rational(42, 10**10)


def test_source_center_requires_no_cutoff_derivative_bound():
    result = bounds.center_bounds(bounds.LENGTH, bounds.POTENTIAL_CAP, bounds.SOURCE_CAP)
    assert result["P_difference"] == 8*bounds.SOURCE_CAP*(
        3+bounds.LENGTH*(1+9*bounds.POTENTIAL_CAP))
    assert result["rhs"] < sp.Rational(3, 10**7)
    assert bounds.center_bounds(bounds.LENGTH, bounds.POTENTIAL_CAP, 0)["rhs"] == 0


def test_pi_lower_witness_is_an_actual_integral_inequality():
    x = sp.Symbol("x", nonnegative=True)
    polynomial = sum((-1)**j*x**(2*j) for j in range(8))
    assert sp.cancel(1/(1+x*x)-polynomial-x**16/(1+x*x)) == 0
    gap = 4*sp.integrate(polynomial, (x, 0, 1))-3
    assert gap > 0
    assert gap == bounds.elementary_margins()["pi_integral_lower_minus_three"]


def test_initial_plateau_fourth_powers_use_the_actual_clock():
    y, delta = sp.symbols("y delta", positive=True)
    f = 1-delta/y**4
    a = y*f**sp.Rational(1, 4)
    h = f**sp.Rational(1, 4)*sp.diff(a, y)/a
    assert sp.simplify(h**4-y**8/(y**4-delta)**3) == 0
    y0 = Q(5, 2)
    a4 = y0**4-Q(1, 10**14)
    h4 = y0**8/a4**3
    margins = bounds.elementary_margins()
    independent = {"initial_a_fourth_lower": a4-Q(12, 5)**4,
                   "initial_a_fourth_upper": Q(13, 5)**4-a4,
                   "initial_h_fourth_lower": h4-Q(3, 8)**4,
                   "initial_h_fourth_upper": Q(5, 12)**4-h4}
    assert all(value > 0 for value in independent.values())
    assert all(str(margins[key]) == str(value) for key, value in independent.items())


@pytest.mark.parametrize("bad", [True, False, 1.0, sp.Float("0.1"), sp.oo, sp.nan,
                                 "inf", "1/0", sp.Symbol("r")])
def test_nonexact_or_nonfinite_inputs_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        bounds.contraction_gate(bad, 1, 0, 1)


@pytest.mark.parametrize("args", [(1, 1, 0, 1), (1, 2, 0, 1),
                                 ("1/2", 1, 2, 1), (0, 0, 0, 0),
                                 (-1, 0, 0, 1), (1, -1, 0, 1)])
def test_failed_banach_or_domain_gates_rejected(args):
    with pytest.raises(ValueError):
        bounds.contraction_gate(*args)


def test_einstein_term_cannot_be_replaced_by_only_mode_response():
    data = bounds.calibration()
    assert data["pairs"]["Einstein_P"] > 10**8*data["mode_response"]["total"]
    assert data["rhs_lipschitz"] > data["mode_response"]["total"]


def test_interval_enlargement_is_not_automatically_a_contraction():
    pair = bounds.pair_bounds(bounds.DELTA, "1/1000", bounds.RADIUS,
                             bounds.POTENTIAL_CAP, bounds.SOURCE_CAP,
                             bounds.BACKGROUND_DERIVATIVE_CAP,
                             bounds.AUXILIARY_Q_CAP, bounds.AUXILIARY_P_CAP)
    with pytest.raises(ValueError, match="strict full-map contraction"):
        bounds.contraction_gate(bounds.INVERSE_CAP,
                                pair["rhs_without_mode_response"], 0, bounds.RADIUS)
