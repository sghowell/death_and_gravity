from fractions import Fraction

import pytest
import sympy as sp
from p8_variable_constraints import (
    action,
    auxiliary,
    domain,
    observables,
    reduction,
    vector,
)


def test_regular_canonical_one_form_primary_and_secondary_identities():
    assert set(reduction.checks().values()) == {0}


def test_secondary_margin_is_continuous_bernstein_not_sampling():
    d = domain.derive()
    assert d["c1_degree"] == (27,)
    assert d["positive_branch_degree"] == (27, 3)
    assert d["clock_F_identity"] == d["margin_identity"] == 0
    assert domain.checks()["c1"]["count"] == 28
    assert domain.checks()["2_to_4_cleared_polynomial"]["count"] == 112


@pytest.mark.parametrize("time,lapse", [
    (0, 4), (sp.Rational(1, 100), 4), (sp.Rational(-1, 100), 3), (0, 1),
])
def test_all_lapses_shifts_and_canonical_momenta_reconstruct(time, lapse):
    assert set(auxiliary.checks(time, lapse).values()) == {0}


def test_physical_center_Cauchy_map_exact_but_not_a_cone_verdict():
    assert set(observables.center_checks().values()) == {0}
    d = observables.center(4)
    K = action.K
    assert sp.factor(d["det_phase_map"]-5*(K*K-4*K+192)/(2397*K*K)) == 0
    assert d["q_q_bracket"] == sp.zeros(3)
    assert (d["symplectic"].T+d["symplectic"]).applyfunc(sp.cancel) == sp.zeros(6)


def test_actual_time_map_is_not_uniform_order_zero_in_momentum():
    d = observables.center_time_map()
    assert d["det1"] == 0
    assert d["u2_K_coefficient"] == sp.Rational(25, 14382)
    assert d["relative_u2_K_coefficient"] == sp.Rational(5, 6)
    assert d["O_jets"][1] != sp.zeros(3, 6)
    assert d["O_jets"][2] != sp.zeros(3, 6)


def test_nonuniform_free_oscillator_control_reverses_frozen_sign_without_changing_solutions():
    d = observables.nonuniform_controls()
    assert d["exact_transformed_free_oscillator"] == 0
    assert d["original_frozen_frequency_squared"] == -d["transformed_center_frozen_frequency_squared"]
    assert d["map_second_derivative_at_center"] == 2*d["original_frozen_frequency_squared"]


def test_actual_physical_observables_are_not_global_Darboux_configurations():
    d = observables.observable_jets(sp.Rational(1, 100), 4)
    bracket = (d["O"][0]*reduction.J*d["O"][0].T).applyfunc(sp.cancel)
    assert bracket[0, 2] != 0
    assert sp.diff(bracket[0, 2], action.K) == 0
    assert bracket[1, 2] == 0


def test_omitting_canonical_time_boundary_changes_secondary_constraint():
    d, r = action.derive(), reduction.derive()
    actual = reduction.jets(0, 4, 0)["D"][0].subs(action.K, 1)
    correction = sp.diff(r["H_time_correction"], reduction.P0, 2).subs(
        {d["u"]: 0, d["c"]: 4, action.K: 1})
    omitted = sp.factor(actual-correction)
    assert actual == -sp.Rational(5, 24)
    assert omitted == -sp.Rational(37, 144)
    assert omitted-actual == -sp.Rational(7, 144)


def test_literal_vector_action_schur_and_clock_normalization():
    assert set(vector.checks().values()) == {0}
    assert set(vector.center_checks().values()) == {0}


def test_vector_inner_limit_is_fixed_K_and_shrinking_time_window():
    assert set(vector.inner_checks().values()) == {0}


def test_vector_formal_high_frequency_sign_controls():
    d = vector.derive()
    c, u, K = d["c"], d["u"], d["K"]
    assert d["kinetic"].subs({u: 0, c: 4, K: 200}) > 0
    assert d["speed_squared"].subs({u: 0, c: 4}) == sp.Rational(3, 2)
    assert d["kinetic"].subs({u: 0, c: 1, K: 200}) < 0
    assert sp.denom(d["kinetic"].subs({u: 0, c: 1})).subs(K, 192) == 0


@pytest.mark.parametrize("time,lapse,order", [
    (True, 4, 0), (0.0, 4, 0), (0, sp.Float(4), 0), (0, sp.oo, 0),
    (0, 2, 0), (sp.Rational(11, 100), 4, 0), (0, 4, True), (0, 4, 1.0), (0, 4, 3),
    (None, 4, 0), (0, sp.I, 0),
])
def test_exact_domain_and_jet_order_boundary(time, lapse, order):
    with pytest.raises((TypeError, ValueError)):
        reduction.jets(time, lapse, order)


def test_exact_fraction_time_is_accepted():
    assert reduction.jets(Fraction(0), Fraction(4), 0)["D"][0] == -sp.Rational(5, 24)
    assert reduction.jets(0, 4, sp.Integer(1))["D"][0] == -sp.Rational(5, 24)
