from fractions import Fraction

import pytest
import sympy as sp
from p8_variable_prepared import analytic, exact, independent, majorants, physical


def test_every_analytic_majorant_has_independent_fraction_replay():
    primary = majorants.calibration()
    separate = independent.analytic_constants()
    assert primary.keys() == separate.keys()
    assert primary == {key: sp.Rational(value) for key, value in separate.items()}
    assert len(majorants.checks()) == 21


def test_every_physical_constant_has_independent_fraction_replay():
    primary = physical.calibration()
    separate = independent.physical_constants()
    assert primary.keys() == separate.keys()
    assert primary == {key: sp.Rational(value) for key, value in separate.items()}
    result = physical.checks()
    assert all(sp.simplify(value) == 0 for value in result["residuals"].values())
    assert all(value > 0 for value in result["strict_margins"].values())


def test_independent_full_polynomial_replay():
    assert independent.checks()["exact_inverse_polynomials"] == 52


def test_taylor_l1_tail_omitting_K_zero_vanishing_is_invalid():
    # A constant allowed by a disk bound need not satisfy the center tail.
    v, delta = exact.DELTA_RADIUS, exact.DELTA_MAX
    assert 3*v/4 > 300*delta**2
    assert (3*v/4)/v**2 == 300


def test_symplectic_delta_cubed_requires_joint_constant_to_vanish():
    d = exact.DELTA_MAX
    wrong_error = d**2*(3*exact.DELTA_RADIUS)/5
    actual_error = physical.calibration()["omega_error"]
    assert wrong_error > actual_error
    assert actual_error == 3*d**3/5


def test_leading_relative_polynomials_and_uniform_remainder_orders():
    d, jets = analytic.derive(), analytic.low_jets()
    u, delta, momentum = d["u"], d["delta"], d["K"]
    x, eps = jets["x"], sp.Symbol("eps", positive=True)
    for parity, power in (("even", 4), ("odd", 3)):
        leading = (delta+8*u**2)*jets[f"q_{parity}"]
        inner = leading.subs({delta: eps**2, u: eps*x})
        assert sp.expand(inner-eps**power*jets[f"H_{parity}"]) == 0
    assert momentum in jets["q_even"].free_symbols
    sample_delta = sp.Rational(1, 10**12)
    for parity, order in (("odd", 5), ("even", 6)):
        large = physical.relative_jet_remainder(sample_delta, sp.Rational(1, 100), parity)
        small = physical.relative_jet_remainder(sample_delta/4, sp.Rational(1, 200), parity)
        assert large == 2**order*small


def test_center_remainder_is_nontrivial_and_not_locked_coefficient():
    delta = exact.DELTA_MAX
    assert physical.center_remainder(delta) == sp.Rational(48001, 10**18)
    assert physical.center_remainder(delta) < sp.Rational(1, 10000)*delta
    assert sp.Rational(14, 33) != sp.Rational(4, 5)
    assert sp.Rational(62, 165)*delta > 2*physical.center_remainder(delta)


def test_source_and_state_errors_are_distinct_and_keep_amplitudes():
    delta = sp.Rational(1, 10**12)
    assert physical.endpoint_inclusion_error(delta) == 200*delta
    assert physical.fixed_light_error(delta, 2) == 17200*delta
    assert physical.fixed_source_error(delta, 2, 3) == delta*(17200+378000000)
    assert physical.fixed_source_error(delta, 2, 0) == physical.fixed_light_error(delta, 2)


@pytest.mark.parametrize("call", [
    lambda: physical.fixed_light_error(Fraction(1, 10**12), -1),
    lambda: physical.fixed_source_error(Fraction(1, 10**12), 1, True),
    lambda: physical.relative_jet_remainder(Fraction(1, 10**12), Fraction(1, 10), "even"),
    lambda: physical.relative_jet_remainder(Fraction(1, 10**12), 0, "both"),
    lambda: independent.fraction(0.1),
    lambda: independent.inverse_monomial(0, True),
])
def test_negative_physical_and_independent_controls(call):
    with pytest.raises((TypeError, ValueError)):
        call()
