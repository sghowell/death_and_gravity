"""No source-formula assumption or independent proof margin is left unattached."""
import sympy as sp
from p8_affine import bridges, connection
from p8_affine import dictionary as d


def test_every_action_and_original_target_interface():
    assert len(bridges.identities()) == 25
    assert all(value == 0 for value in bridges.identities().values())


def test_every_continuous_domain_premise_is_discharged():
    assert len(bridges.checks()) == 14
    assert all(bool(value) is True for value in bridges.checks().values())


def test_closed_endpoint_margin_not_overstated():
    point = d.tube_point(physical_X=sp.Rational(9, 10))
    assert point["additional_quotient_factor"] == sp.Rational(4, 5)
    assert connection.calibration()["CD_exceptional_factor_lower"] == sp.Rational(4, 5)


def test_hidden_extra_kernel_cannot_be_missed_by_displayed_delta():
    singular = connection.exceptional()
    assert singular["quotient_rank"] == 57
    p, s = connection.P, connection.S
    c = connection.cd_substitution()[connection.C]
    quartic = connection.cd_substitution()[connection.F4]
    delta = sp.factor(2*p+c*s**2+2*quartic*s**4)
    assert delta.subs(p, singular["p"]) == 1
