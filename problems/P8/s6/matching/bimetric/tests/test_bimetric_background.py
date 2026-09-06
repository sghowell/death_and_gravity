import pytest
import sympy as sp
from p8_bimetric import background as bg


@pytest.mark.parametrize("check", [bg.variation_checks, bg.bianchi_checks, bg.bounce_checks,
                                  bg.NEC_violating_point_checks])
def test_exact_two_lapse_and_bounce_identities(check):
    assert all(sp.simplify(value) == 0 for value in check().values())


def test_obstruction_has_no_mass_gap_escape():
    data = bg.bounce_equations()
    assert data["CD_required_null_stress"].is_negative
    assert data["minimum_Hdot_error"].is_positive
    assert sp.simplify(bg.TAU**2*data["minimum_Hdot_error"]-4) == 0
    assert sp.simplify(data["summed_null"].subs(bg.RHOP, 0)+2*(bg.MG2+bg.MF2)*bg.HDOT) == 0
    assert data["CD_half_tau_window_endpoint_H_magnitude"] == sp.Rational(8, 5)/bg.TAU
    assert data["CD_half_tau_window_minimum_Hdot"] == sp.Rational(48, 25)/bg.TAU**2


def test_positive_branch_constraint_and_regular_minkowski_control():
    data = bg.bounce_equations()
    assert data["f_constraint_at_stationarity"].subs(bg.Y, 1) == 0
    assert sp.factor(data["f_constraint_at_stationarity"]) == -bg.NU*(bg.Y-1)*(bg.Y**2+bg.Y+1)/bg.Y**3
    assert data["f_lapse_from_acceleration"].subs(bg.HDOT, 0) == 1
    assert data["summed_null"].subs({bg.HDOT: 0, bg.RHOP: 0}) == 0


def test_omission_and_forbidden_matter_controls():
    controls = bg.controls()
    assert controls["discarding_f_lapse_loses_ratio_constraint"] == -7*bg.NU/8
    assert controls["discarding_f_acceleration_loses_Planck_term"] == 2*bg.MF2*bg.HDOT
    assert controls["forbidden_NEC_control_null_stress"] == -sp.Rational(2, 5)
    assert controls["zero_interaction_loses_branch_factor"] == 0


def test_kinetic_nec_independent_of_potential():
    data = bg.equations()
    assert sp.simplify(data["rho"]+data["pressure"]-bg.KIN/bg.Ng**2) == 0
    assert not sp.simplify(data["rho"]+data["pressure"]).has(bg.POT)
