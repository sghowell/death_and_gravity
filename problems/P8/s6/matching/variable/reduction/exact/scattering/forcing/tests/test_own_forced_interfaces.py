"""Actual inherited source/domain and independent constant interfaces."""
from fractions import Fraction

import sympy as sp
from p8_own_forced import audit, bridges, loading, phase


def test_eighteen_literal_source_and_independent_constant_interfaces():
    assert len(bridges.identities()) == 18
    assert all(value == 0 for value in bridges.identities().values())


def test_fourteen_actual_branch_domain_premises():
    assert len(bridges.checks()) == 14
    assert all(value is True for value in bridges.checks().values())


def test_strict_limiting_load_uses_a_uniform_nonzero_margin():
    values = audit.constants()
    assert values["sharp_loaded_Y_floor_per_eta"]-Fraction(9, 1280) == Fraction(3989, 30720000)
    assert values["sharp_loaded_Y_u_floor_per_eta"] > loading.calibration()["Y_u_load_lower_over_eta"]


def test_fixed_physical_clock_preserves_original_Q_not_a_chosen_readout():
    k, tau = sp.symbols("k tau", positive=True)
    ku = sp.Symbol("k_u", real=True)
    normalized = phase.endpoint_limit(k, ku)
    physical = normalized*sp.diag(1, tau)
    assert sp.simplify(physical.inv()[0, 0]-sp.sqrt(sp.Rational(loading.A)/k)) == 0
    assert physical.inv()[0, 1] == 0
    assert sp.simplify(physical.det()-tau*k/phase.RHO) == 0


def test_two_errors_and_physical_kinetic_factor_are_both_required():
    values = audit.constants()
    exact = (2-2*Fraction(4, 399))*Fraction(9, 1280)*Fraction(2, 5)
    assert exact == values["Q_pair_separation_per_eta"] == Fraction(237, 42560)
    assert exact-Fraction(1, 200) == Fraction(121, 212800) > 0


def test_no_central_convergence_or_matter_source_premise_is_added():
    assert loading.calibration()["matter_source_response"] is False
    assert phase.calibration()["central_remainder_convergence_assumed"] is False
    assert phase.calibration()["subsequential_reference_limits_asserted"] is False
    assert phase.calibration()["matter_source_or_EFT_verdict"] is False
