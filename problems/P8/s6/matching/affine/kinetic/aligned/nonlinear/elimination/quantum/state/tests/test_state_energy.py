"""Physical lapse, temporal constraint, pressure and vacuum-pole controls."""
import sympy as sp
from p8_vector_state import energy


def test_actual_energy_constraint_and_canonical_map():
    assert all(value == 0 for value in energy.checks().values())


def test_omitting_clock_mass_variation_does_not_match_the_local_pole():
    d = energy.actual_coefficients()
    assert energy.zero_order_pole()["omitting_clock_variation_pole_difference"] == 20/(9*d["h"])


def test_zero_momentum_canonical_energy_and_pressure_are_regular():
    d = energy.canonical()
    assert sp.factor((d["rhoL"]-d["rhoT"]).subs(energy.q, 0)) == 0
    # Isotropic pressure is obtained by averaging all three homogeneous
    # coordinate directions; a single longitudinal orientation has its
    # own zero-momentum limit but agrees with the transverse formula here.
    assert sp.factor((d["pL"]-d["pT"]).subs(energy.q, 0)) == 0


def test_mass_independent_control_has_ordinary_Proca_energy():
    d = energy.canonical()
    expectedT = ((energy.dv-d["rateT"]*energy.v)**2+d["omega2"]*energy.v**2)/(2*energy.a**3)
    expectedL = ((energy.dv-d["rateL"]*energy.v)**2+d["omega2"]*energy.v**2)/(2*energy.a**3)
    assert sp.factor(d["rhoT"].subs({energy.alpha: 0, energy.beta: 0})-expectedT) == 0
    assert sp.factor(d["rhoL"].subs({energy.alpha: 0, energy.beta: 0})-expectedL) == 0
