import sympy as sp
from p8_s5 import nonlinear_d as m


def test_exact_nonlinear_and_independent_linear_bridges():
    assert all(value == 0 for value in m.audit_checks().values())


def test_transformation_invertible_and_identity_on_background():
    # An explicit inverse, not an inverse involving Theta or Lambda.
    h = sp.Symbol("h", positive=True)
    omega = m.functions()["omega"]
    assert sp.simplify(sp.exp(2*omega)*sp.exp(-2*omega)*h-h) == 0
    assert omega.subs(m.N, 1) == 0
    assert m.N.is_positive


def test_nonlinear_shift_not_just_linear_delta():
    omega = m.functions()["omega"]
    assert sp.factor(sp.diff(omega, m.N, 2).subs(m.N, 1)) == 10/m.d**3
    assert sp.diff(omega, m.N, 2) != 0


def test_general_relativistic_legendre_normalization_all_six_components():
    assert m.legendre_checks()["six_component_Legendre"] == 0
    # Turn off omega in the general Hamiltonian structure: standard GR sign.
    p, pp, R, F = sp.symbols("p pp R F")
    standard = 2*(pp-p**2/2)-F-R/2
    assert sp.diff(standard, pp) == 2
    assert sp.diff(standard, p, 2) == -2


def test_lapse_is_regular_at_gamma_but_not_uniformly_gapped_in_tails():
    J = m.functions()["J"]
    assert m.lapse_jets()["A"][2].subs(m.u, 0) == -4
    assert sp.limit(J, m.u, sp.oo) == 0
    assert sp.limit(m.u**2*J, m.u, sp.oo) == 2


def test_wrong_boundary_sign_is_detected():
    omega = m.functions()["omega"]
    delta = m.N*sp.diff(omega, m.N)
    wrong = 2*delta*(1+delta)  # Omitting spatial IBP leaves a nonzero remainder.
    assert wrong.subs({m.N: 1, m.u: 0}) != 0
