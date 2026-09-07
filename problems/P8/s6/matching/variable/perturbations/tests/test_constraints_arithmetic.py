import pytest
import sympy as sp
from p8_variable_constraints import action, domain, independent, reduction, vector


def test_independent_polynomial_derivation_matches_all_Bernstein_coefficients():
    primary = domain.derive()
    other = independent.margin_coefficients()
    assert primary["c1_coefficients"] == other["c1"]
    assert primary["positive_branch_coefficients"] == other["positive"]


@pytest.mark.parametrize("fixture", independent.scalar_fixtures())
def test_independent_quadratic_polarization_replays_all_phase_jet_coefficients(fixture):
    d = reduction.jets(fixture["u"], fixture["c"], 1)
    point = {action.K: fixture["K"]}
    for prefix in ("M", "H"):
        for order in range(2):
            assert d[prefix][order].subs(point) == sp.Matrix(fixture[f"{prefix}{order}"])
    for order in range(2):
        assert d["p0"][order].subs(point) == sp.Matrix([fixture["p0"][order]])
        assert d["D"][order].subs(point) == fixture["D"][order]


@pytest.mark.parametrize("fixture", independent.vector_fixtures())
def test_independent_literal_shift_solve_and_Taylor_pump(fixture):
    d = vector.derive()
    point = {d["u"]: fixture["u"], d["c"]: fixture["c"], d["K"]: fixture["K"]}
    for key, primary in (("kinetic", "kinetic"), ("bare", "omega_bare_squared"),
                         ("pump", "normalization_pump"), ("speed_squared", "speed_squared"),
                         ("mass_squared", "mass_algebraic_squared")):
        assert sp.cancel(d[primary].subs(point)) == fixture[key]


def test_independent_matrix_and_fixture_input_controls():
    with pytest.raises(ValueError):
        independent.inverse([[1, 1], [1, 1]])
    with pytest.raises(TypeError):
        independent.inverse([[1.0]])
    with pytest.raises(ValueError):
        independent.scalar_jet_fixture(0, 4, 0)
    with pytest.raises(TypeError):
        independent.scalar_jet_fixture(True, 4, 1)
