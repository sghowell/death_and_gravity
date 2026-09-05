import sympy as sp
from p8_matching import horndeski
from p8_s5 import nonlinear_d as model


def test_independent_ADM_and_tensor_map_derivation():
    result = horndeski.identities()
    assert all(sp.simplify(value) == 0 for value in result["residuals"].values())
    assert result["root_bracket_signs"] == [-1, 6]


def test_exact_two_jacobian_crossings_are_inside_finite_time():
    result = horndeski.identities()
    root = result["positive_crossing_time_over_tau"]
    assert root.is_positive
    assert (1-root).is_positive
    lam = result["clock_J_over_C"]
    assert lam.subs(model.u, 0) == -1
    assert lam.subs(model.u, 1) == sp.Rational(3, 4)
    assert sp.simplify(lam.subs(model.u, root)) == 0
    assert sp.simplify(lam.subs(model.u, -root)) == 0


def test_regular_background_metric_is_not_regular_field_map():
    result = horndeski.derive()
    assert result["exact_residuals"]["normalized_C"] == "0"
    assert result["exact_residuals"]["normalized_W"] == "0"
    assert result["exact_residuals"]["positive_root"] == "0"
    assert "original witness" in result["scope"]


def test_s5_lapse_removal_was_not_a_horndeski_map():
    # S5 changes only the spatial metric, not the lapse: W=1. Its spatial
    # transformation is regular (J=1), but fails the additional Horndeski condition.
    X = sp.Symbol("X", positive=True)
    C = horndeski.identities()["normalized_C_example"]
    missing_condition = 1-C-X*sp.diff(C, X)
    assert sp.simplify(missing_condition.subs({X: 1, model.u: 0})-2) == 0
