"""Independent ordered Gaussian solves and causal variational controls."""
import pytest
import sympy as sp
from p8_aligned_elimination import operators


def test_ordered_operator_and_variational_identities():
    for value in operators.checks().values():
        assert all(item == 0 for item in value) if isinstance(value, sp.MatrixBase) else value == 0
    assert all(value is True for value in operators.proof_checks().values())


@pytest.mark.parametrize("zeta", (sp.Rational(1, 100), sp.Rational(1, 2), 2))
def test_direct_stationary_solve_keeps_both_source_terms(zeta):
    data = operators.dense()
    action = data["action"].subs(data["z"], zeta)
    W = data["W"]
    solved = sp.solve([sp.diff(action, value) for value in W], list(W))
    actual = sp.factor(action.subs(solved, simultaneous=True))
    assert sp.factor(actual-data["effective"].subs(data["z"], zeta)) == 0
    omitted_center = (data["S"].T*data["M"]*data["S"])[0]/2
    assert sp.factor(actual-omitted_center-data["effective"].subs(data["z"], zeta)) != 0


def test_zero_curl_is_exact_auxiliary_control():
    data = operators.dense()
    assert data["solution"].subs(data["z"], 0) == data["S"]
    assert data["effective"].subs(data["z"], 0) == 0


def test_unforced_source_does_not_delete_the_vector_Hessian():
    data = operators.dense()
    zero_source = data["action"].subs(dict.fromkeys(data["S"], 0))
    assert sp.hessian(zero_source, list(data["W"])) == data["K"]
    assert data["K"].det() != 1


def test_retarded_one_copy_variation_contains_an_advanced_piece():
    data = operators.causal_control()
    assert data["wrong_retarded_variation_difference"] == sp.Matrix([data["source"][1], -data["source"][0]])
