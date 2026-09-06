import pytest
import sympy as sp
from p8a_preparation import actual_map


@pytest.mark.parametrize("name", list(actual_map.identities()))
def test_full_equation_identity(name):
    assert actual_map.identities()[name] == 0


@pytest.mark.parametrize("name", list(actual_map.difference_identities()))
def test_full_difference_identity(name):
    assert actual_map.difference_identities()[name] == 0


def test_named_prescription_is_fixed():
    assert actual_map.local_coefficient(2) == -sp.Rational(19, 60)
    assert sp.simplify(actual_map.local_coefficient(3)
                       +sp.Rational(19, 60)+sp.log(sp.Rational(3, 2))/2) == 0


def test_naive_switched_radiation_pressure_is_not_conserved():
    controls = actual_map.negative_controls()
    assert controls["naive_radiation_pressure_conservation_defect"] != 0
    assert controls["free_radiation_past_required_c_at_a_one"] == -sp.Rational(1, 960)
    assert controls["trace_only_constraint_integration_constant"] != 0


def test_zero_source_removes_both_components():
    x = sp.Symbol("x")
    a = sp.Function("a")(x)
    assert actual_map.weighted_source(a, sp.S.Zero, x) == {
        "density": 0, "pressure": 0, "trace": 0}
