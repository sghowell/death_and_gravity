import pytest
import sympy as sp
from p8_bimetric_general import background as bg
from p8_bimetric_general import obstruction as ob


@pytest.mark.parametrize("check", [bg.variation_checks, bg.bianchi_checks, ob.branch_checks,
                                  ob.window_checks, ob.control_checks])
def test_exact_arbitrary_beta_equations_and_case_split(check):
    assert all(sp.simplify(value) == 0 for value in check().values())


def test_case_inertias_and_separate_matter_weights_are_positive():
    data = ob.stationary_equations()
    for key in ("dynamic_positive_inertia", "root_positive_inertia", "second_matter_positive_weight"):
        assert data[key].is_positive is True
    assert sp.simplify(data["dynamic_weighted_null"].subs(ob.NFNULL, 0)
                       +2*(bg.MG2+bg.MF2*ob.Y**2)*ob.HDOT+ob.NGNULL) == 0


@pytest.mark.parametrize("parameters", [
    {bg.BETAS[1]: -2, bg.BETAS[2]: 1, bg.BETAS[3]: 0},
    {bg.BETAS[1]: 1, bg.BETAS[2]: -1, bg.BETAS[3]: 1},
    {bg.BETAS[1]: 0, bg.BETAS[2]: 0, bg.BETAS[3]: 0},
])
def test_simple_double_and_everywhere_zero_roots_need_no_division(parameters):
    assert bg.interaction_polynomial(ob.Y).subs(parameters).subs(ob.Y, 1) == 0
    actual = ob.stationary_equations()["dynamic_g_null"].subs(parameters).subs(ob.Y, 1)
    assert sp.simplify(actual-ob.stationary_equations()["root_g_null"]) == 0


def test_forbidden_nec_points_and_omitted_f_equation():
    points = ob.control_points()
    for name in ("forbidden_g_NEC", "forbidden_f_NEC", "regular_ratio_two"):
        assert all(value == 0 for value in points[name].values())
    omitted = points["omitted_f_acceleration"]
    assert omitted == {"EL_Ng": 0, "EL_Nf": 0, "EL_a": 0, "EL_b": sp.Rational(3, 5)}
    assert ob.controls()["forbidden_g_null_stress"] < 0
    assert ob.controls()["forbidden_f_null_stress"] < 0


def test_root_chart_and_formal_negative_lapse_controls_are_not_parent_witnesses():
    data = ob.controls()
    assert data["algebraic_branch_dynamic_derivative_inference_error"] == -sp.Rational(1, 4)
    assert data["zero_lapse_has_zero_metric_determinant"] == 0
    assert data["formal_negative_c_loses_positive_second_matter_weight"] == -1
    assert data["formal_negative_c_can_make_weighted_NEC_sum_negative"] == -4
    assert data["regular_stationary_ratio_two_is_not_y_one"] == 1
