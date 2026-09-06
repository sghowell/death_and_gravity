import pytest
import sympy as sp
from p8a_existence import log_inverse as inv


def test_all_logarithmic_identities():
    assert all(sp.simplify(value) == 0 for value in inv.identities().values())


def test_dyadic_positive_pole_and_cut_are_all_retained():
    data = inv.dyadic_norm_bound(1, 512)
    assert data["length"] == sp.Rational(1, 2**1024)
    assert data["split"] == 2**512
    assert data["norm_upper"] == sum(data[key] for key in (
        "pole_norm_upper", "low_cut_norm_upper", "high_cut_norm_upper"))
    assert 0 < data["norm_upper"] < sp.Rational(1, 3)
    assert inv.calibration()["abstract_conditional_example"]["distance_upper"] == sp.Rational(3, 20)


@pytest.mark.parametrize("value", [True, False, 0.1, float("inf"), sp.oo, -1, sp.sqrt(2), sp.Symbol("x")])
def test_inexact_or_invalid_inputs_rejected(value):
    with pytest.raises((TypeError, ValueError), match="exact nonnegative rational"):
        inv.exact_nonnegative(value)


@pytest.mark.parametrize("b,n", [(0, 0), (1, 2), (10, 20), (sp.Rational(1, 2), 512), (0, sp.Rational(3, 2))])
def test_dyadic_domain_is_explicit(b, n):
    with pytest.raises(ValueError):
        inv.dyadic_norm_bound(b, n)


@pytest.mark.parametrize("args", [(1, 1, 0, 1), (2, 1, 0, 1), (sp.Rational(1, 2), 1, 1, 1), (0, 0, 0, 0)])
def test_failed_conditional_contractions_are_rejected(args):
    with pytest.raises(ValueError):
        inv.contraction_conditions(*args)


def test_zero_functional_case_and_closed_ball_equality():
    result = inv.contraction_conditions(1, 0, 1, 1)
    assert result["q_upper"] == 0
    assert result["self_map_margin"] == 0
    assert result["distance_upper"] == 1


def test_parameter_floor_monotonicity_and_shorter_interval():
    assert inv.dyadic_norm_bound(2, 512)["norm_upper"] > inv.dyadic_norm_bound(1, 512)["norm_upper"]
    assert inv.dyadic_norm_bound(1, 256)["norm_upper"] > inv.dyadic_norm_bound(1, 512)["norm_upper"]


def test_controls_do_not_claim_smoothness_or_discard_the_pole():
    controls = inv.controls()
    assert controls["omitted_pole_Laplace_defect"] != 0
    assert controls["wrong_signed_cut_jump"] != 0
    assert "not_C1" in " ".join(controls)
