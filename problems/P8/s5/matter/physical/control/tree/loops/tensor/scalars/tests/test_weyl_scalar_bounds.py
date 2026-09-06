import pytest
import sympy as sp
from p8_m1_weyl_scalar import bounds as b


def test_exact_laurent_rows_and_every_constant_margin():
    assert all(value == 0 for value in b.row_checks().values())
    assert b.row_majorants() == {"unitary": {"coordinate_row_over_q": 138, "momentum_row": 172},
                                 "gamma": {"coordinate_row_over_q": 432, "momentum_row": 25}}
    assert all(value >= 0 for value in b.validate_constants().values())


def test_large_scale_conditional_example_and_zero():
    estimate = b.estimates(1, 10**324)
    assert estimate["epsilon_upper"] == sp.Rational(10)**(-648)
    assert estimate["interaction_exponent_upper"] == sp.Rational(10)**(-534)
    assert estimate["relative_same_reduced_data_evolution_error"] == sp.Rational(10)**(-533)
    assert estimate["relative_first_Born_remainder"] == sp.Rational(10)**(-1067)
    assert estimate["branch_phase_map_norm"] == sp.Rational(10)**(-564)
    assert estimate["relative_same_pulled_old_phase_error"] == sp.Rational(10)**(-532)
    assert b.estimates(0, 1)["chosen_branch_phase_inverse_norm"] == 1


def test_scalar_gate_does_not_borrow_the_tensor_gate():
    with pytest.raises(ValueError, match="10\\^-116"):
        b.estimates(1, 10**15)
    assert b.estimates(1, 10**58)["interaction_exponent_upper"] == sp.Rational(1, 100)
    controls = b.controls()
    assert controls.pop("zero_matching_coefficient_has_zero_bound") == 0
    assert all(value > 0 for value in controls.values())


@pytest.mark.parametrize("coefficient,scale", [(1.0, 10**100), (1, 1e100), (sp.oo, 1),
                                               (sp.nan, 1), (1, 0), (-1, 10**100),
                                               (sp.Symbol("c"), 10**100)])
def test_unsafe_bound_inputs_fail_closed(coefficient, scale):
    with pytest.raises(ValueError):
        b.estimates(coefficient, scale)


def test_generic_energy_norm_rank_one_inequality_by_exact_squares():
    # Independent anisotropic, mixed-channel direct check of the rank-one
    # singular-value expression. No eigensystem assumption is used in proof.
    q = sp.Integer(101)
    w = sp.diag(71, 137)
    ry, rp = sp.Matrix([5, -7]), sp.Matrix([11, 13])
    dual_squared = (ry.T*w.inv()*ry+rp.T*rp)[0]
    image_squared = (rp.T*w*rp+ry.T*ry)[0]
    assert image_squared <= sp.Rational(3, 2)*q*dual_squared
    exact_generator_norm_squared = (sp.Rational(8, 3)*q**2)**2*dual_squared*image_squared
    assert exact_generator_norm_squared <= (sp.Rational(8, 3)*q**2)**2*sp.Rational(3, 2)*q*dual_squared**2


def test_neumann_and_time_ordered_remainder_enclosures():
    mu, eta = sp.Rational(1, 100), sp.Rational(1, 10**32)
    exponential_upper = 1/(1-mu)
    same_pulled_error = b.G_NORM*((exponential_upper-1)+eta*(exponential_upper+1))/(1-eta)
    assert same_pulled_error < 6*mu+8*eta
    assert b.G_NORM*mu/(1-mu) < b.EVOLUTION_ERROR*b.EPS_MAX
    assert b.G_NORM*mu**2/(2*(1-mu)) < b.BORN_REMAINDER*b.EPS_MAX**2
