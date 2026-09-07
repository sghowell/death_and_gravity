from fractions import Fraction as Q

import pytest
import sympy as sp
from p8_trimetric import extension, independent, matter, model, vacuum


@pytest.mark.parametrize("module", (model, matter, vacuum, extension))
def test_exact_source_and_vacuum_identities(module):
    assert all(sp.simplify(value) == 0 for value in module.checks().values())


def test_independent_fraction_matrix_jets_and_polynomial_engine():
    report = independent.checks()
    assert report["full_matrix_first_jet_directions_checked"] == 64
    assert len(report["coefficientwise_identities"]) == 3


def test_exact_matrix_jet_ring_rejects_inexact_inputs_and_singular_inverses():
    assert independent.Dual(Q(2), Q(1))/independent.Dual(Q(3), Q(2)) == independent.Dual(Q(2, 3), Q(-1, 9))
    for value in (True, 0.1, float("inf"), sp.Float("0.1")):
        with pytest.raises(TypeError):
            independent.Dual(value)
    with pytest.raises(ZeroDivisionError):
        independent.Dual(1)/independent.Dual(0, 1)
    with pytest.raises(ZeroDivisionError):
        independent.inverse([[0]*4 for _ in range(4)])
    assert independent.inverse([[int(i == j) for j in range(4)] for i in range(4)]) == independent.identity()


def test_dropping_full_metric_equations_gives_false_original_vacuum():
    controls = vacuum.controls()
    assert controls["dropping_e_f_equations_false_vacuum_u_Euler_00"] == 0
    assert controls["original_regular_g_flat_Euler_00"] == -2
    assert controls["original_regular_f_flat_Euler_00"] == -2
    assert controls["only_one_nonzero_link_still_fails_g_Euler_00"] == -2
    assert controls["B_zero_still_fails_g_Euler_00"] == -2


def test_no_flat_obstruction_does_not_divide_B_or_use_the_u_source_equation():
    u = sp.diag(2, 3, 4, 5)
    out = model.euler_maps(sp.eye(4), sp.eye(4), u, b=0, pg=0, pf=sp.Rational(-2),
                           matter_gradient=sp.ones(4))
    assert out["E_e"] == sp.zeros(4)
    assert u.T*out["E_v"] == 4*u.det()*sp.eye(4)


def test_both_zero_links_are_outside_the_theorem_and_lose_the_auxiliary_inverse():
    u = sp.diag(*sp.symbols("u0:4", positive=True))
    out = model.euler_maps(sp.eye(4), sp.eye(4), u, b=0, pg=0, pf=0, epsilon=0)
    assert out["E_e"] == out["E_v"] == out["E_u"] == sp.zeros(4)
    assert vacuum.controls()["both_links_zero_B_zero_flat_equations"] == 0


def test_chosen_parent_has_no_beta2_only_effective_potential():
    d = vacuum.derive()
    betas = d["effective_beta_values"]
    assert sp.cancel(betas[1]*betas[3]-betas[2]**2) == 0
    assert vacuum.controls()["pure_beta2_middle_minor"] == -1
    assert vacuum.controls()["old_beta2_with_vacuum_shift_middle_minor"] == -sp.Rational(1, 4)
    assert sp.cancel(d["effective_b"]/d["effective_a"]-d["p_f"]/d["p_g"]) == 0
    assert d["source_paper_a"] == -d["effective_a"]


def test_constant_matter_potential_does_not_repair_full_flat_e_f_equations():
    val = sp.Symbol("V", real=True)
    u = sp.eye(4)
    source = matter.canonical_source(u, sp.zeros(4, 1), val)
    out = model.euler_maps(u, u, u, matter_gradient=source["J_u"], b=-6, pg=1, pf=1)
    assert out["E_e"] == out["E_v"] == -2*sp.eye(4)
    assert out["E_u"] == -model.EPS*val*sp.eye(4)


def test_formal_vacuum_cancellation_does_not_cancel_exact_action():
    d = vacuum.derive()
    x = {d["B"]: -6, d["epsilon"]: 1, d["V"]: -4}
    assert d["truncated_constant_V_action_over_h0_volume"].subs(x) == 0
    assert d["exact_constant_V_action_over_h0_volume"].subs(x) == -sp.Rational(27, 16)
    assert d["B_with_constant_matter_potential"].subs(x) == -8


def test_physical_matter_metric_and_sourced_lorentz_constraint_are_not_frozen():
    c = matter.controls()
    assert c["leading_e_v_antisymmetric_01"] == 0
    assert c["positive_canonical_timelike_Y"] == sp.Rational(5, 9)
    assert c["required_e_v_antisymmetric_epsilon_coefficient_01"] == -sp.Rational(1, 6)
    assert c["actual_h_correction_01"] == -sp.Rational(1, 2)


def test_formal_auxiliary_expansion_rejects_B_zero_and_wrong_dimensions():
    with pytest.raises(ValueError):
        matter.first_order(sp.eye(4), sp.eye(4), sp.zeros(4, 1), b=0)
    with pytest.raises(ValueError):
        matter.canonical_source(sp.eye(4), sp.zeros(3, 1))
    with pytest.raises(ValueError):
        model.euler_maps(sp.eye(3), sp.eye(4), sp.eye(4))


def test_beta4_extension_is_a_real_flat_countercontrol_with_positive_mass():
    d = extension.derive()
    assert d["relative_mass_squared"].is_positive
    assert extension.controls()["positive_vacuum_mass_squared_G_F_q_1"] == 2
    assert extension.controls()["omitting_beta4g_restores_nonzero_g_Euler_00"] == -2
    assert extension.controls()["negative_link_changes_mass_sign"] == -2


def test_extended_quartic_potential_is_not_only_an_fp_fit():
    d = extension.derive()
    delta = sp.Symbol("delta", real=True)
    # A traceful, non-TT direction sees the exact quartic as well.
    e, v, u = (1+delta)*sp.eye(4), (1-delta)*sp.eye(4), sp.eye(4)
    action = model.potential_density(e, v, u, b=d["B"], pg=d["p_g"], pf=d["p_f"], bg=d["beta4g"], bf=d["beta4f"])
    assert sp.expand(action) == 24*d["q"]*delta**2+4*d["q"]*delta**4


def test_unequal_planck_masses_do_not_identify_physical_source_with_massless_mode():
    d = extension.derive()
    massive = sp.Matrix([1/d["G"], -1/d["F"]])
    projection = (d["physical_source_covector"].T*massive)[0]
    assert sp.factor(projection) == (d["F"]-d["G"])/(2*d["F"]*d["G"])
    assert projection.subs(d["F"], d["G"]) == 0
