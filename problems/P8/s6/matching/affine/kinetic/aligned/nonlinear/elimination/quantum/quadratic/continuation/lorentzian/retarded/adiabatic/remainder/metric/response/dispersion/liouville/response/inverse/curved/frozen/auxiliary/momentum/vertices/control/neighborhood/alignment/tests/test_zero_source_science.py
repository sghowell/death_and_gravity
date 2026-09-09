"""Scientific new-action checks, without recursive parent replay."""
import json

import pytest
import sympy as sp
from p8_affine import verify as original
from p8_auxiliary_neighborhood import model as old_model
from p8_coupled_vertices import coefficients as old_coefficients
from p8_offclock_scalar import datum, principal, spatial
from p8_zero_source import (
    affine,
    audit,
    bounds,
    coefficients,
    model,
    neighborhood,
    response,
    tree,
)


def test_all_new_action_exact_residuals():
    rows=audit.residuals()
    assert len(rows)==94
    assert sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in rows.values())==406
    for value in rows.values():
        assert all(x==0 for x in value) if isinstance(value,sp.MatrixBase) else value==0


def test_all_new_action_gates():
    assert len(audit.gates())==31
    assert all(value is True for value in audit.gates().values())


def test_literal_source_change_is_not_zero_in_old_action():
    g=old_model.generic()
    before=model.data()["old_Hamiltonian"]
    after=model.data()["new_Hamiltonian"]
    assert sp.factor(after-before)!=0
    assert not after.has(g["d"],g["c"])
    assert model.checks()["literal_normal_mass_source_and_contact_replacement"]==0


def test_actual_Q_second_jet_not_discarded():
    d=model.source_jets()
    h=old_model.coefficients()["background"]["h"]
    assert sp.factor(d["Q_second_lapse_derivative"]+3*sp.diff(h,old_model.u)/h**3)==0
    assert d["second"]!=0
    assert d["zero"]==d["first"]==0


def test_full_covariant_source_is_normal_not_only_homogeneous():
    rows=affine.checks()
    assert rows["actual_original_source_is_covariantly_normal"]==sp.zeros(4,1)
    assert rows["actual_original_source_has_no_spatial_component"]==sp.zeros(4,1)


def test_all64_new_connection_Euler_equations_and_unchanged_Hessian():
    rows=affine.checks()
    assert rows["all64_new_stationary_Euler_equations"]==sp.zeros(64,1)
    assert rows["connection_Hessian_change_is_identically_zero"]==sp.zeros(4)
    assert rows["retained_gradient_lifts_to_all64_connection_sources"]==sp.zeros(64,1)


def test_projective_and_56_complement_annihilation():
    rows=affine.checks()
    assert rows["projective_source_annihilation"]==sp.zeros(4,1)
    assert rows["complement_source_annihilation"]==sp.zeros(60,1)
    assert rows["new_stationary_shifted_one_form_has_zero_normal_source"]==sp.zeros(4,1)


def test_wrong_affine_counterterm_sign_would_fail():
    d=affine.data()
    residual=affine.geometry.clean(d["unchanged_full_connection_Hessian"]*d["new_stationary_connection"]
                                  +d["new_connection_source"]-2*d["added_connection_source"])
    assert residual!=sp.zeros(64,1)


def test_retained_mass_convention_is_not_reversed():
    D=affine.data()["unchanged_retained_mass_matrix"]
    point=D.subs(affine.connection.P,sp.Rational(1,2))
    assert point==sp.diag(1,-1,-1,-1)


def test_new_Hamiltonian_even_vector_parity_and_unchanged_light_restriction():
    new=model.data()["new_Hamiltonian"]
    assert sp.factor(new.subs(old_model.j,-old_model.j)-new)==0
    assert sp.factor((new-model.data()["old_Hamiltonian"]).subs(old_model.j,0))==0
    assert new.has(old_model.vector,old_model.electric,old_model.magnetic)


def test_original_boundary_rows_and_actual_new_gauss_row_retained():
    old=old_coefficients.jets()
    new=coefficients.jets()
    assert len(new)==14
    assert all(len(row)==5 for row in new.values())
    assert sum(any(value!=0 for value in row) for row in new.values())==12
    for powers,row in old.items():
        assert sp.factor(new[powers][0]-row[0])==0
        if sum(powers)==1:
            assert sp.factor(new[powers][1]-row[1])==0


@pytest.mark.parametrize("order",range(5))
def test_actual_new_gauss_coefficient_derivatives_and_positive_root_branch(order):
    assert coefficients.checks()[f"actual_new_Gauss_square_lapse_derivative_{order}"]==0


def test_previous_sharper_bounds_not_assumed_unchanged():
    d=bounds.coefficients()
    assert len(d["new_over_old_bound_ratios"])==5
    assert all(1<value<10**9 for value in d["new_over_old_bound_ratios"])
    assert d["new_universal_Cauchy_lapse_derivative_bounds"]==tuple(
        10000*sp.factorial(n)*100**n for n in range(5))


def test_new_complex_and_Newton_bounds():
    d=bounds.auxiliary()
    assert d["new_Hamiltonian_piece_upper_before_N"]["source_Gauss"]==0
    assert d["new_Hamiltonian_complex_absolute_upper"]<10000
    c=d["new_auxiliary_contraction_bounds"]
    assert c["Newton_contraction_upper"]<sp.Rational(1,10)
    assert c["Newton_closed_disc_image_radius_upper"]<sp.Rational(1,10**14)
    assert c["temporal_vector_over_input_radius_upper"]==4


@pytest.mark.parametrize("radius",(0,sp.Rational(1,10**30),sp.Rational(1,10**24)))
def test_new_temporal_root_radius_endpoint_and_interior(radius):
    d=bounds.response_radius(radius)
    assert d["temporal_vector_absolute_upper"]==4*radius
    assert d["lapse_displacement_upper"]==800000000*radius


def test_positive_physical_majorant_scaling_keeps_stationary_quartic():
    d=bounds.scaling()
    assert d["cubic_maximum_coefficient_ratio_degree"]==4
    assert d["quartic_maximum_coefficient_ratio_degree"]==6
    assert d["all_abstract_cubic_coefficients_nonnegative"] is True
    assert d["all_abstract_quartic_coefficients_nonnegative"] is True
    assert sp.Poly(d["abstract_quartic_positive_majorant"],d["ratio_symbol"]).nth(6)!=0


def test_same_new_named_scale_has_exact_finite_change_bounds():
    d=tree.data()
    assert d["new_named_M_tau"]==10**400
    assert d["new_minus_old_cubic_transition_norm_upper"]<sp.Rational(1,10**14)
    assert d["new_minus_old_connected_quartic_tree_norm_upper"]<sp.Rational(1,10**150)
    assert d["new_cubic_transition_bound"]<sp.Rational(1,1000)
    assert d["new_connected_quartic_tree_bound"]<sp.Rational(1,1000)


def test_common_scale_change_keeps_full_canonical_power_counting():
    center=tree.at_scale(tree.SCALE)
    twice=tree.at_scale(2*tree.SCALE)
    assert twice["new_minus_old_cubic_transition_norm_upper"]==center["new_minus_old_cubic_transition_norm_upper"]/2
    assert twice["new_minus_old_connected_quartic_tree_norm_upper"]==center["new_minus_old_connected_quartic_tree_norm_upper"]/4


def test_wrong_old_sharper_vertex_estimate_is_detectably_different():
    R=bounds.RATIO
    assert R**8>R**6>R**4>1
    old=tree.old.data()
    assert tree.data()["new_cubic_transition_bound"]==R**4*old["cubic_transition_numerator"]/tree.SCALE
    assert tree.data()["new_connected_quartic_tree_bound"]==R**8*old["quartic_connected_tree_numerator"]/tree.SCALE**2


def test_fixed_light_vector_operator_and_covariance_are_identical():
    rows=response.checks()
    assert rows["full_fixed_light_vector_Hessian_unchanged"]==sp.zeros(2)
    assert rows["unchanged_covariance_inverse_equation"]==sp.zeros(2)
    assert rows["noncommuting_Gaussian_completion_keeps_operator_order"]==0


def test_source_functional_and_mean_have_different_light_orders():
    rows=response.checks()
    assert all(rows[f"changed_effective_light_action_has_no_degree_{n}"]==0 for n in range(4))
    assert all(rows[f"changed_coherent_mean_has_no_degree_{n}"]==sp.zeros(2,1) for n in range(2))
    assert response.gates()["source_functional_not_declared_unchanged_at_quartic_light_order"] is True


def test_noncommuting_gaussian_shortcut_is_not_accepted():
    M=sp.Matrix([[2,1],[1,3]])
    K=sp.Matrix([[4,0],[0,5]])
    assert M*K!=K*M
    A=M+K
    assert M*A.inv()*M!=M*M*A.inv()


def test_source_free_phase_correction_cancels_the_old_specific_mixing():
    assert neighborhood.checks()["literal_new_central_phase_correction"]==0
    assert principal.blocks()["momentum_coordinate_block"]!=sp.zeros(3)
    assert neighborhood.data()["new_momentum_coordinate_block"]==sp.zeros(3)
    assert neighborhood.checks()["new_longitudinal_pivot_is_actual_temporal_mass_not_old_joint_kappa"]==0


def test_finite_Euler_time_jets_not_deleted_to_remove_q_fourth_term():
    checks=neighborhood.checks()
    assert checks["arbitrary_remaining_Euler_time_jets_have_no_q_fourth_determinant"]==0
    assert checks["new_bare_coordinate_determinant_has_no_q_fourth_term"]==0
    assert neighborhood.data()["new_Lagrangian_coordinate_Hessian"][0,0]==0


@pytest.mark.parametrize("point",(1-datum.RADIUS,1+datum.RADIUS))
def test_new_temporal_pivot_and_vector_limit_at_exact_offclock_anchors(point):
    d=neighborhood.data()
    assert d["new_joint_temporal_auxiliary_pivot"].subs(old_model.N,point)<0
    assert d["new_longitudinal_high_frequency_pivot"].subs(spatial.t,sp.sqrt(point))>0


def test_mass_ratio_issue_and_old_same_scale_quantum_boundary_not_hidden():
    ratio=neighborhood.data()["new_vector_mass_ratio_on_bounce_slice"]
    assert ratio.subs(old_model.N,1)==1
    assert ratio.subs(old_model.N,1-datum.RADIUS)>1
    assert sp.diff(ratio,old_model.N).subs(old_model.N,1)==-sp.Rational(8,81)
    assert response.gates()["old_same_scale_quadratic_frozen_diagnostic_is_not_cured"] is True
    assert response.gates()["no_numeric_old_pole_bracket_transferred_to_new_scale"] is True


def test_all_exact_input_controls_rejected_even_after_cache_population():
    tree.at_scale(tree.SCALE)
    bounds.response_radius(0)
    assert audit.controls()["rejected_inputs"]==34


def test_json_exact_symbolic_roundtrip_and_float_rejection():
    d=original.serialize({"source":model.source_jets(),"bounds":bounds.coefficients(),
                          "scaling":bounds.scaling(),"tree":tree.data(),"response":response.data()})
    assert json.loads(json.dumps(d))==d
    for value in (1.0,sp.Float(1),sp.oo,sp.nan):
        with pytest.raises((TypeError,ValueError)):
            original.serialize(value)
