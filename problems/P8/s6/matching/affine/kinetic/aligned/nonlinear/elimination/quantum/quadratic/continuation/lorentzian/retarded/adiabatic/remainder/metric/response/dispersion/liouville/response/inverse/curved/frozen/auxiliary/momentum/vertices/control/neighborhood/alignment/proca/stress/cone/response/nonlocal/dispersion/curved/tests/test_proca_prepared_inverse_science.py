"""Independent actual ordinary-Proca prepared tree-plus-loop inverse tests."""
import json

import pytest
import sympy as sp
from p8_affine import verify as original
from p8_proca_prepared_inverse import audit, coordinates, matching, tree, volterra
from p8_proca_rank_one_inverse import kernel


def test_complete_new_coupled_inverse_exact_residuals():
    rows=audit.residuals()
    assert len(rows)==116
    assert sum(len(value) if isinstance(value,sp.MatrixBase) else 1 for value in rows.values())==154
    for value in rows.values():
        assert all(x==0 for x in value) if isinstance(value,sp.MatrixBase) else value==0


def test_analytic_and_scope_gates():
    assert len(audit.gates())==25
    assert all(value is True for value in audit.gates().values())


@pytest.mark.parametrize("sector",("T","L"))
def test_literal_minimal_mode_time_change_checks_canonical_rate_and_frequency(sector):
    rows=coordinates.mode_checks()
    assert rows[sector+"_pure_time_change_frequency"]==0
    assert rows[sector+"_pure_time_change_canonical_rate"]==0
    assert rows[sector+"_exact_covariance_time_change"]==sp.zeros(3)


@pytest.mark.parametrize("sector,component",(("T","energy"),("T","pressure"),("L","energy"),("L","pressure")))
def test_actual_physical_readouts_transform_as_scalars(sector,component):
    assert coordinates.mode_checks()[sector+"_"+component+"_exact_time_change_scalar_readout"]==sp.zeros(1,3)


@pytest.mark.parametrize("sector",("T","L"))
def test_all_adiabatic_subtraction_orders_obey_time_change(sector):
    rows=coordinates.subtraction_and_local_checks()
    for component in ("energy","pressure"):
        for order in range(3):
            assert rows[sector+"_"+component+"_adiabatic_time_change_"+str(order)]==0


@pytest.mark.parametrize("order",(0,1,2))
def test_actual_finite_local_time_variations_and_background_Ward(order):
    rows=coordinates.subtraction_and_local_checks()
    for component in ("energy","pressure"):
        assert rows[component+"_actual_finite_local_time_change_"+str(order)]==0
    assert rows["actual_finite_local_background_Ward_"+str(order)]==0


def test_fixed_profile_contacts_are_retained_in_both_adapted_currents():
    data=coordinates.data()
    eta,scale,u=coordinates.eta,coordinates.scale,coordinates.u
    assert data["fixed_profile_plus_vector_local_time_channel"].coeff(sp.diff(scale,u))==-3*coordinates.S
    assert data["fixed_profile_plus_vector_local_scalar_cross_channel"].coeff(sp.diff(eta,u))==3*coordinates.S
    assert coordinates.checks()["actual_frozen_profile_normalized_stress_contact"]==sp.zeros(2)
    assert data["profile_and_state_reselected"] is False


def test_actual_varied_Ward_and_weighted_cross_adjoint():
    rows=coordinates.checks()
    for key in ("full_fixed_profile_time_channel_from_actual_varied_Ward",
                "local_cross_channels_weighted_adjoint","local_time_profile_density_variation",
                "local_scalar_profile_density_variation"):
        assert rows[key]==0


def test_time_channel_is_local_without_setting_pure_scale_energy_to_zero():
    data=coordinates.data()
    assert data["local_time_channel_maximum_derivative_order"]==2
    assert data["quantum_nonlocal_channel_count"]==1
    assert volterra.checks()["time_constraint_has_no_quantum_fourth_derivative"]==0
    assert volterra.checks()["time_constraint_has_no_quantum_third_scale_derivative"]==0


def test_literal_tree_action_and_matter_charge_reduction_are_replayed():
    rows=tree.checks()
    assert rows["actual_tree_charge_reduction_replayed"]==sp.zeros(2,1)
    assert rows["actual_literal_tree_action_replayed"]==0
    assert rows["independent_adapted_Euler_equals_transformed_actual_tree"]==sp.zeros(2,1)


@pytest.mark.parametrize("row,column,order,expected",(
    (0,0,4,-6*tree.delta**2),(0,1,3,6*tree.delta),
    (1,0,3,-6*tree.delta),(1,1,2,6)))
def test_actual_highest_local_tree_coefficients_in_all_four_entries(row,column,order,expected):
    data=tree.data()
    assert data["complete_local_derivative_coefficients"][str(row)+str(column)][order]==expected


def test_full_derivative_inventory_and_both_cross_terms_are_kept():
    data=tree.data()
    assert data["maximum_local_derivative_orders"]==sp.Matrix([[4,3],[3,2]])
    assert data["complete_local_derivative_coefficients"]["00"][3]!=0
    assert data["highest_time_channel_coefficient"]==-6*tree.delta**2


def test_continuous_clock_map_and_nonzero_time_coefficient_bounds():
    data=tree.data()
    assert data["clock_delta_lower"]==sp.Rational(32,125)
    assert data["clock_delta_upper"]==sp.Rational(1,2)
    assert data["highest_time_channel_absolute_lower"]==sp.Rational(6144,15625)
    assert data["instantaneous_time_channel_inverse_absolute_upper"]==sp.Rational(15625,6144)
    assert tree.checks()["actual_clock_delta_monotone_on_positive_half"]==0


def test_matter_reconstruction_uses_original_negative_charge_sign():
    eta,w,u=coordinates.eta,coordinates.scale,tree.u
    expected=-3*tree.parent.ell*(w+tree.H*eta-tree.delta*sp.diff(eta,u))-tree.parent.w*sp.diff(eta,u)
    assert sp.expand(tree.data()["matter_reconstruction_rate"]-expected)==0


@pytest.mark.parametrize("sector,offset",(("T",3),("L",1)))
def test_new_current_pairs_before_physical_dimension(sector,offset):
    assert (matching.high(sector)["pair"]-sp.Matrix([0,2*(matching.D-offset)])).applyfunc(sp.factor)==sp.zeros(2,1)


def test_first_dimensional_jet_is_not_discarded():
    rows=matching.checks()
    assert rows["new_transverse_high_square_and_first_dimensional_jet"]==sp.zeros(2,1)
    assert rows["new_longitudinal_dimensional_pair_jet_retained"]==0
    assert rows["new_longitudinal_dimensional_square_jet_retained"]==0


@pytest.mark.parametrize("sector",("T","L"))
def test_new_curved_all_momentum_fourth_subtractions_without_freezing_H(sector):
    rows=matching.checks()
    for i in range(2):
        for j in range(2):
            assert rows[sector+f"_new_full_dimensional_all_momentum_curved_fourth_contact_{i}_{j}"]==0


def test_new_scalar_log_and_finite_local_normalizations():
    data=matching.data()
    assert data["normalized_scalar_log_weight"]==4
    assert data["normalized_scalar_fourth_pole"]==2
    assert data["fixed_new_scalar_finite_local_fourth_coefficient"]==-4
    assert matching.checks()["new_scalar_high_kernel_four_primitive_normalization"]==0


def test_actual_dimensional_diagonal_and_connection_identities():
    rows=matching.checks()
    for key in ("generic_current_identity_proper_momentum_scale_cancels_for_every_dimension",
                "generic_current_identity_leading_geometry_diagonal",
                "generic_current_identity_leading_geometry_first_lag_jet",
                "generic_current_identity_generic_subleading_connection_not_dropped"):
        assert rows[key]==0


def test_new_diagonal_inverse_keeps_one_instantaneous_channel():
    assert volterra.checks()["actual_rank_one_loop_completed_by_classical_time_channel"]==sp.eye(2)-sp.eye(2)
    assert kernel.data()["instantaneous_inverse"]==0
    assert volterra.data()["time_channel_inverse_continuous_absolute_upper"]>0


def test_physical_lapse_reconstruction_uses_bounded_local_first_row():
    rows=volterra.checks()
    for order in range(4):
        assert rows["physical_lapse_reconstruction_bounded_first_row_kernel_"+str(order)]==0
    assert rows["prepared_force_time_transform_primitive_sign"]==0


def test_continuous_original_force_and_physical_scale_bounds():
    data=volterra.data()
    assert data["I4_transformed_force_over_original_C0_upper"]==sp.Rational(13,30)
    assert data["I3_transformed_time_force_over_original_C0_upper"]==sp.Rational(47,30)
    assert data["physical_scale_over_adapted_C0_upper"]==sp.Rational(13,5)
    assert data["no_numerical_inverse_norm_or_stability_claim"] is True


def test_gravitational_coupling_formula_keeps_mass_fixed():
    first,second=volterra.coupling(1),volterra.coupling(2)
    assert first["fixed_mass_time_product"]==second["fixed_mass_time_product"]==1000
    assert first["gamma"]==4*second["gamma"]


def test_native_input_guards_after_caches_are_populated():
    matching.high("L")
    volterra.coupling(1)
    assert audit.controls()["rejected_inputs"]==39


def test_exact_json_roundtrip_and_nonfinite_rejection():
    data=original.serialize({"coordinate":coordinates.data(),"tree":tree.data(),"matching":matching.data(),"inverse":volterra.data()})
    assert json.loads(json.dumps(data))==data
    for value in (1.0,sp.Float(1),sp.oo,sp.nan):
        with pytest.raises((TypeError,ValueError)):
            original.serialize(value)
