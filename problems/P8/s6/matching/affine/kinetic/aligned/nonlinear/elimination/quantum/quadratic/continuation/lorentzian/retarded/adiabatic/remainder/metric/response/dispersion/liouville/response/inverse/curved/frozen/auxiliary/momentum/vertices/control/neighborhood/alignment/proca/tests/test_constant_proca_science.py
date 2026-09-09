"""Focused scientific checks for the separately named constant-mass candidate."""
import json

import pytest
import sympy as sp
from p8_affine import verify as original
from p8_constant_proca import affine, audit, bounds, model, propagation, quantum
from p8_offclock_scalar import datum, principal, spatial
from p8_zero_source import coefficients as previous_coefficients


def test_all_exact_constant_mass_residuals():
    rows=audit.residuals()
    assert len(rows)==85
    assert sum(len(v) if isinstance(v,sp.MatrixBase) else 1 for v in rows.values())==1330
    for value in rows.values():
        assert all(x==0 for x in value) if isinstance(value,sp.MatrixBase) else value==0


def test_all_domain_and_scope_gates():
    assert len(audit.gates())==23
    assert all(value is True for value in audit.gates().values())


def test_all64_new_affine_stationary_equations():
    rows=affine.checks()
    assert rows["new_all64_stationary_Euler_equations"]==sp.zeros(64,1)
    assert rows["new_stationary_trace_is_same_original_B_dphi"]==sp.zeros(4,1)


def test_connection_Hessian_change_is_literal_and_nonzero():
    d=affine.data()
    assert d["added_full_connection_Hessian"]!=sp.zeros(64)
    assert affine.checks()["literal_mass_counterterm_Hessian"]==sp.zeros(4)
    assert affine.checks()["literal_mass_counterterm_source"]==sp.zeros(64,1)


def test_mass_update_requires_its_source_contact_terms():
    d=affine.data()
    assert d["added_full_connection_source"]!=sp.zeros(64,1)
    wrong=affine.geometry.clean(d["new_full_connection_Hessian"]*d["source_free_stationary_connection"]
                               +d["new_full_connection_source"]-d["added_full_connection_source"])
    assert wrong!=sp.zeros(64,1)


def test_new_projective_and_complement_blocks():
    rows=affine.checks()
    assert rows["new_Hessian_preserves_all_projective_columns"]==sp.zeros(64,4)
    assert rows["new_source_preserves_projective_annihilation"]==sp.zeros(4,1)
    assert rows["retained_map_still_annihilates_all_56_complement_directions"]==sp.zeros(4,60)
    assert rows["new_mass_has_no_complement_retained_cross_block"]==sp.zeros(60,4)


def test_new_inverse_trace_response_and_determinant_lemma():
    rows=affine.checks()
    assert rows["new_full_inverse_trace_lift_Euler_equations"]==sp.zeros(64,4)
    assert rows["new_retained_response_exact_constant_Lorentz_matrix"]==sp.zeros(4)
    assert rows["exact_quotient_determinant_lemma_ratio"]==0
    assert affine.data()["new_retained_mass"]==sp.diag(1,-1,-1,-1)


def test_exact_new_Hamiltonian_and_temporal_solution():
    d=model.data()
    g=model.old.generic()
    assert not d["new_Hamiltonian"].has(g["gt"],g["gs"],g["d"],g["c"])
    assert d["new_temporal_constraint_solution"]==-model.old.j/g["U"]
    assert model.checks()["literal_constant_mass_Hamiltonian_difference"]==0


def test_all_homogeneous_light_data_and_vector_parity_unchanged():
    assert model.checks()["unchanged_all_homogeneous_light_Hamiltonians"]==0
    assert model.checks()["new_complete_Hamiltonian_vector_parity"]==0


def test_complete_quadratic_coefficient_bridge():
    before,after=previous_coefficients.jets(),model.jets()
    assert len(after)==14
    assert sum(any(v!=0 for v in row) for row in after.values())==12
    for key in before:
        assert sp.factor(after[key][0]-before[key][0])==0


@pytest.mark.parametrize("order",range(5))
def test_actual_new_Hamiltonian_mass_coefficient_jets(order):
    checks=model.checks()
    assert checks[f"literal_new_Gauss_coefficient_derivative_{order}"]==0
    assert checks[f"literal_new_spatial_vector_coefficient_derivative_{order}"]==0


def test_vector_squared_lapse_derivative_is_not_falsely_transferred():
    key=model.powers(model.old.VARIABLES.index(model.old.vector),1)
    assert sp.factor(model.jets()[key][1]-previous_coefficients.jets()[key][1])!=0


def test_continuous_new_classical_auxiliary_bounds():
    d=bounds.auxiliary()
    c=d["constant_mass_auxiliary_bounds"]
    assert c["Newton_contraction_upper"]<sp.Rational(1,10)
    assert c["Newton_closed_disc_image_radius_upper"]<sp.Rational(1,10**14)
    assert c["temporal_vector_over_input_radius_upper"]==2
    assert 2*sum(d["same_new_Hamiltonian_piece_triangle_upper"].values())<10000


@pytest.mark.parametrize("radius",(0,sp.Rational(1,10**30),sp.Rational(1,10**24)))
def test_new_temporal_root_radius_domain(radius):
    result=bounds.response_radius(radius)
    assert result["input_radius"]==radius
    assert result["temporal_vector_absolute_upper"]==2*radius


def test_exact_common_scale_tree_comparison():
    d=bounds.tree()
    assert d["common_named_M_tau"]==10**400
    assert d["constant_mass_minus_source_free_cubic_upper"]<sp.Rational(1,10**14)
    assert d["constant_mass_minus_source_free_connected_quartic_upper"]<sp.Rational(1,10**150)
    assert d["constant_mass_minus_source_free_cubic_upper"]==2*d["constant_mass_cubic_block_upper"]
    assert d["constant_mass_minus_source_free_connected_quartic_upper"]==2*d["constant_mass_connected_quartic_tree_upper"]


def test_full_canonical_scale_powers():
    a,b=bounds.at_scale(10**400),bounds.at_scale(2*10**400)
    assert b["constant_mass_minus_source_free_cubic_upper"]==a["constant_mass_minus_source_free_cubic_upper"]/2
    assert b["constant_mass_minus_source_free_connected_quartic_upper"]==a["constant_mass_minus_source_free_connected_quartic_upper"]/4


def test_full_massive_Proca_symbol_and_inverse():
    checks=propagation.checks()
    assert checks["full_massive_Proca_inverse_left"]==sp.zeros(4)
    assert checks["full_massive_Proca_inverse_right"]==sp.zeros(4)
    assert checks["three_massive_polarizations_in_full_determinant"]==0


def test_Proca_divergence_constraint_not_four_unconstrained_modes():
    d=propagation.symbol()
    assert d["full_massive_plane_wave_operator"]!=d["constrained_wave_operator"]
    assert propagation.checks()["full_Proca_operator_implies_divergence_constraint"]==sp.zeros(1,4)
    assert propagation.checks()["unconstrained_principal_Proca_symbol_is_singular"]==0


def test_actual_nearby_vector_canonical_block_and_physical_metric():
    d=propagation.neighborhood()
    assert d["unchanged_zero_scalar_vector_cross_block"]==sp.zeros(3)
    assert propagation.checks()["physical_lapse_and_metric_give_the_correct_central_dispersion"]==0
    assert d["physical_squared_spatial_momentum"]==principal.q/spatial.t**2


@pytest.mark.parametrize("point",(1-datum.RADIUS,1+datum.RADIUS))
def test_positive_new_vector_coefficients_at_actual_offclock_anchors(point):
    d=propagation.neighborhood()
    subs={spatial.t:sp.sqrt(point),principal.q:sp.Integer(10**8)}
    assert d["constant_mass_momentum_Hessian"][2,2].subs(subs)>0
    assert d["constant_mass_coordinate_Hessian"][2,2].subs(subs)>0
    assert d["new_temporal_auxiliary_pivot"].subs(model.old.N,point)<0


def test_same_clock_quantum_modes_do_not_give_same_energy_readout():
    for row in quantum.matrices().values():
        assert row["on_clock_Hamiltonian_difference"]==sp.zeros(2)
        assert row["canonical_energy_change"]!=sp.zeros(2)
        assert row["canonical_pressure_change"]==sp.zeros(2)


def test_both_functional_lapse_operator_jets_really_change():
    for row in quantum.matrices().values():
        assert row["old_minus_new_first_lapse_operator_jet"]!=sp.zeros(2)
        assert row["old_minus_new_second_lapse_operator_jet"]!=sp.zeros(2)


def test_constant_mass_local_energy_and_pressure_matching():
    c=quantum.local_coefficients()
    assert c[0]["energy"]==-sp.Rational(5,2)
    assert c[0]["pressure"]==sp.Rational(5,2)
    for n in range(3):
        assert quantum.checks()[f"actual_constant_mass_finite_local_energy_{n}"]==0
        assert quantum.checks()[f"unchanged_covariantly_matched_local_pressure_{n}"]==0


def test_actual_new_local_Ward_conservation():
    for n in range(3):
        assert quantum.checks()[f"ordinary_local_quantum_Ward_conservation_{n}"]==0


def test_state_and_fixed_old_profile_not_reselected_or_silently_recomputed():
    d=quantum.boundary()
    assert d["state_not_reselected"] is True
    assert d["actual_new_first_and_higher_mass_lapse_jets"]==0
    assert quantum.gates()["old_tadpole_and_frozen_quantum_diagnostic_not_transferred"] is True
    assert quantum.gates()["no_new_quantum_stress_profile_norm_claim"] is True


def test_all_inexact_and_outside_domain_inputs_rejected_after_cache_population():
    bounds.at_scale(10**400)
    bounds.response_radius(0)
    assert audit.controls()["rejected_inputs"]==34


def test_exact_json_roundtrip_and_float_rejection():
    d=original.serialize({"symbol":propagation.symbol(),"mass":affine.data()["retained_mass_change"],
                          "jets":model.jets(),"tree":bounds.tree(),"quantum":quantum.boundary()})
    assert json.loads(json.dumps(d))==d
    for value in (1.0,sp.Float(1),sp.oo,sp.nan):
        with pytest.raises((TypeError,ValueError)):
            original.serialize(value)
