"""Independent central acceleration, recentered root and local evolution tests."""
import json

import pytest
import sympy as sp
from p8_affine import verify as original
from p8_proca_nearby_bounce import audit, existence, jets, taylor


def test_every_literal_identity_and_scalar_count():
    rows=audit.residuals()
    assert len(rows)==45
    assert all(not isinstance(v,sp.MatrixBase) for v in rows.values())
    assert all(v==0 for v in rows.values())


def test_all_continuous_and_analytic_proof_gates():
    assert len(audit.gates())==45
    assert all(v is True for v in audit.gates().values())


@pytest.mark.parametrize("order,key",((1,"primitive_first_time_integrand"),(3,"primitive_third_time_integrand")))
def test_actual_primitive_integrand_derivatives(order,key):
    expected=sp.diff(jets.model.coefficients()["boundary_primitive_s"],jets.u,order).subs(jets.u,0)
    assert sp.factor(expected-jets.data()[key])==0


@pytest.mark.parametrize("key",("actual_boundary_primitive_phi","actual_boundary_primitive_phi_phi_phi"))
def test_fixed_primitive_basepoint_is_not_moved_with_the_family(key):
    assert jets.data()[key].subs(jets.N,1)==0


def test_new_third_primitive_jet_has_the_independent_polynomial():
    N=jets.N
    expected=63*(16*N**sp.Rational(3,2)+9*N**4-30*N*N+5)/(5*N**sp.Rational(3,2))
    assert sp.factor(jets.data()["actual_boundary_primitive_phi_phi_phi"]-expected)==0


def test_original_clock_values_recovered():
    datum=existence.datum(0)
    assert datum["positive_conserved_matter_charge"]==sp.Rational(1,10)
    assert datum["central_physical_scale"]==1
    assert datum["central_lapse_acceleration"]==0
    assert datum["central_physical_Hubble_proper_acceleration"]==4


def test_fixed_phase_constraint_is_taken_before_family_substitution():
    assert existence.checks()["actual_central_fixed_phase_lapse_constraint"]==0
    bg=jets.background.data()
    assert bg["fixed_phase_lapse_Hessian"]!=sp.diff(
        jets.data()["positive_matter_density_squared"],jets.N)


def test_second_constraint_path_is_independently_replayed():
    rows=jets.checks()
    assert rows["independent_fixed_phase_constraint_path_Taylor"]==0
    assert rows["actual_lapse_acceleration_satisfies_twice_differentiated_constraint"]==0


def test_lapse_acceleration_has_the_full_polynomial_and_nonzero_offclock_value():
    N=jets.N
    numerator=3*N*(N*N-1)*(13921850250011*N**6-93228677254991*N**4+129756291744989*N*N-31453170750009)
    denominator=250000*(16874985*N**4-14244998*N*N-8625003)
    assert sp.factor(jets.data()["actual_lapse_second_time_derivative"]-numerator/denominator)==0
    assert existence.datum(existence.SHIFT_MAX)["central_lapse_acceleration"]<0


def test_omitting_the_third_primitive_jet_is_a_nonzero_negative_control():
    error=jets.data()["missing_third_primitive_jet_lapse_acceleration_error"]
    assert error.subs(jets.N,1+existence.SHIFT_MAX)!=0
    assert jets.checks()["missing_third_primitive_jet_negative_control"]==0


def test_physical_proper_acceleration_keeps_lapse_and_conformal_accelerations():
    d,N=jets.data(),jets.N
    expected=(d["actual_hat_Hubble_first_time_derivative"]+sp.Rational(3,2)*(1-N*N)
              +d["actual_lapse_second_time_derivative"]/(2*N))/N**2
    assert sp.factor(d["actual_physical_Hubble_proper_time_derivative"]-expected)==0
    assert sp.factor(d["actual_physical_Hubble_proper_time_derivative"]
                     -d["actual_hat_Hubble_first_time_derivative"])!=0


@pytest.mark.parametrize("key",("matter_density_squared","lapse_Hessian_times_N_seven_halves",
                              "lapse_acceleration","physical_Hubble_proper_acceleration"))
def test_exact_continuous_enclosures_ordered_and_nontrivial(key):
    box=jets.boxes()[key]
    assert box["lower"]<box["upper"]
    assert isinstance(box["lower"],sp.Rational)
    assert isinstance(box["upper"],sp.Rational)


def test_entire_central_interval_has_positive_charge_pivot_and_bounce_margin():
    box=jets.boxes()
    assert sp.Rational(1,200)<box["matter_density_squared"]["lower"]
    assert box["matter_density_squared"]["upper"]<sp.Rational(11,100)**2
    assert box["lapse_Hessian_times_N_seven_halves"]["upper"]<-2*(1+existence.SHIFT_MAX)**4
    assert box["physical_Hubble_proper_acceleration"]["lower"]>sp.Rational(39999,10000)


@pytest.mark.parametrize("shift",(sp.Rational(1,10**8),sp.Rational(1,10**6)))
def test_new_centers_have_constraint_compatible_nonzero_charge_and_positive_scale(shift):
    d=existence.datum(shift)
    assert d["positive_conserved_matter_charge_squared"]>0
    assert sp.simplify(d["positive_conserved_matter_charge"]**2-d["positive_conserved_matter_charge_squared"])==0
    assert d["central_physical_scale"]>1
    assert d["central_physical_Hubble_proper_acceleration"]>sp.Rational(39999,10000)
    assert d["initial_vector_field_and_momentum_zero"] is True


def test_complex_primitive_bound_is_fresh_not_the_old_real_time_derivative_bound():
    d=existence.coefficient_domain()
    assert d["actual_primitive_absolute_upper"]==sp.Rational(2,99)
    assert d["fresh_complex_time_primitive_phi_upper"]==sp.Rational(200,99)
    assert d["fresh_complex_time_primitive_phi_upper"]>2
    assert d["no_unjustified_extension_of_old_real_time_Hamiltonian_bound"] is True


def test_joint_complex_Hamiltonian_triangle_includes_full_primitive():
    d=existence.coefficient_domain()
    expected=2*(2*sp.Rational(5,4)**2+2*1002+2*sp.Rational(200,99)+sp.Rational(9,25)**2)
    assert d["homogeneous_Hamiltonian_triangle_upper"]==expected
    assert expected<d["declared_joint_holomorphic_Hamiltonian_upper"]


def test_recentered_outer_domain_fits_the_original_literal_branch():
    d=existence.coefficient_domain()
    assert d["outer_complex_lapse_radius_about_new_center"]+existence.SHIFT_MAX<sp.Rational(1,100)
    assert 2*d["outer_complex_time_radius"]<=d["original_primitive_complex_time_radius"]


def test_all_six_Cauchy_derivative_constants_include_factorials():
    d=existence.lapse_root()
    assert d["H_Nu_upper"]==800000000
    assert d["H_Np_and_H_Nell_upper"]==32000000
    assert d["H_NNN_upper"]==3840000000000
    assert d["H_NNu_upper"]==640000000000
    assert d["H_NNp_and_H_NNell_upper"]==25600000000
    assert d["inverse_exact_center_lapse_pivot_upper"]==sp.Rational(1,2)


def test_newton_fixed_point_and_pivot_are_uniform_on_the_whole_parameter_box():
    d=existence.lapse_root()
    assert d["Newton_contraction_upper"]<sp.Rational(1,10)
    assert d["Newton_closed_disc_image_radius_upper"]<existence.ROOT_RADIUS
    assert d["actual_root_displacement_upper"]<sp.Rational(1,10**15)
    assert d["remaining_lapse_pivot_absolute_lower"]>sp.Rational(19,10)


def test_literal_canonical_reduction_and_charge_reconstruction():
    rows=existence.checks()
    for key in ("literal_homogeneous_hat_scale_equation","literal_homogeneous_trace_density_equation",
                "literal_nonzero_matter_charge_conservation","literal_physical_free_matter_charge"):
        assert rows[key]==0


def test_actual_Hamiltonian_parity_and_odd_Hubble_are_checked():
    rows=existence.checks()
    assert rows["actual_time_and_phase_even_Hamiltonian"]==0
    assert rows["actual_time_and_phase_odd_hat_Hubble"]==0


def test_picard_domain_is_invariant_and_contractive_with_two_phase_columns():
    d=existence.evolution()
    assert d["phase_flow_Lipschitz_upper"]==2*sp.Integer(10)**29
    assert d["Picard_contraction_upper"]==sp.Rational(1,5)
    assert d["Picard_image_radius_upper"]==sp.Rational(5,10**26)
    assert d["Picard_image_radius_upper"]<d["solution_phase_disc_radius"]


def test_shorter_time_has_stronger_estimates_without_new_state_selection():
    first,second=existence.evolution(),existence.evolution(existence.TIME_WINDOW/2)
    assert second["Picard_contraction_upper"]*2==first["Picard_contraction_upper"]
    assert second["reconstructed_log_hat_scale_absolute_upper"]*2==first["reconstructed_log_hat_scale_absolute_upper"]
    assert second["quantitative_time_dependent_scalar_cone_or_time_advance_bound"] is False


def test_no_quantum_profile_or_response_is_transferred():
    assert jets.data()["quantum_profile_and_response_included"] is False
    assert existence.datum(0)["not_the_S6_82_profile_corrected_quantum_background"] is True


def test_native_domain_guards_reject_every_declared_bad_input_after_valid_calls():
    existence.datum(0)
    existence.evolution()
    assert audit.controls()["rejected_inputs"]==100


def test_exact_report_data_roundtrip_and_nonfinite_rejection():
    data=original.serialize({"jets":jets.data(),"boxes":jets.boxes(),"root":existence.lapse_root(),
                             "evolution":existence.evolution(),"domain":existence.coefficient_domain(),"enlarged":taylor.bounds(),
                             "central_coefficients":taylor.coefficients()})
    assert json.loads(json.dumps(data))==data
    for value in (1.0,sp.Float(1),sp.oo,sp.nan):
        with pytest.raises((TypeError,ValueError)):
            original.serialize(value)


def test_all_enlarged_literal_and_remainder_identities():
    rows=taylor.checks()
    assert len(rows)==18
    assert all(value==0 for value in rows.values())


@pytest.mark.parametrize("family,count",(("pivot",9),("center_force",7),("H",3),("Hp",2),("Hell",3)))
def test_each_actual_central_jet_family_has_its_nonzero_coefficients_and_exact_bounds(family,count):
    data=taylor.coefficients()[family]
    assert len(data["coefficients"])==count
    for row in data["coefficients"].values():
        assert row["enclosure"]["absolute_upper"]>=0
        assert sum(row["multi_degree"])<=data["degree"]
        assert row["coefficient_before_central_charge_substitution"]!=0


def test_radial_Cauchy_tail_has_no_missing_multivariate_multiplicity():
    b=taylor.bounds()
    assert b["joint_radial_Cauchy_ratio"]==sp.Rational(1,2500)
    assert b["center_force_radial_Cauchy_ratio"]==sp.Rational(1,5000)
    assert b["pivot_degree_two_tail_upper"]==sp.Rational(512,2499)
    assert b["center_force_degree_three_tail_upper"]==sp.Rational(1,156218750)
    assert b["Hamiltonian_degree_one_tail_upper"]==sp.Rational(4,2499)
    assert b["phase_derivative_degree_one_tail_upper"]==sp.Rational(32,2499)


def test_sharp_bounds_keep_the_full_analytic_tail_and_improve_the_declared_domain():
    b=taylor.bounds()
    assert b["lapse_root_radius"]==sp.Rational(1,10**6)
    assert b["complex_phase_parameter_radius"]==sp.Rational(1,10**5)
    assert b["lapse_pivot_change_upper"]>b["pivot_degree_two_tail_upper"]>0
    assert b["actual_Hamiltonian_absolute_upper"]<9
    assert b["actual_Hp_absolute_upper"]<sp.Rational(1,50)
    assert b["actual_Hell_absolute_upper"]<sp.Rational(1,5)


def test_enlarged_lapse_map_is_strictly_invariant_with_nonzero_remaining_pivot():
    b=taylor.bounds()
    assert b["Newton_contraction_upper"]<sp.Rational(1,8)
    assert b["Newton_closed_disc_image_radius_upper"]<taylor.ROOT_RADIUS
    assert b["actual_root_displacement_upper"]<sp.Rational(4,10**7)
    assert b["remaining_lapse_pivot_absolute_lower"]>sp.Rational(7,4)


def test_enlarged_actual_reduced_flow_not_the_original_coarse_flow_majorant():
    b=taylor.bounds()
    assert b["actual_trace_evolution_absolute_upper"]<10
    assert b["actual_density_evolution_absolute_upper"]<10
    assert b["declared_phase_flow_absolute_upper"]==10
    assert b["declared_phase_flow_absolute_upper"]<existence.FLOW_BOUND


def test_enlarged_picard_time_and_reconstruction_extend_the_old_same_solution():
    d=taylor.evolution()
    assert d["time_half_width"]==sp.Rational(1,10**7)
    assert d["Picard_image_radius_upper"]==sp.Rational(1,10**6)
    assert d["Picard_contraction_upper"]==sp.Rational(2,5)
    assert d["log_hat_scale_absolute_upper"]<sp.Rational(1,10**9)
    assert d["same_literal_classical_solution_extends_smaller_proved_interval"] is True
    assert d["not_a_quantitative_time_dependent_cone_or_UV_verdict"] is True


def test_enlarged_domain_guards_are_applied_before_cached_bounds():
    taylor.evolution()
    assert taylor.controls()["rejected_inputs"]==34
    smaller=taylor.evolution(taylor.TIME_WINDOW/2)
    assert smaller["Picard_contraction_upper"]==sp.Rational(1,5)
