"""Independent exact finite-frequency, observable and normalized-packet tests."""
import json

import pytest
import sympy as sp
from p8_affine import verify as original
from p8_proca_nearby_packets import (
    audit,
    domains,
    finite,
    modes,
    normal_form,
    observable,
    packets,
)


def test_all_exact_identities_and_scalar_entries():
    rows=audit.residuals()
    assert len(rows)==45
    assert sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in rows.values())==273
    assert all(v==(sp.zeros(*v.shape) if isinstance(v,sp.MatrixBase) else 0) for v in rows.values())


def test_all_analytic_normalform_and_packet_gates():
    assert len(audit.gates())==51
    assert all(value is True for value in audit.gates().values())


@pytest.mark.parametrize("key",tuple(finite.checks()))
def test_every_exact_finite_frequency_canonical_identity(key):
    value=finite.checks()[key]
    assert value==sp.zeros(*value.shape)


def test_finite_q_construction_never_inverts_a_Legendre_Hessian():
    d=finite.data()
    assert d["A0"]*finite.parent.formula()["kinetic"]==sp.eye(2) or (
        d["A0"]*finite.parent.formula()["kinetic"]-sp.eye(2)).applyfunc(sp.factor)==sp.zeros(2)
    assert finite.parent.formula()["all_finite_q_Legendre_regular"] is False
    assert d["A2"]!=sp.zeros(2)
    assert d["E1"]!=sp.zeros(2)
    assert d["C0"]!=sp.zeros(2)


@pytest.mark.parametrize("index",range(4))
def test_each_action_normalized_leading_diagonal_transport_is_zero(index):
    assert modes.data()["B0"][index,index]==0


def test_complete_mode_basis_keeps_both_frequency_signs_and_action_signs():
    d=modes.data()
    assert list(d["Lambda"].diagonal())==[modes.wc,modes.wm,-modes.wc,-modes.wm]
    assert d["Krein_signature"]==sp.diag(-1,-1,1,1)
    assert modes.checks()["positive_and_negative_action_normalization"]==sp.zeros(4)


def test_actual_complex_time_disc_is_not_a_loose_independent_phase_derivative_bound():
    d=domains.analytic_domain()
    assert d["actual_complex_solution_time_radius"]==sp.Rational(1,10**7)
    assert d["inner_real_time_half_width"]==sp.Rational(1,2*10**7)
    assert d["actual_complex_solution_phase_image_upper"]==sp.Rational(1,10**6)
    assert d["actual_solution_rescaled_trace_upper"]<sp.Rational(3,10**5)


@pytest.mark.parametrize("key",("clock_physical_speed_squared_excess",
    "on_constraint_rescaled_lapse_Hessian","constraint_matter_square","M"))
def test_complex_modulus_annuli_are_separated_not_real_order_assumptions(key):
    row=domains.analytic_domain()["full_complex_rational_annuli"][key]
    assert isinstance(row["lower"],sp.Rational)
    assert 0<row["lower"]<row["upper"]


def test_complex_annulus_handles_negative_centers_exactly():
    row=domains.complex_annulus({"numerator":{"center":-5,"deviation_upper":1},
                                "denominator":{"center":-3,"deviation_upper":1}})
    assert row=={"lower":sp.Integer(1),"upper":sp.Integer(3)}
    assert all(isinstance(value,sp.Rational) for value in row.values())


def test_complex_denominator_and_float_controls_are_not_skipped_after_valid_call():
    valid={"numerator":{"center":2,"deviation_upper":0},"denominator":{"center":3,"deviation_upper":0}}
    assert domains.complex_annulus(valid)["lower"]==sp.Rational(2,3)
    with pytest.raises(ValueError,match="zero"):
        domains.complex_annulus({"numerator":{"center":1,"deviation_upper":0},"denominator":{"center":1,"deviation_upper":1}})
    with pytest.raises(TypeError,match="exact"):
        domains.complex_annulus({"numerator":{"center":1.0,"deviation_upper":0},"denominator":{"center":2,"deviation_upper":0}})


def test_conformal_and_scale_bounds_retain_full_analytic_tails():
    d=domains.analytic_domain()
    assert d["e_minus_one_upper"]>d["e_fourth_minus_one_upper"]>0
    assert d["e_minus_one_upper"]<sp.Rational(4,10**6)
    assert d["R_minus_one_upper"]>d["actual_complex_log_scale_upper"]>0
    assert d["R_minus_one_upper"]<sp.Rational(2,10**9)


@pytest.mark.parametrize("key",("A0","A1","A2","E0","E1","C0"))
def test_every_finite_frequency_Laurent_piece_has_explicit_entry_majorants(key):
    matrix=domains.coefficient_majorants()[key]
    assert matrix.shape==(2,2)
    assert all(0<=value<100 for value in matrix)


def test_actual_basis_and_inverse_majorants_are_the_explicit_root_triangle_sums():
    d=domains.basis_majorants()
    assert d["S"]["infinity_norm_upper"]==85*sp.sqrt(2)/2
    assert d["S_inverse"]["infinity_norm_upper"]==89*sp.sqrt(2)/4
    assert max(row["infinity_norm_upper"] for row in d.values())<domains.BASIS_NORM


def test_signed_frequency_gaps_are_uniformly_separated_in_the_complex_disc():
    d=domains.analytic_domain()
    assert d["same_sign_frequency_gap_lower"]>domains.GAP
    assert d["opposite_sign_frequency_gap_lower"]>domains.GAP
    assert d["matter_coordinate_frequency_modulus_lower"]>sp.Rational(1,2)
    assert d["matter_coordinate_frequency_modulus_upper"]<2


def test_time_Cauchy_derivatives_keep_second_factorial():
    d=normal_form.constants()
    assert d["basis_and_inverse_first_derivative_upper"]==2*10**10
    assert d["basis_second_derivative_upper"]==8*10**17
    assert d["L0_first_derivative_upper"]==2*10**11
    assert d["leading_transport_triangle_upper"]<3*10**13
    assert d["leading_transport_first_derivative_triangle_upper"]<2*10**21


def test_near_identity_normalform_keeps_gap_and_change_time_derivatives():
    rows=normal_form.checks()
    assert rows["off_diagonal_transport_removed_by_actual_gap_matrix"]==sp.zeros(4)
    assert rows["exact_time_dependent_near_identity_remainder"]==sp.zeros(4)
    assert rows["gap_derivative_is_retained_in_change_derivative"]==0
    d=normal_form.constants()
    assert d["off_diagonal_change_derivative_triangle_upper"]>0
    assert d["exact_transformed_remainder_triangle_upper"]>normal_form.Z_DOT_BOUND


def test_normalform_error_has_no_exponential_of_carrier_or_crude_transport_norm():
    low=normal_form.error(normal_form.MIN_MOMENTUM)
    high=normal_form.error(2*normal_form.MIN_MOMENTUM)
    assert low["integrated_remainder_upper"]==sp.Rational(1,4)
    assert low["phase_propagator_norm"]==1
    assert low["no_exponential_of_carrier_or_B0_norm"] is True
    assert high["modal_phase_approximation_error_upper"]*2==low["modal_phase_approximation_error_upper"]
    assert high["full_U_error_upper"]<low["full_U_error_upper"]


def test_relational_observable_and_formal_probe_have_the_exact_gauge_and_Ward_identities():
    rows=observable.checks()
    assert rows["linear_relational_matter_observable_is_gauge_invariant"]==0
    assert rows["unitary_gauge_is_the_reconstructed_matter_coordinate"]==0
    assert rows["formal_clock_and_matter_probe_Noether_forces_cancel"]==0
    assert rows["vanishing_background_observable_has_no_linear_metric_probe_force"]==0


def test_actual_clock_pole_in_matter_observable_is_nonzero_and_not_cancelled():
    assert observable.checks()["full_principal_matter_response_has_both_distinct_positive_poles"]==0
    assert observable.data()["clock_principal_pole_weight_lower"]==sp.Rational(691200000,1000002000001)
    assert observable.data()["clock_principal_pole_weight_lower"]>sp.Rational(1,2000)


@pytest.mark.parametrize("exponent",(72,75,120))
def test_concrete_real_three_dimensional_packet_has_controlled_band_errors_and_tails(exponent):
    d=packets.packet(exponent)
    assert d["envelope_bandwidth"]**3==d["classical_carrier"]
    assert d["comoving_frequency_band"][0]>normal_form.MIN_MOMENTUM
    assert d["total_relative_field_error_upper"]<sp.Rational(1,10**27)
    assert d["leading_real_packet_outside_cube_L2_mass_upper"]<sp.Rational(1,10**9)
    assert d["exact_packet_normalized_outside_cube_L2_mass_upper"]<sp.Rational(1,10**8)


def test_real_fourier_carriers_are_disjoint_and_choose_conjugate_outgoing_modes():
    d=packets.packet()
    assert d["classical_carrier"]>3*d["envelope_bandwidth"]
    assert (d["positive_carrier_outgoing_clock_mode_index"],d["negative_carrier_conjugate_clock_mode_index"])==(2,0)
    assert (d["positive_carrier_outgoing_matter_mode_index"],d["negative_carrier_conjugate_matter_mode_index"])==(3,1)


def test_full_sinc_normalization_tail_and_three_dimensional_diffraction_are_checked():
    rows=packets.checks()
    for key in ("flat_spectral_envelope_inverse_Fourier_integral","one_dimensional_spectral_envelope_unit_norm",
                "three_dimensional_product_envelope_unit_norm","full_two_sided_sinc_tail_integral",
                "spatial_dispersion_difference_exact_rationalization"):
        assert rows[key]==0
    d=packets.packet()
    assert d["diffraction_relative_error_upper"]==sp.Rational(4,10**31)


def test_final_clock_cube_is_outside_matter_future_of_initial_cube_with_positive_margin():
    d=packets.packet()
    assert d["classical_ray_separation_lower"]==sp.Rational(99,4*10**14)
    assert d["comparison_cube_half_width"]==sp.Rational(99,32*10**14)
    assert d["clock_cube_margin_beyond_matter_future_of_initial_cube"]>0


def test_normalized_exact_packet_mass_not_unnormalized_approximation_mass():
    d=packets.packet()
    assert d["exact_packet_normalized_outside_cube_L2_mass_upper"]==8*(
        d["leading_real_packet_outside_cube_L2_mass_upper"]+d["total_relative_field_error_upper"]**2)
    assert d["exact_packet_normalized_outside_cube_L2_mass_upper"]>d["leading_real_packet_outside_cube_L2_mass_upper"]>0


def test_finite_band_tails_are_not_compact_preparation_or_EFT_authority():
    d=packets.packet()
    assert d["compact_spatial_support_or_compact_local_preparation"] is False
    assert d["interacting_EFT_band_or_quantum_resolution"] is False
    assert d["UV_exclusion_or_original_P8_closure"] is False
    assert observable.data()["no_compact_source_retarded_preparation_or_quantum_detector_theorem"] is True


@pytest.mark.parametrize("carrier",(normal_form.MIN_MOMENTUM,normal_form.MIN_MOMENTUM+1))
def test_observable_error_requires_the_whole_radial_band_not_only_carrier(carrier):
    with pytest.raises(ValueError,match="complete packet band"):
        observable.relative_error(carrier)


def test_all_native_domain_guards_after_valid_calls():
    domains.analytic_domain()
    normal_form.error(sp.Integer(10)**72)
    packets.packet()
    assert audit.controls()["rejected_inputs"]==137


def test_exact_report_data_roundtrip_and_nonfinite_rejection():
    data=original.serialize({"finite":finite.data(),"basis":modes.data(),"domain":domains.analytic_domain(),
        "constants":normal_form.constants(),"observable":observable.data(),"packet":packets.packet()})
    assert json.loads(json.dumps(data))==data
    for value in (1.0,sp.Float(1),sp.oo,sp.nan):
        with pytest.raises((TypeError,ValueError)):
            original.serialize(value)
