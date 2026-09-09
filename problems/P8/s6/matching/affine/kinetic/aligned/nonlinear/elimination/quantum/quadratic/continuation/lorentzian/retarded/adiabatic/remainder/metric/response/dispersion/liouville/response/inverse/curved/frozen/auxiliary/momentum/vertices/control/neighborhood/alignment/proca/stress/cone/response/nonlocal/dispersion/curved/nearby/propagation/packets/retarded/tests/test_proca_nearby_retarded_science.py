"""Independent native source, retarded transfer, shell and compact-probe tests."""
import json

import pytest
import sympy as sp
from p8_affine import verify as original
from p8_proca_nearby_retarded import (
    audit,
    fronts,
    remainder,
    source,
    support,
    transport,
)


def test_all_literal_algebra_and_scalar_entries():
    rows=audit.residuals()
    assert len(rows)==30
    assert sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in rows.values())==78
    assert all(v==(sp.zeros(*v.shape) if isinstance(v,sp.MatrixBase) else 0) for v in rows.values())


def test_all_source_response_and_distribution_gates():
    assert len(audit.gates())==30
    assert all(value is True for value in audit.gates().values())


@pytest.mark.parametrize("key",tuple(audit.residuals()))
def test_each_exact_identity_individually(key):
    value=audit.residuals()[key]
    assert value==(sp.zeros(*value.shape) if isinstance(value,sp.MatrixBase) else 0)


def test_original_hamiltonian_chart_really_has_nonzero_q_squared_piece():
    rows=source.data()["polynomial_q_coefficients"]
    assert len(rows)==3
    assert rows[2]!=sp.zeros(4)
    assert all(not value.has(source.parent.q) for matrix in rows for value in matrix)


def test_original_chart_is_regular_at_zero_q_without_a_Legendre_inverse():
    d=source.data()
    zero=d["weighted_first_order_generator"].subs(source.parent.q,0)
    assert not zero.has(sp.zoo,sp.oo,sp.nan)
    assert d["no_small_frequency_Legendre_inverse_or_division_by_Hubble"] is True
    assert d["q_zero_extension_used_only_for_bounded_frequency_Fourier_estimate"] is True


def test_both_shift_and_volume_derivatives_are_retained_in_canonical_map():
    Dd=source.data()["old_to_packet_map_dot"]
    assert Dd.has(source.parent.formula()["B_dot"])
    assert Dd.has(source.parent.H)
    assert source.checks()["full_time_dependent_canonical_map_reproduces_exact_packet_equations"]==sp.zeros(4)


def test_local_probe_source_contains_physical_volume_and_only_matter_canonical_force():
    d=source.data()
    assert d["physical_volume_density"]==source.R**3*source.N*source.e**3
    assert d["original_chart_source_vector"]==sp.Matrix([0,0,0,source.N*source.e**3])
    assert d["packet_chart_source_vector"]==source.R**sp.Rational(3,2)*d["original_chart_source_vector"]


def test_output_inverse_frequency_and_volume_are_not_dropped():
    row=source.data()["packet_chart_observable_row"]
    assert row==sp.Matrix([[0,1/(source.k*source.R**sp.Rational(3,2)),0,0]])


def test_both_endpoint_errors_remain_in_uniform_transfer_bound():
    d=transport.constants()
    assert d["both_endpoint_near_identity_matrix_error_coefficient"]==(
        6*transport.normal_form.Z_BOUND+4*transport.domains.T*transport.normal_form.REMAINDER_BOUND)
    assert d["near_identity_change_at_source_is_inverted_not_omitted"] is True
    assert d["scalar_error_is_absolute_and_uniform_in_phase_zeros"] is True


@pytest.mark.parametrize("k",(transport.normal_form.MIN_MOMENTUM,sp.Integer(10)**72))
def test_physical_transfer_and_scalar_errors_have_different_inverse_frequency_orders(k):
    low,high=transport.error(k),transport.error(2*k)
    assert high["declared_physical_transfer_error_upper"]*2==low["declared_physical_transfer_error_upper"]
    assert high["declared_scalar_multiplier_error_upper"]*4==low["declared_scalar_multiplier_error_upper"]
    assert low["actual_physical_transfer_error_upper"]<low["declared_physical_transfer_error_upper"]


def test_clock_front_amplitude_contains_both_endpoint_densities():
    d=fronts.data()
    assert d["clock_shell_amplitude"].has(fronts.source_parameters[0],fronts.target_parameters[0])
    assert not d["matter_shell_amplitude"].has(fronts.source_parameters[0],fronts.target_parameters[0])
    assert d["clock_shell_amplitude_lower"]==sp.Rational(1,10**8)
    assert d["matter_shell_amplitude_lower"]==sp.Rational(1,10**4)


def test_no_sine_zero_is_used_as_a_relative_error_denominator():
    d=transport.error(transport.normal_form.MIN_MOMENTUM)
    assert isinstance(d["declared_scalar_multiplier_error_upper"],sp.Rational)
    assert not d["declared_scalar_multiplier_error_upper"].free_symbols


def test_complete_low_frequency_matrix_has_all_three_nonnegative_majorants():
    d=remainder.low_frequency()
    assert len(d["complete_original_generator_coefficient_majorants"])==3
    assert all(value>=0 for mat in d["complete_original_generator_coefficient_majorants"] for value in mat)
    assert d["hat_q_upper_on_compact_frequency_ball"]==4*transport.normal_form.MIN_MOMENTUM**2


def test_huge_low_frequency_existence_bound_is_finite_and_not_rounded():
    d=remainder.low_frequency()
    B=d["complete_original_generator_infinity_norm_upper"]
    assert isinstance(B,sp.Rational) and B>0
    assert d["original_chart_fundamental_matrix_norm_upper"]==sp.exp(transport.domains.T*B)
    assert d["spatial_L2_squared_upper"].is_finite is True
    assert not d["spatial_L2_squared_upper"].has(sp.Float)


def test_high_frequency_L2_remainder_uses_three_dimensional_measure_and_Plancherel():
    d=remainder.data()
    C,K=transport.SCALAR_ERROR_COEFFICIENT,transport.normal_form.MIN_MOMENTUM
    assert d["high_frequency_multiplier_L2_squared_upper"]==4*sp.pi*C**2/K
    assert d["high_frequency_spatial_remainder_L2_squared_upper"]==C**2/(2*sp.pi**2*K)


def test_the_low_frequency_ball_is_not_discarded_from_the_response():
    d=remainder.data()
    low=d["low_frequency_complete_bound"]
    M=low["complete_exact_minus_two_front_multiplier_upper"]
    K=low["compact_frequency_ball_radius"]
    assert low["spatial_L2_squared_upper"]==K**3*M**2/(6*sp.pi**2)
    assert M>8*remainder.AMPLITUDE_UPPER*transport.domains.T


def test_spherical_shell_has_sine_over_k_not_an_extra_radius_factor():
    S,k=sp.symbols("S k",positive=True)
    assert sp.limit(sp.sin(k*S)/k,k,0)==S
    assert remainder.checks()["spherical_delta_shell_forward_Fourier_has_correct_radius_and_four_pi"]==0


def test_signed_normal_coordinate_does_not_simplify_delta_to_zero():
    assert remainder.checks()["surface_delta_pairing_does_not_acquire_normal_slab_width"]==0
    assert remainder.checks()["thin_normal_test_slab_L2_measure_scales_to_zero"]==0


def test_default_compact_neighborhoods_have_strict_matter_spacelike_margin():
    d=support.neighborhood()
    assert d["clock_minus_matter_ray_radius_lower"]==sp.Rational(99,8*10**14)
    assert d["common_time_and_Euclidean_space_neighborhood_radius"]==sp.Rational(99,8*10**16)
    assert d["every_pair_matter_spacelike_margin_lower"]==sp.Rational(4653,4*10**16)
    assert d["source_lower_time_margin"]>0
    assert d["detector_upper_time_margin"]>0


@pytest.mark.parametrize("start,end",((-support.domains.T/8,support.domains.T/8),(0,support.domains.T/4)))
def test_other_exact_interior_time_pairs_preserve_full_support_separation(start,end):
    d=support.neighborhood(start,end)
    assert d["every_pair_matter_spacelike_margin_lower"]>0
    assert d["every_pair_positive_clock_time_margin_lower"]>0


def test_bad_time_and_momentum_inputs_are_rejected_even_after_valid_calls():
    support.neighborhood()
    transport.error(transport.normal_form.MIN_MOMENTUM)
    assert audit.controls()["rejected_inputs"]==85
    with pytest.raises(ValueError,match="interior"):
        support.neighborhood(-support.domains.INNER_HALF_WIDTH,0)
    with pytest.raises(TypeError,match="exact"):
        support.neighborhood(-1.0,1.0)


def test_front_and_compact_probe_claims_do_not_assert_a_finite_EFT_or_quantum_transfer():
    d=support.data()
    assert d["unbounded_classical_frequency_front_is_not_finite_EFT_causality"] is True
    assert d["default_neighborhood"]["no_explicit_compact_probe_frequency_band_or_UV_matching_error"] is True
    assert audit.controls()["no_parent_action_state_profile_or_scientific_library_changed"] is True


def test_exact_large_rational_and_exponential_bounds_serialize_without_floats():
    d=json.loads(json.dumps(original.serialize(remainder.low_frequency())))
    assert d["bound_is_finite_not_a_numerical_stability_assertion"] is True
    assert "exp(" in d["original_chart_fundamental_matrix_norm_upper"]
