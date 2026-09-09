"""Native complete new-action support, complex-frequency and fresh scalar CCR tests."""
import json

import pytest
import sympy as sp
from p8_affine import verify as certificate
from p8_proca_retuned_margin import central, model
from p8_proca_retuned_support import audit, bridge, domain, growth, momentum, quantum


def zero(value):
    return all(v==0 for v in value) if isinstance(value,sp.MatrixBase) else value==0


def test_named_and_scalar_exact_identity_counts():
    rows=audit.residuals()
    assert len(rows)==64
    assert sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in rows.values())==433
    assert all(zero(value) for value in rows.values())


@pytest.mark.parametrize("name",tuple(audit.residuals()))
def test_every_native_exact_identity(name):
    assert zero(audit.residuals()[name])


def test_all_fresh_coefficient_support_and_state_gates():
    rows=audit.gates()
    assert len(rows)==45
    assert all(value is True for value in rows.values())


def test_source_bridge_uses_new_constraint_and_literal_deformed_lapse_pivot():
    d=bridge.data()
    assert d["new_constraint_matter_square"]==model.system()["constraint_matter_square"].as_expr()
    assert d["new_constraint_matter_square"]!=model.old.system()["constraint_matter_square"].as_expr()
    assert len(d["full_original_polynomial_q_coefficient_matrices"])==3
    assert d["observable_is_new_relational_matter_linearization"]==bridge.principal.chi


def test_complete_packet_generator_keeps_all_four_lower_orders():
    d=bridge.finite.data()
    k=bridge.finite.k
    generator=bridge.data()["generic_complete_packet_Laurent_generator"]
    expected=k*d["J"]+d["L0"]+d["L1"]/k+d["L2"]/k**2+d["L3"]/k**3
    assert generator==expected
    assert all(d[name]!=sp.zeros(4) for name in ("L0","L1","L2","L3"))


def test_new_moving_gradient_keeps_actual_lapse_velocity_not_frozen_time():
    s=model.system()
    assert s["lapse_flow"]!=0
    assert s["gradient_lapse_velocity_part"]!=0
    assert zero(bridge.checks()["new_complete_moving_gradient_matches_new_rational_principal_formula"])


def test_new_basis_inverse_and_whole_disc_norms():
    b=domain.basis()
    assert (b["exact_new_action_basis"]*b["exact_new_action_inverse"]).applyfunc(sp.factor)==sp.eye(4)
    assert b["exact_basis_majorants"]["S"]["infinity_norm_upper"]==411*sp.sqrt(2)/80
    assert b["exact_basis_majorants"]["S_inverse"]["infinity_norm_upper"]==89*sp.sqrt(2)/32
    assert b["Cauchy_basis_time_derivative_infinity_norm_upper"]==160000000


def test_complete_low_complex_root_disc_has_explicit_original_generator_bound():
    d=domain.coefficient_majorants()
    assert len(d["complete_original_polynomial_q_coefficient_majorants"])==3
    assert d["low_scalar_root_unit_disc_generator_infinity_norm_upper"]==sp.Rational(28169,64)
    assert growth.data()["low_original_scalar_response_triangle_upper"]<4


def test_actual_new_action_coefficient_bounds_include_fixed_phase_and_matter_roots():
    d=domain.actual_bounds()
    assert d["new_alpha_modulus_upper"]<1
    assert d["actual_beta_modulus_upper"]<1
    assert d["new_complex_annuli"]["on_constraint_rescaled_lapse_Hessian"]["lower"]>2
    assert d["new_complex_annuli"]["constraint_matter_square"]["lower"]>0
    assert d["actual_source_output_prefactor_squared_upper"]<4


def test_full_complex_growth_has_exact_matter_radius_and_finite_lower_order_exponent():
    d=growth.data()
    assert d["declared_complete_lower_order_Euclidean_norm_upper"]==4000000000
    assert d["complete_scalar_source_output_prefactor"]==320
    assert d["full_interval_entire_multiplier_uniform_prefactor"]==4+320*sp.exp(400)
    radius=d["actual_matter_null_support_radius"]
    assert isinstance(radius,sp.Integral)
    assert radius.function==d["matter_frequency_uses_actual_new_lapse_conformal_and_hat_scale"]
    assert len(radius.function.atoms(sp.Function))==2


@pytest.mark.parametrize("pair",((-growth.HALF,growth.HALF),(-growth.HALF,0),(0,growth.HALF),(0,0)))
def test_valid_ordered_times_and_coincident_time_are_included(pair):
    d=growth.time_pair(*pair)
    assert d["coordinate_duration"]==pair[1]-pair[0]
    assert d["support_radius_is_actual_null_integral_not_this_upper_estimate"] is True
    assert d["actual_matter_null_radius_upper"]>=0


BAD=(True,False,"1",1.0,sp.Float(1),None,[],{},sp.oo,-sp.oo,sp.zoo,sp.nan,sp.I)


@pytest.mark.parametrize("bad",BAD)
def test_invalid_time_types_rejected_after_valid_evaluation(bad):
    assert growth.time_pair(0,0)["coordinate_duration"]==0
    for pair in ((bad,0),(0,bad)):
        with pytest.raises((TypeError,ValueError)):
            growth.time_pair(*pair)


@pytest.mark.parametrize("pair",((0,-growth.HALF),(-domain.T,0),(0,domain.T)))
def test_unordered_or_outside_certified_time_pair_rejected(pair):
    with pytest.raises(ValueError):
        growth.time_pair(*pair)


@pytest.mark.parametrize("bad",BAD+(0,-1,-sp.Rational(1,2)))
def test_invalid_complex_null_example_size_rejected(bad):
    with pytest.raises((TypeError,ValueError)):
        momentum.null_cone_example(bad)


def test_arbitrarily_large_complex_null_vector_requires_regular_original_chart():
    d=momentum.null_cone_example(10**30)
    assert d["bilinear_square"]==0
    assert d["Euclidean_modulus_squared"]==2*10**60
    assert d["covered_by_original_polynomial_generator_not_a_singular_chart"] is True


def test_complex_root_bound_does_not_identify_bilinear_and_Hermitian_squares():
    vector=sp.Matrix([2,sp.I,0])
    assert vector.dot(vector)==3
    assert (vector.conjugate().T*vector)[0]==5
    assert sp.sqrt(3)<sp.sqrt(5)
    assert len(momentum.checks())==7


def test_new_principal_gate_is_not_available_for_the_weaker_actual_fast_bounce():
    weak=central.diagnostic(sp.Rational(1,10000))
    assert weak["positive_matter_charge_possible"] is True
    assert weak["strict_classical_bounce_acceleration"] is True
    assert weak["clock_physical_squared_speed"]>1
    assert 0<domain.actual_bounds()["new_physical_clock_squared_speed_upper"]<1


def test_new_density_state_has_new_evolution_and_positive_Gram_factor():
    d=quantum.data()
    C,B=d["new_ordered_initial_covariance"],d["new_positive_initial_Gram_factor"]
    assert (C-B*B.conjugate().T).applyfunc(sp.factor)==sp.zeros(4)
    assert d["new_reference_time"]==-domain.T/2
    assert d["same_Gaussian_prescription_applied_to_new_evolution_not_old_state_transfer"] is True
    assert "U_new" in d["new_covariance_transport"]


def test_new_canonical_growth_is_polynomial_and_not_old_response_bound():
    d=quantum.data()
    assert d["new_real_packet_transfer_infinity_norm_upper"]==320*sp.exp(400)
    assert d["new_high_frequency_density_transfer_coefficient"]==161280*sp.exp(400)
    assert d["new_low_frequency_density_transfer_norm_upper"]<128
    assert d["new_all_frequency_density_transfer_degree"]==3
    assert d["new_all_order_Fourier_derivative_degree"]=="3+7*multi_index_order"


def test_new_Kubo_normalization_and_commutator_sign():
    d=quantum.data()
    hbar,kappa=quantum.state.hbar,quantum.state.kappa
    assert sp.factor(d["new_relational_matter_commutator_multiplier"]
        +sp.I*hbar*d["new_two_time_density_transfer"][1,3]/kappa)==0
    assert sp.factor(d["new_normalized_Kubo_prefactor"]*d["new_relational_matter_commutator_multiplier"]
        -d["new_normalized_probe_classical_response_multiplier"])==0
    assert d["new_normalized_Kubo_prefactor"]==kappa*d["new_unscaled_probe_Kubo_prefactor"]


def test_all_invalid_inputs_counted_and_equal_time_not_rejected():
    d=audit.controls()
    assert d["new_domain_rejected_inputs"]==32
    assert d["new_Gaussian_inputs"]["rejected_inputs"]==17
    assert d["rejected_inputs"]==49
    assert growth.time_pair(0,0)["coordinate_duration"]==0


def test_all_new_support_state_and_source_data_serialize_without_float_rounding():
    blocks=(bridge.data(),domain.actual_bounds(),domain.basis(),domain.coefficient_majorants(),
        growth.data(),momentum.data(),quantum.data())
    serialized=json.loads(json.dumps(certificate.serialize(blocks)))
    assert len(serialized)==7
    assert serialized[4]["full_interval_entire_multiplier_uniform_prefactor"]=="4 + 320*exp(400)"
