"""Independent canonical normalization, Gaussian positivity and quantum-source tests."""
import json

import pytest
import sympy as sp
from p8_affine import verify as original
from p8_proca_compact_response import probes
from p8_proca_scalar_ccr import audit, canonical, growth, response, state


def zero(value):
    return all(entry==0 for entry in value) if isinstance(value,sp.MatrixBase) else value==0


def test_native_residual_and_scalar_counts():
    rows=audit.residuals()
    assert len(rows)==22
    assert sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in rows.values())==175
    assert all(zero(value) for value in rows.values())


@pytest.mark.parametrize("key",tuple(audit.residuals()))
def test_each_exact_identity_individually(key):
    assert zero(audit.residuals()[key])


def test_all_written_proof_and_scope_gates():
    assert len(audit.gates())==12
    assert all(value is True for value in audit.gates().values())


def test_literal_density_momenta_keep_all_three_volume_powers():
    d=canonical.data()
    assert d["normalized_to_density_phase"]==sp.diag(1,1,canonical.R**3,canonical.R**3)
    assert d["canonical_density_Hamiltonian"].has(canonical.piv,canonical.pic)
    assert not d["canonical_density_Hamiltonian"].has(canonical.parent.p,canonical.parent.parent.P)


def test_omitting_volume_derivative_breaks_canonical_generator():
    d=canonical.data()
    D=d["normalized_to_density_phase"]
    wrong=D*canonical.parent.data()["weighted_first_order_generator"]*D.inv()
    assert not zero((wrong-d["canonical_density_generator"]).applyfunc(sp.factor))
    Omega=d["constant_symplectic_form"]
    assert not zero((wrong*Omega+Omega*wrong.T).applyfunc(sp.factor))


def test_packet_symplectic_form_is_k_scaled_not_unit_canonical():
    d=canonical.data()
    E,Omega=d["canonical_density_to_packet_map"],d["constant_symplectic_form"]
    assert zero((E*Omega*E.T-canonical.k*Omega).applyfunc(sp.factor))
    assert not zero((E*Omega*E.T-Omega).applyfunc(sp.factor))


def test_density_source_keeps_physical_volume_factor():
    d=canonical.data()
    assert d["normalized_probe_source_vector"][3]==canonical.parent.data()["physical_volume_density"]
    assert d["normalized_probe_source_vector"][3].has(canonical.R**3)


def test_original_density_generator_has_no_inverse_momentum_chart():
    d=canonical.data()
    q=canonical.parent.parent.q
    assert all(sp.denom(sp.cancel(value)).has(q) is False for value in d["canonical_density_generator"])
    assert all(sp.degree(value,q)<=2 for value in d["canonical_density_generator"] if value!=0)


def test_explicit_map_growth_keeps_both_unequal_frequency_powers():
    d=growth.data()
    assert d["high_frequency_density_transfer_polynomial_degree"]==3
    assert d["declared_packet_from_density_norm_coefficient"]==24
    assert d["declared_density_from_packet_norm_coefficient"]==21
    assert d["high_frequency_density_transfer_triangle_coefficient"]==3024000000
    assert d["declared_high_frequency_density_transfer_coefficient"]==10**10


def test_full_low_frequency_bound_is_finite_without_exponential_rounding():
    bound=growth.data()["low_frequency_density_transfer_norm_upper"]
    assert bound.has(sp.exp)
    assert bound.is_finite is True
    assert not bound.has(sp.Float)


@pytest.mark.parametrize("square,expected",((0,1),(3,2),(sp.Rational(5,4),sp.Rational(3,2)),
    (sp.Symbol("q_squared",nonnegative=True),sp.sqrt(1+sp.Symbol("q_squared",nonnegative=True)))))
def test_exact_all_momentum_Gaussian_weight(square,expected):
    assert state.frequency_weight(square)==expected


@pytest.mark.parametrize("bad",(True,False,"1",1.0,sp.Float(1),None,[],{},-1,
    -sp.Rational(1,2),sp.oo,-sp.oo,sp.zoo,sp.nan,sp.I,
    sp.Symbol("unconstrained"),sp.Symbol("only_real",real=True)))
def test_invalid_Gaussian_weights_rejected_even_after_valid_call(bad):
    assert state.frequency_weight(0)==1
    with pytest.raises((TypeError,ValueError),match="exact nonnegative finite"):
        state.frequency_weight(bad)


def test_initial_covariance_positive_Gram_and_uncertainty_saturation():
    d=state.data()
    C,B=d["ordered_Cauchy_covariance"],d["positive_Cauchy_covariance_Gram_factor"]
    assert C.shape==(4,4)
    assert B.shape==(4,2)
    assert zero((C-B*B.conjugate().T).applyfunc(sp.factor))
    assert C.rank()==2
    assert d["symmetric_Cauchy_covariance"][0,0]*d["symmetric_Cauchy_covariance"][2,2]==state.hbar**2/(4*state.kappa**2)


def test_wrong_covariance_momentum_sign_breaks_required_CCR():
    C=state.data()["ordered_Cauchy_covariance"]
    Omega=canonical.data()["constant_symplectic_form"]
    assert not zero((C.T-C-sp.I*state.hbar*Omega/state.kappa).applyfunc(sp.factor))


def test_covariance_keeps_action_and_Planck_normalizations():
    d=state.data()
    assert d["commutator_normalization"]==state.hbar/state.kappa
    assert d["initial_reference_time"]==-sp.Rational(1,2*10**7)
    assert d["Hadamard_stress_tensor_and_interacting_quantum_matching_not_asserted"] is True


def test_Kubo_sign_and_volume_restore_positive_classical_response():
    d=response.data()
    result=sp.expand(d["normalized_probe_Kubo_prefactor"]*d["relational_matter_commutator_multiplier"])
    assert result==d["normalized_probe_classical_response_multiplier"]
    assert result!=-d["normalized_probe_classical_response_multiplier"]
    assert not result.has(state.hbar,state.kappa)


def test_unscaled_probe_retains_inverse_action_scale():
    d=response.data()
    result=sp.expand(d["unscaled_probe_Kubo_prefactor"]*d["relational_matter_commutator_multiplier"])
    assert result==d["normalized_probe_classical_response_multiplier"]/state.kappa
    assert result.has(state.kappa)


def test_compact_quantum_lower_is_the_same_parent_signal():
    d=response.data()
    lower=d["positive_compact_smeared_commutator_magnitude_lower"]
    assert lower==3*state.hbar*probes.SIGNAL_LOWER/(4*state.kappa)
    assert lower.is_positive is True
    assert d["positive_compact_classical_response_lower"]==3*probes.SIGNAL_LOWER/4


def test_control_count_and_no_renormalized_state_or_parent_transfer():
    assert audit.controls()["rejected_inputs"]==17
    assert audit.controls()["no_parent_source_or_scientific_library_changed"] is True
    assert state.data()["new_scalar_Gaussian_choice_not_a_transfer_of_old_Proca_state"] is True
    assert response.data()["no_interacting_parent_or_nonperturbative_quantum_gravity_locality_transfer"] is True


def test_exact_quantum_data_serializes_without_float_or_state_rounding():
    encoded=json.loads(json.dumps(original.serialize(state.data())))
    assert encoded["commutator_normalization"]==str(state.hbar/state.kappa)
    assert len(encoded["ordered_Cauchy_covariance"])==4
    assert "exp(" in json.dumps(original.serialize(growth.data()))
