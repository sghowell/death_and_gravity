"""Native complete global support, chart-pole controls and new scalar CCR checks."""
import json

import pytest
import sympy as sp
from p8_affine import verify as certificate
from p8_proca_global_support import (
    audit,
    charts,
    coverage,
    model,
    modes,
    phase,
    quantum,
    transfer,
)


def zero(value):
    return all(v==0 for v in value) if isinstance(value,sp.MatrixBase) else value==0


def test_named_and_scalar_identity_counts():
    rows=audit.residuals()
    assert len(rows)==55
    assert sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in rows.values())==349
    assert all(zero(v) for v in rows.values())


@pytest.mark.parametrize("name",tuple(audit.residuals()))
def test_each_exact_global_identity(name):
    assert zero(audit.residuals()[name])


def test_all_global_coverage_degree_and_pole_gates():
    rows=audit.gates()
    assert len(rows)==20
    assert all(v is True for v in rows.values())


def test_actual_global_clock_uses_new_total_margin_not_old_nearby_solution():
    d=model.background()
    assert d["new_total_margin"]==sp.Rational(1,200)
    assert d["a"]==(1+model.u**2)**2
    assert d["delta_J"]==1/(50*(1+model.u**2)**6)
    assert d["ell"]==1/(10*(1+model.u**2)**6)


def test_full_density_generator_is_polynomial_without_Theta_or_Lambda_poles():
    d=phase.data()
    assert len(d["complete_polynomial_q_coefficients"])==3
    o=model.old
    for v in d["regular_density_generator"]:
        den=sp.denom(sp.cancel(v))
        assert not den.has(o.theta,o.lam,o.q)
    assert d["normalized_physical_probe_force"]==sp.Matrix([0,0,0,o.a**3])


@pytest.mark.parametrize("name",("unitary","gamma"))
def test_complete_chart_maps_intertwine_and_have_verified_poles(name):
    d=transfer.data()[name]
    assert zero(d["complete_intertwining_residual"])
    assert d["forward_map_maximum_frequency_degree"]==2
    assert d["inverse_map_maximum_frequency_degree"]==1
    assert d["complete_friction_maximum_q_degree"]==0
    assert d["complete_force_remainder_maximum_q_degree"]==0
    assert d["unexpected_denominator_factors"]==[]
    assert d["all_endpoint_and_remainder_denominator_factors"]


def test_unitary_complete_remainders_are_exactly_independent_of_q():
    d=charts.data()["unitary"]
    assert d["complete_velocity_friction"].diff(model.q)==sp.zeros(2)
    assert d["complete_bounded_force_remainder"].diff(model.q)==sp.zeros(2)


def test_gamma_finite_frequency_pole_is_not_an_original_phase_pole():
    o=model.old
    center={symbol:value.subs(model.u,0) for symbol,value in model.substitution().items()}
    center[o.q]=sp.Rational(152,25)
    velocity=charts.data()["gamma"]["alpha_over_a_cubed"].subs(center)
    assert any(v.has(sp.zoo,sp.nan) for v in velocity)
    regular=phase.data()["regular_density_generator"].subs(center)
    assert all(v.is_finite is True for v in regular)


def test_global_complex_cutoff_separates_all_gamma_velocity_poles():
    d=coverage.data()
    assert d["gamma_complex_q_threshold_strict_upper"]<2000
    assert d["declared_global_complex_scalar_root_threshold"]==64
    assert 64**2/2>2000
    assert d["gamma_scale_squared_upper"]<2
    assert d["gamma_half_matter_mixing_squared_upper"]==sp.Rational(1,800)


@pytest.mark.parametrize("pair",((-10,10),(-10,-1),(1,10),(-1,0),(0,1),
    (-coverage.CUT,coverage.CUT),(0,0),(-10**6,10**6)))
def test_exact_global_time_partition_is_contiguous_and_at_most_three_segments(pair):
    d=coverage.time_pair(*pair)
    segments=d["ordered_chart_segments"]
    assert 1<=len(segments)<=3
    assert segments[0]["start"]==pair[0]
    assert segments[-1]["end"]==pair[1]
    for index,segment in enumerate(segments):
        assert segment["start"]<=segment["end"]
        if index:
            assert segments[index-1]["end"]==segment["start"]
        if segment["chart"]=="gamma":
            assert -sp.Rational(1,4)<=segment["start"]<=segment["end"]<=sp.Rational(1,4)
        else:
            assert segment["end"]<=-sp.Rational(1,8) or segment["start"]>=sp.Rational(1,8)
    assert d["support_radius_is_actual_integral_not_a_segmentwise_overestimate"] is True


def test_actual_global_null_radius_adds_exactly_across_both_overlaps():
    full=coverage.time_pair(-2,3)["actual_physical_matter_null_radius"]
    parts=coverage.time_pair(-2,-coverage.CUT)["actual_physical_matter_null_radius"]
    parts+=coverage.time_pair(-coverage.CUT,coverage.CUT)["actual_physical_matter_null_radius"]
    parts+=coverage.time_pair(coverage.CUT,3)["actual_physical_matter_null_radius"]
    assert sp.simplify(full-parts)==0
    assert coverage.time_pair(0,0)["actual_physical_matter_null_radius"]==0


def test_global_null_primitive_and_finite_conformal_length_are_actual_not_clock_duration():
    d=model.background()
    assert sp.factor(sp.diff(d["physical_matter_null_primitive"],model.u)-1/d["a"])==0
    assert coverage.time_pair(0,1)["actual_physical_matter_null_radius"]==sp.pi/8+sp.Rational(1,4)
    assert coverage.data()["whole_global_clock_conformal_time_length"]==sp.pi/2


BAD=(True,False,sp.true,sp.false,"1",1.0,sp.Float(1),None,[],{},
    sp.oo,-sp.oo,sp.zoo,sp.nan,sp.I,sp.sqrt(2),sp.Symbol("unfixed"))


@pytest.mark.parametrize("bad",BAD)
def test_invalid_exact_global_endpoints_rejected_after_valid_use(bad):
    assert len(coverage.time_pair(-10,10)["ordered_chart_segments"])==3
    for pair in ((bad,0),(0,bad)):
        with pytest.raises((TypeError,ValueError)):
            coverage.time_pair(*pair)


@pytest.mark.parametrize("pair",((1,0),(0,-1),(coverage.CUT,-coverage.CUT)))
def test_reversed_global_retarded_times_rejected(pair):
    with pytest.raises(ValueError):
        coverage.time_pair(*pair)


@pytest.mark.parametrize("name",("unitary","gamma"))
def test_principal_velocity_basis_has_no_inverse_clock_matter_gap(name):
    d=modes.data()
    c=d["charts"][name]
    assert sp.simplify(c["basis_determinant"]+4*sp.sqrt(d["clock_squared_speed"])/d["physical_scale"]**2)==0
    assert c["no_inverse_clock_matter_frequency_gap"] is True
    assert (c["exact_velocity_mode_basis"]*c["exact_velocity_mode_inverse"]).applyfunc(sp.simplify)==sp.eye(4)


def test_new_global_Gaussian_covariance_and_reference_time_are_distinct():
    d=quantum.data()
    C,B=d["new_global_ordered_initial_covariance"],d["new_global_initial_positive_Gram_factor"]
    assert (C-B*B.conjugate().T).applyfunc(sp.factor)==sp.zeros(4)
    assert d["new_global_reference_time"]==0
    assert "U_global(t,0,k)" in d["new_global_covariance_transport"]
    assert d["new_state_uses_global_evolution_not_either_frozen_nearby_state"] is True


def test_global_Kubo_and_density_CCR_sign_use_physical_volume():
    d=quantum.data()
    assert sp.factor(d["new_global_normalized_Kubo_prefactor"]*d["new_global_relational_commutator_multiplier"]
        -d["new_global_normalized_classical_response_multiplier"])==0
    assert d["new_global_normalized_Kubo_prefactor"]==quantum.state.kappa*d["new_global_unscaled_Kubo_prefactor"]
    C=d["new_global_ordered_initial_covariance"]
    assert (C-C.T-sp.I*quantum.state.hbar*phase.data()["constant_symplectic_form"]/quantum.state.kappa).applyfunc(sp.factor)==sp.zeros(4)


def test_global_schwartz_growth_keeps_all_three_segment_and_Duhamel_weights():
    d=quantum.data()
    assert d["compact_strip_transfer_polynomial_degree"]==9
    assert d["all_orders_Fourier_derivative_degree"]=="9+13*multi_index_order"


def test_all_invalid_global_endpoint_and_state_inputs_counted():
    d=audit.controls()
    assert d["global_endpoint_rejected_inputs"]==37
    assert d["new_global_Gaussian_inputs"]["rejected_inputs"]==17
    assert d["rejected_inputs"]==54


def test_all_complete_global_blocks_serialize_without_float_rounding():
    blocks=(model.background(),phase.data(),charts.data(),transfer.data(),
        coverage.data(),modes.data(),quantum.data())
    data=json.loads(json.dumps(certificate.serialize(blocks)))
    assert len(data)==7
    assert data[0]["new_total_margin"]=="1/200"
