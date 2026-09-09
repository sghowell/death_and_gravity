"""Independent new minimal-Proca prepared response and fixed-profile tests."""
import json

import pytest
import sympy as sp
from p8_affine import verify as original
from p8_proca_local_response import bounds as local_bounds
from p8_proca_nonlocal_response import (
    audit,
    envelopes,
    evolution,
    source,
    tadpole,
    tail,
    tangent,
)


def test_all_new_prepared_response_exact_residuals():
    rows=audit.residuals()
    assert len(rows)==112
    assert sum(len(value) if isinstance(value,sp.MatrixBase) else 1 for value in rows.values())==133
    for value in rows.values():
        assert all(x==0 for x in value) if isinstance(value,sp.MatrixBase) else value==0


def test_all_continuous_response_and_scope_gates():
    assert len(audit.gates())==41
    assert all(value is True for value in audit.gates().values())


@pytest.mark.parametrize("sector",("T","L"))
def test_literal_constant_mass_frequency_and_moving_canonical_sources(sector):
    data=source.data(sector)
    assert data["r"]==source.n[0]-source.z*source.v[0]
    assert data["rg"].coeff(source.n[0])==-sp.Rational(1,2)
    assert all(value==0 for value in source.ZERO.values())
    rows=source.checks()
    assert rows[sector+"_literal_minimal_frequency_variation"]==0
    assert rows[sector+"_literal_minimal_canonical_map_variation"]==0


@pytest.mark.parametrize("sector",("T","L"))
def test_actual_readout_contacts_and_physical_normalization(sector):
    data=source.readout(sector,"energy")
    assert data["A"]==data["B"]==1
    assert data["delta_A"]==data["delta_B"]==-source.n[0]-3*source.v[0]
    for component in ("energy","pressure"):
        for key in ("A","B"):
            assert source.checks()[sector+"_"+component+"_literal_contact_"+key]==0


@pytest.mark.parametrize("sector",("T","L"))
def test_covariance_Ward_identity_is_checked_before_and_after_variation(sector):
    rows=source.checks()
    assert rows[sector+"_ordinary_mode_covariance_Ward"]==sp.zeros(1,3)
    assert rows[sector+"_actual_varied_mode_covariance_Ward"]==sp.zeros(1,3)


@pytest.mark.parametrize("sector",("T","L"))
def test_actual_new_energy_adiabatic_baseline_and_independent_variation(sector):
    for order in range(3):
        assert source.checks()[sector+"_actual_new_energy_adiabatic_baseline_"+str(order)]==0
        assert source.checks()[sector+"_actual_varied_adiabatic_Ward_"+str(order)]==0
    for order in (2,4):
        assert tangent.low_checks()[sector+"_frozen_dimensional_variation_order_"+str(order)]==0


@pytest.mark.parametrize("sector",("T","L"))
def test_reference_residual_has_no_lower_inverse_frequency_power(sector):
    data=tangent.reference(sector)
    assert data["varied_residual_low_coefficients"]=={j:0 for j in range(4)}
    assert all(value==0 for value in envelopes.reference(sector)["box_reconstructions"])


@pytest.mark.parametrize("sector",("T","L"))
def test_full_box_tenth_source_coefficients_differ_from_old_mass_response(sector):
    data=tangent.reference(sector)
    expected=-(1-source.z)/512 if sector=="T" else -(1+source.z)/512
    assert data["tenth_source_derivative_coefficients"]["N"]==0
    assert sp.factor(data["tenth_source_derivative_coefficients"]["Z"]-expected)==0
    assert data["tenth_source_derivative_fixtures"]=={"N":0,"Z":-sp.Rational(1,512)}


@pytest.mark.parametrize("sector,component,expected",(
    ("T","energy",774182467),("L","energy",754679053),
    ("T","pressure",476096177),("L","pressure",566221945)))
def test_all_four_integrable_subtracted_reference_tails(sector,component,expected):
    data=tail.reference_tail(sector,component)
    assert data["low_tail_numerator_coefficients"]=={j:0 for j in range(3)}
    assert data["varied_reference_tail_over_inverse_frequency_fifth_upper"]==expected
    assert data["all_majorants_nonnegative"] is True


def test_recomputed_transport_constants_and_all_three_polarizations():
    data=evolution.constants()
    common=data["common"]
    assert common["delta_mixing_ninth_inverse_frequency_coefficient"]==8296904784986635
    assert common["delta_mixing_eighth_inverse_frequency_coefficient"]==260404337448359
    for sector in ("T","L"):
        for component in ("energy","pressure"):
            assert data["by_sector"][sector]["physical_readouts"][component]["delta_reference_bilinear_over_nu_squared_upper"]==30


def test_nonzero_initial_mixing_and_phase_interference_are_retained():
    data=evolution.constants()["frozen_initial_and_evolved_mixing"]
    assert data["initial_mixing_envelopes"][6]==1813229
    assert data["evolution_mixing_envelope"]==5347035781757616
    assert evolution.algebra_checks()["nonzero_initial_mixing_phase_control"]==0


def test_exact_radial_integrals_not_a_momentum_grid():
    rows=evolution.algebra_checks()
    for power in (4,7,8):
        assert rows["radial_mixing_integral_"+str(power)]==0
    assert rows["radial_reference_tail_integral_five"]==0


def test_complete_vector_bounds_use_the_new_local_physical_response():
    data=evolution.bound(local_bounds.SCALE)
    local=local_bounds.at_scale(local_bounds.SCALE)["response_bounds"]["physical_finite_local_response"]
    for component in ("energy","pressure"):
        assert data["finite_local_component_bounds"][component]==local[component]["C4_to_C0_upper"]
        assert data["complete_metric_response_C10_to_C0_component_bounds"][component]<sp.Rational(1,10**790)


def test_actual_fixed_profile_contacts_and_normalized_outputs():
    data=tadpole.physical_vertices()
    rho,p=data["rho"],data["pressure"]
    assert data["fixed_profile_density_gradient"]==sp.Matrix([rho,-3*p])
    assert data["fixed_profile_density_hessian"]==sp.Matrix([[-rho-p,3*rho],[3*rho,-9*p]])
    assert data["fixed_profile_physical_stress_jacobian"]==sp.Matrix([[rho+p,0],[rho+p,0]])


def test_total_chart_pullback_cancels_second_map_contact_only_after_profile():
    data=tadpole.physical_vertices()
    assert data["checks"]["total_second_metric_map_contact_cancels_only_after_fixed_profile"]==0
    assert data["checks"]["actual_background_cancelled_clock_current_map"]==sp.zeros(2)
    assert data["fixed_profile_density_gradient"]!=sp.zeros(2,1)


def test_continuous_actual_clock_chart_C10_lift():
    data=tadpole.chart_lift()
    assert data["joint_input_C10_upper"]==sp.Rational(46090764897,8)
    assert data["joint_output_C0_upper"]==3
    assert data["output_bound_uses_omega_N_at_most_one_half"] is True
    assert len(data["Leibniz_C10_row_bounds"])==11
    assert all(row["reconstruction"]==0 for row in data["coefficient_derivative_envelopes"].values())


def test_full_fixed_profile_plus_vector_bounds_in_both_charts():
    data=tadpole.bound(local_bounds.SCALE)
    assert data["M_tau"]==10**400
    assert data["fixed_mass_time_product"]==1000
    assert data["joint_background_cancelled_physical_C10_to_C0_upper"]<sp.Rational(1,10**790)
    assert data["joint_background_cancelled_clock_C10_to_C0_upper"]<sp.Rational(1,10**779)
    assert data["no_no_loss_or_coupled_inverse_inference"] is True


def test_scale_change_does_not_rescale_the_fixed_mass():
    a,b=tadpole.bound(local_bounds.SCALE),tadpole.bound(2*local_bounds.SCALE)
    assert a["fixed_mass_time_product"]==b["fixed_mass_time_product"]==1000
    assert a["joint_background_cancelled_clock_C10_to_C0_upper"]==4*b["joint_background_cancelled_clock_C10_to_C0_upper"]


def test_tenth_derivative_is_controlled_but_eleventh_is_rejected():
    assert source.linear_bound(source.v[10])["highest_source_derivative"]==10
    with pytest.raises(ValueError,match="C10"):
        source.linear_bound(source.v[11])


def test_native_guards_after_scientific_cache_population():
    source.data("T")
    source.adiabatic("T","energy",1)
    tangent.reference("T")
    envelopes.reference("T")
    evolution.bound(local_bounds.SCALE)
    assert audit.controls()["rejected_inputs"]==146


def test_exact_response_json_roundtrip_and_float_rejection():
    data=original.serialize({"full":tadpole.bound(local_bounds.SCALE),"transport":evolution.constants()})
    assert json.loads(json.dumps(data))==data
    for value in (1.0,sp.Float(1),sp.oo,sp.nan):
        with pytest.raises((TypeError,ValueError)):
            original.serialize(value)

