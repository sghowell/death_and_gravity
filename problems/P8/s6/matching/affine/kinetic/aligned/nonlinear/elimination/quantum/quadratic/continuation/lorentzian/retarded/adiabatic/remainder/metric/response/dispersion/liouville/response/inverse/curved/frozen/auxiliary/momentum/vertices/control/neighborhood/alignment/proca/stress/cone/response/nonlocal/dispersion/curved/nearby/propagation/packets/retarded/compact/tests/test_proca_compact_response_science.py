"""Independent exact compact-probe, signal, spatial-tail and conditional-budget tests."""
import json

import pytest
import sympy as sp
from p8_affine import verify as original
from p8_proca_compact_response import audit, band, bumps, matching, probes


def test_all_exact_scalar_identities():
    rows=audit.residuals()
    assert len(rows)==29
    assert all(value==0 for value in rows.values())


def test_all_declared_geometry_signal_and_conditional_gates():
    assert len(audit.gates())==29
    assert all(value is True for value in audit.gates().values())


@pytest.mark.parametrize("key",tuple(audit.residuals()))
def test_each_native_identity_individually(key):
    assert audit.residuals()[key]==0


@pytest.mark.parametrize("square,expected",((0,1),(1,0),(2,0),(sp.Rational(1,4),sp.exp(sp.Rational(-1,3)))))
def test_exact_bump_center_plateau_boundary_and_exterior(square,expected):
    assert bumps.beta_squared(square)==expected


def test_normalizers_are_actual_integrals_not_replaced_by_one():
    d=bumps.data()
    assert isinstance(d["one_dimensional_normalizer"],sp.Integral)
    assert d["three_dimensional_normalizer"].has(sp.Integral)
    assert d["actual_compact_source_formula"].has(bumps.A1,bumps.A3,bumps.rho)
    assert bumps.checks()["source_change_of_variables_keeps_all_four_normalization_powers"]==0


def test_bump_recurrence_keeps_the_boundary_denominator_and_exponential():
    d=bumps.data()
    assert len(d["first_boundary_recurrence_polynomials"])==5
    assert d["first_boundary_recurrence_polynomials"][1]==-2*bumps.r
    assert bumps.checks()["polynomial_times_boundary_exponential_has_zero_limit"]==0


def test_radial_and_detector_bumps_keep_global_piecewise_zero_extension():
    d=bumps.data()
    assert d["global_one_dimensional_bump"].has(sp.Piecewise)
    assert d["actual_compact_source_formula"].has(sp.Piecewise)
    assert d["actual_compact_detector_formula"].has(sp.Piecewise)
    assert d["spatial_radial_bumps_are_functions_of_squared_radius"] is True


def test_bump_domain_rejects_nonfinite_negative_unknown_and_floating_squares():
    bumps.beta_squared(0)
    assert bumps.controls()["rejected_inputs"]==16
    with pytest.raises(ValueError,match="finite nonnegative"):
        bumps.beta_squared(-1)
    with pytest.raises(TypeError,match="floating"):
        bumps.beta_squared(sp.Float(0))


def test_parent_supports_and_actual_clock_radii_are_not_changed():
    d=probes.formal()
    assert d["parent_support_neighborhood_radius"]==sp.Rational(99,8*10**16)
    assert d["detector_cap_and_time_width"]==d["parent_support_neighborhood_radius"]/10
    assert d["clock_ray_radius_lower"]==probes.T/16
    assert d["clock_ray_radius_upper"]==3*probes.T


def test_actual_width_retains_full_symbolic_low_frequency_exponential():
    d=probes.actual()
    assert d["complete_spatial_remainder_L2_upper"].has(sp.exp)
    assert d["exact_normal_width"]==probes.EPSILON_NUMERATOR/d["complete_spatial_remainder_L2_upper"]**2
    assert d["exact_source_time_and_spatial_width"]==d["exact_normal_width"]/20
    assert not d["exact_normal_width"].has(sp.Float)


def test_formal_probe_width_is_instantiated_by_the_actual_complete_remainder():
    formal=probes.formal()
    actual=probes.actual()
    substituted=formal["normal_width"].subs(probes.REMAINDER_NORM,actual["complete_spatial_remainder_L2_upper"])
    assert substituted==actual["exact_normal_width"]


def test_rational_width_upper_enforces_every_geometric_restriction():
    assert probes.EPSILON_UPPER<probes.A/100
    assert probes.A<probes.S_MIN/4
    assert probes.EPSILON_UPPER<probes.S_MIN/4
    assert 6*probes.A<probes.H


def test_positive_signal_lower_contains_all_four_bumps_and_clock_sphere_measure():
    d=probes.formal()
    assert d["clock_cap_detector_pointwise_lower"]==sp.Rational(1,16)
    assert d["clock_cap_solid_angle_lower"]==sp.pi*probes.A**2/(16*probes.S_MAX**2)
    assert probes.SIGNAL_LOWER==probes.AMP_LOWER*probes.S_MIN*probes.A**3/(1024*probes.S_MAX**2)
    assert probes.SIGNAL_LOWER>0


def test_detector_norm_scales_with_the_square_root_of_the_normal_width():
    d=probes.formal()
    assert sp.simplify(d["detector_L1_time_L2_space_norm_upper"]**2-
                       64*probes.A**4*d["normal_width"])==0


def test_source_norm_scales_with_minus_three_halves_of_its_own_width():
    d=probes.formal()
    assert sp.simplify(d["source_L1_time_L2_space_norm_upper"]**2-
                       64/d["source_time_and_spatial_width"]**3)==0


def test_complete_remainder_pairing_is_exactly_budgeted_by_quarter_signal():
    d=probes.formal()
    assert sp.simplify(probes.REMAINDER_NORM*d["detector_L1_time_L2_space_norm_upper"])==probes.SIGNAL_LOWER/4


def test_spatial_tail_contains_complete_source_detector_norm_product():
    d=band.data()
    assert d["compact_probe_L1_L2_norm_product_upper"]==6400*probes.A**2/d["normal_width"]
    assert 20**3<100**2


def test_spatial_cutoff_is_finite_but_not_rounded_or_called_an_EFT_cutoff():
    d=band.data()
    assert d["finite_comoving_spatial_cutoff"].is_finite is True
    assert (d["finite_comoving_spatial_cutoff"]-d["minimum_classical_high_frequency_domain"]).is_positive is True
    assert d["spatial_cutoff_not_a_temporal_frequency_or_interacting_EFT_cutoff"] is True
    assert d["projected_source_not_claimed_compact"] is True


def test_full_classical_and_low_spatial_band_pairings_are_both_positive():
    d=band.data()
    assert d["positive_full_classical_pairing_lower"]==3*probes.SIGNAL_LOWER/4
    assert d["positive_low_spatial_frequency_pairing_lower"]==probes.SIGNAL_LOWER/2


@pytest.mark.parametrize("low,tail",((sp.Rational(1,8),sp.Rational(1,8)),(0,sp.Rational(1,4)),
                                    (sp.Rational(1,10),sp.Rational(1,5))))
def test_distinct_exact_conditional_parent_error_budgets(low,tail):
    d=matching.budget(low,tail)
    assert d["conditional_parent_full_response_lower_fraction"]==sp.Rational(1,2)-low-tail
    assert d["conditional_parent_full_response_lower"]>0
    assert d["actual_matching_and_parent_tail_bounds_are_not_established"] is True


def test_zero_margin_conditional_budget_is_rejected_after_a_valid_call():
    matching.budget()
    with pytest.raises(ValueError,match="less than one half"):
        matching.budget(sp.Rational(1,4),sp.Rational(1,4))


def test_low_matching_and_parent_high_tail_are_independent_requirements():
    d=matching.data()
    assert len(d["required_common_parent"])==6
    assert d["low_frequency_agreement_alone_allows_high_frequency_cancellation"] is True
    assert d["algebraic_cancellation_control_is_not_a_constructed_UV_theory"] is True
    assert matching.checks()["low_band_agreement_only_has_an_exact_cancellation_control"]==0


def test_all_negative_controls_are_rejected():
    assert audit.controls()["rejected_inputs"]==51
    assert audit.controls()["no_parent_source_or_scientific_library_changed"] is True


def test_exact_actual_exponential_widths_serialize_without_float_conversion():
    d=json.loads(json.dumps(original.serialize(probes.actual())))
    assert "exp(" in d["exact_normal_width"]
    assert d["finite_symbolic_exponential_not_numerically_evaluated"] is True
