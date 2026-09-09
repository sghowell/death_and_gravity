"""Independent exact complex-jet, complete-frequency and new compact-probe tests."""
import json

import pytest
import sympy as sp
from p8_affine import verify as original
from p8_proca_compact_response import bumps
from p8_proca_scalar_ccr import state
from p8_proca_sharp_response import audit, basis, jets, probes, transfer


def zero(value):
    return all(v==0 for v in value) if isinstance(value,sp.MatrixBase) else value==0


def test_exact_named_and_scalar_identity_counts():
    rows=audit.residuals()
    assert len(rows)==25
    assert sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in rows.values())==70
    assert all(zero(value) for value in rows.values())


@pytest.mark.parametrize("name",tuple(audit.residuals()))
def test_each_exact_identity(name):
    assert zero(audit.residuals()[name])


def test_all_actual_complex_transfer_and_probe_gates():
    assert len(audit.gates())==76
    assert all(value is True for value in audit.gates().values())


@pytest.mark.parametrize("name",tuple(jets.PARTIAL_CAPS))
def test_every_logarithmic_partial_cap_covers_native_whole_box(name):
    actual=jets.native_bounds()["logarithmic_partial_upper"][name]
    assert all(a<=b for a,b in zip(actual,jets.PARTIAL_CAPS[name],strict=True))


def test_logarithmic_modulus_bound_is_invariant_under_fraction_sign():
    value=(2+jets.parent.uf)/(3+jets.parent.Nf)
    assert jets.logarithmic_partial_upper(value)==jets.logarithmic_partial_upper(-value)
    assert jets.modulus_upper(value)==jets.modulus_upper(-value)


@pytest.mark.parametrize("case",("zero","time_zero","lapse_zero","pole"))
def test_whole_box_logarithmic_zeros_and_rational_poles_rejected(case):
    p=jets.parent
    values={"zero":p.FIELD.zero,"time_zero":p.uf,
        "lapse_zero":p.Nf-p.FIELD.from_expr(jets.box.LAPSE_CENTER),"pole":1/p.uf}
    with pytest.raises(ValueError):
        (jets.modulus_upper if case=="pole" else jets.logarithmic_partial_upper)(values[case])


def test_logarithmic_velocity_retains_lapse_trace_and_conformal_motion():
    d=jets.data()
    assert d["declared_lapse_velocity_upper"]==sp.Rational(7,50000)
    assert d["declared_trace_velocity_upper"]==9
    assert d["total_logarithmic_moduli"]["clock_physical_speed_squared"]==sp.Rational(363,100000)
    assert d["clock_action_amplitude_log_velocity_triangle_upper"]<sp.Rational(1,400)
    assert d["signed_mode_frequency_velocity_triangle_upper"]<sp.Rational(1,100)


def test_small_basis_bounds_are_exact_Laurent_row_sums():
    d=basis.data()
    expected={"S":sp.Rational(331,80)*sp.sqrt(2),
        "S_inverse":sp.Rational(73,32)*sp.sqrt(2),
        "S_dot":sp.Rational(200407,20000000)*sp.sqrt(2)}
    for name,value in expected.items():
        assert d["basis_majorants"][name]["infinity_norm_upper"]==value
    assert all(d["complete_basis_time_derivative"].has(v) for v in (basis.lc,basis.lm,basis.elld))


def test_actual_scalar_prefactor_retains_both_volumes_lapse_and_cube_conformal():
    d=basis.data()
    assert d["actual_two_endpoint_scalar_prefactor_squared_upper"]==sp.Rational(101,99)**3*sp.Rational(101,100)**8
    assert d["actual_two_endpoint_scalar_prefactor_squared_upper"]<4


@pytest.mark.parametrize("bad",(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,
    sp.zoo,None,sp.Symbol("unbounded"),sp.sqrt(2),[],0,-1,transfer.MIN_MOMENTUM-1))
def test_new_high_frequency_domain_rejects_bad_values_after_valid_use(bad):
    assert transfer.momentum(transfer.MIN_MOMENTUM)==10**14
    with pytest.raises((TypeError,ValueError)):
        transfer.error(bad)


@pytest.mark.parametrize("k",(10**14,10**15,sp.Rational(10**30+1,2)))
def test_new_high_frequency_constants_hold_on_their_own_domain(k):
    d=transfer.error(k)
    assert d["near_identity_difference_upper"]<sp.Rational(1,2)
    assert d["integrated_normal_form_remainder_upper"]<=sp.Rational(1,4)
    assert d["uniform_two_endpoint_scalar_response_error_upper"]==sp.Integer(10)**16/sp.Rational(k)**2


def test_new_time_Cauchy_bound_uses_complex_B0_not_real_point_samples():
    d=transfer.data()
    assert d["Cauchy_B0_first_derivative_upper"]==transfer.B0_UPPER/(jets.domains.T/2)
    assert d["actual_Z_dot_triangle_upper"]>transfer.B0_DOT_UPPER/jets.domains.GAP


def test_uniform_transfer_does_not_exponentiate_large_k_polynomial():
    d=transfer.data()
    assert d["all_k_at_least_one_integrated_generator_upper"]==sp.Rational(1,10)
    assert d["exact_interaction_picture_transfer_norm_upper"]==sp.Rational(10,9)
    assert d["declared_all_k_at_least_one_packet_transfer_norm_upper"]==48
    assert d["unit_low_frequency_ball_original_generator_norm_upper"]==sp.Rational(28169,64)


def test_all_four_literal_L_matrices_are_bounded_with_their_actual_contacts():
    rows=transfer.coefficient_majorants()
    assert set(rows)=={"L0","L1","L2","L3"}
    assert all(max(sum(m[i,j] for j in range(4)) for i in range(4))<10**4 for m in rows.values())

def test_remainder_has_three_nonoverlapping_complete_frequency_pieces():
    d=probes.data()
    K,C=transfer.MIN_MOMENTUM,transfer.SCALAR_ERROR
    assert d["complete_remainder_low_unit_ball_L2_squared_upper"]==25/(6*sp.pi**2)
    assert d["complete_remainder_intermediate_band_L2_squared_upper"]==110**2*(K-1)/(2*sp.pi**2)
    assert d["complete_remainder_high_tail_L2_squared_upper"]==C*C/(2*sp.pi**2*K)
    assert d["complete_remainder_L2_squared_rational_upper"]<10**18


def test_new_widths_are_exact_rational_and_not_the_old_literal_pair():
    assert isinstance(probes.EPSILON,sp.Rational)
    assert isinstance(probes.RHO,sp.Rational)
    assert probes.EPSILON==probes.parent.EPSILON_NUMERATOR/10**18
    assert probes.RHO==probes.EPSILON/20
    assert probes.EPSILON>probes.parent.EPSILON_UPPER
    assert probes.data()["same_actual_background_and_retarded_kernel_but_different_probe_widths"] is True


def test_new_actual_piecewise_formulas_keep_normalizers_and_set_all_new_widths():
    d=probes.data()
    source=d["actual_new_compact_source_formula"]
    detector=d["actual_new_compact_detector_formula"]
    assert all(source.has(v) for v in (bumps.A1,bumps.A3,sp.Piecewise))
    assert detector.has(sp.Piecewise)
    assert not source.has(bumps.rho)
    assert not detector.has(sp.Symbol("normal_width",positive=True),sp.Symbol("cap_width",positive=True))
    assert all(value.has(sp.Integral) for value in d["unchanged_actual_bump_normalizers"])


def test_new_geometry_is_rechecked_at_larger_epsilon_not_old_upper():
    assert probes.EPSILON<probes.parent.A/100
    assert probes.RHO<probes.parent.H
    assert probes.parent.A**2<sp.Rational(3,4)*(probes.parent.S_MIN-probes.EPSILON)**2


def test_new_full_and_low_band_response_lowers_keep_strict_tail_budget():
    d=probes.data()
    D=probes.parent.SIGNAL_LOWER
    assert d["new_complete_remainder_pairing_upper"]==D/4
    assert d["positive_new_full_classical_pairing_lower"]==3*D/4
    assert d["new_source_only_high_spatial_tail_upper"]<D/4
    assert d["positive_new_low_spatial_band_pairing_lower"]==D/2


def test_rational_new_cutoff_has_claimed_decimal_upper_without_float():
    assert isinstance(probes.CUTOFF,sp.Rational)
    assert probes.CUTOFF<10**94
    assert probes.CUTOFF>10**93
    assert not probes.CUTOFF.has(sp.Float,sp.exp)


def test_new_compact_quantum_response_keeps_hbar_and_action_scale():
    d=probes.data()
    assert d["positive_new_compact_quantum_commutator_magnitude_lower"]==3*state.hbar*probes.parent.SIGNAL_LOWER/(4*state.kappa)
    assert d["finite_rational_width_and_cutoff_not_interacting_EFT_scales"] is True
    assert d["projected_source_not_compact_and_temporal_frequency_not_truncated"] is True


def test_all_invalid_exact_domains_are_rejected():
    d=audit.controls()
    assert d["groups"]["complex_jets"]["rejected_inputs"]==32
    assert d["groups"]["frequency_domain"]["rejected_inputs"]==34
    assert d["rejected_inputs"]==66


def test_every_new_exact_data_block_serializes_without_float_rounding():
    blocks=(jets.data(),basis.data(),transfer.data(),probes.data())
    serialized=json.loads(json.dumps(original.serialize(blocks)))
    assert len(serialized)==4
    assert serialized[3]["new_exact_normal_width"]==str(probes.EPSILON)
