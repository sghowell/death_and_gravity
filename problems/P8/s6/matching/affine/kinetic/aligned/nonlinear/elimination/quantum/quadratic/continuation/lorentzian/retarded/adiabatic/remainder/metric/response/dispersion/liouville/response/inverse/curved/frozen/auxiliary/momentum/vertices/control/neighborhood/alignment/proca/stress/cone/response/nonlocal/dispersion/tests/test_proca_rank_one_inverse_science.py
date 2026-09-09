"""Independent rank-one ordinary-Proca block and scalar inverse tests."""
import json

import pytest
import sympy as sp
from p8_affine import verify as original
from p8_proca_rank_one_inverse import audit, cut, kernel, spectral


def test_all_new_exact_residuals():
    rows=audit.residuals()
    assert len(rows)==50
    assert sum(len(value) if isinstance(value,sp.MatrixBase) else 1 for value in rows.values())==67
    for value in rows.values():
        assert all(x==0 for x in value) if isinstance(value,sp.MatrixBase) else value==0


def test_all_analytic_domain_and_scope_gates():
    assert len(audit.gates())==21
    assert all(value is True for value in audit.gates().values())


@pytest.mark.parametrize("sector,coefficient",(("T",1-spectral.z),("L",1+spectral.z)))
def test_actual_minimal_pair_vertices_have_no_lapse_channel(sector,coefficient):
    assert spectral.pair(sector)==sp.Matrix([0,2*coefficient])
    assert spectral.checks()[sector+"_new_pair_from_actual_minimal_Hamiltonian"]==sp.zeros(2,1)


@pytest.mark.parametrize("sector",("T","L"))
def test_actual_canonical_contact_is_retained_and_not_called_zero(sector):
    data=spectral.contact(sector)
    assert data["actual_fixed_canonical_contact"]!=sp.zeros(2)
    assert data["static_contact_identity"]==sp.zeros(2)
    assert (data["actual_fixed_canonical_contact"]+data["zero_frequency_bubble"]
            -data["static_vacuum_Hessian"]).applyfunc(sp.factor)==sp.zeros(2)


@pytest.mark.parametrize("sector",("T","L"))
def test_new_flat_adiabatic_frequency_taylor_coefficients(sector):
    rows=spectral.checks()
    for output in ("N","Z"):
        for source in range(2):
            for order in (2,4):
                assert rows[sector+f"_new_adiabatic_frequency_Taylor_{output}_{source}_{order}"]==0


def test_rank_one_gram_and_strict_positive_weight():
    data=spectral.data()
    assert data["actual_physical_pair_Gram_matrix"].det()==0
    assert sp.factor(data["positive_pair_polynomial"]-3*(spectral.z-sp.Rational(1,3))**2)==sp.Rational(8,3)
    assert spectral.checks()["actual_new_curved_local_fourth_coefficient"]==sp.zeros(2)


def test_massive_closed_form_and_three_subtractions_are_independent_checks():
    rows=spectral.checks()
    for key in ("new_three_moment_closed_scalar_block","same_exact_three_subtraction_identity",
                "same_three_subtracted_denominator_and_radial_measure"):
        assert rows[key]==0


@pytest.mark.parametrize("order",(1,2,3,4))
def test_low_frequency_integral_coefficients_by_finite_polynomial_sum(order):
    weight=sp.Poly(spectral.data()["positive_radial_weight"]*(1-spectral.y**2)**(order-1),spectral.y)
    moment=sum(coefficient/sp.Integer(power[0]+1) for power,coefficient in weight.terms())
    assert spectral.data()["low_frequency_H_coefficients_in_p_over_mass_squared"][order]==(-1)**(order-1)*moment/4**order


def test_threshold_gap_not_inferred_from_high_frequency_log():
    data=spectral.data()
    assert data["H_at_zero"]==4
    assert data["H_at_threshold"]==sp.Rational(16,15)
    assert data["H_slope_at_zero"]==sp.Rational(9,35)/spectral.m**2
    assert spectral.checks()["threshold_value_by_independent_radial_integral"]==0
    assert spectral.checks()["spurious_large_frequency_truncation_root"]==0


@pytest.mark.parametrize("ratio",(0,1,sp.Rational(9,2)))
def test_positive_axis_uses_full_massive_formula_and_removable_zero(ratio):
    data=spectral.positive_axis(ratio)
    assert data["H_continuous_lower"]==4
    assert data["range_inverse_absolute_upper"]==sp.Rational(1,4)
    if ratio==0:
        assert data["H"]==4
        assert data["range_inverse"]==-sp.Rational(1,4)


def test_half_plane_sign_and_real_monotonicity_exact_integrands():
    rows=spectral.checks()
    assert rows["strict_half_plane_imaginary_sign_integrand"]==0
    assert rows["strict_real_gap_monotonicity_integrand"]==0


@pytest.mark.parametrize("fraction",(sp.Rational(1,10),sp.Rational(1,2),sp.Rational(9,10)))
def test_open_cut_density_keeps_exact_positive_denominator(fraction):
    data=cut.fraction(fraction)
    D,U=data["real_positive_H_on_upper_cut"],data["positive_imaginary_H_over_pi"]
    rho=data["positive_range_inverse_cut_density"]
    assert U>0
    assert data["real_cut_H_continuous_lower"]>sp.Rational(16,15)
    assert sp.simplify(sp.expand(rho*(D**2+sp.pi**2*U**2)-U))==0


def test_continuous_cut_gap_and_exact_endpoint_constants():
    data=cut.data()
    assert data["uniform_first_sheet_H_real_part_lower"]==sp.Rational(16,15)
    assert data["range_inverse_absolute_upper_on_first_sheet"]==sp.Rational(15,16)
    assert data["density_over_sqrt_z_at_threshold"]==sp.Rational(675,512)
    assert data["density_times_log_tau_over_mass_squared_squared_at_infinity"]==sp.Rational(1,2)
    assert cut.checks()["positive_atanh_remainder_derivative"]==0


def test_no_hidden_pole_or_instantaneous_scalar_inverse_remainder():
    assert cut.data()["no_instantaneous_range_inverse"] is True
    assert cut.checks()["new_inverse_high_frequency_limit"]==0
    assert kernel.data()["instantaneous_inverse"]==0
    assert kernel.checks()["new_no_hidden_high_frequency_instantaneous_term"]==0


def test_exact_spectral_moments_and_primitive_normalization():
    assert cut.data()["spectral_static_moment"]==sp.Rational(1,4)
    assert cut.data()["spectral_inverse_squared_frequency_moment"]==sp.Rational(9,560)/spectral.m**2
    assert kernel.data()["primitive_lower"]==-sp.Rational(1,2)
    assert kernel.data()["primitive_upper"]==0
    assert kernel.checks()["new_kernel_primitive_normalization"]==0


def test_retarded_kernel_L1_claim_has_no_invented_numeric_L1_bound():
    data=kernel.data()
    assert data["L1_half_line_kernel"] is True
    assert data["short_window_C0_inverse_norm_tends_to_zero"] is True
    assert data["no_numerical_L1_norm_claim"] is True
    assert kernel.checks()["new_causal_sine_Laplace_normalization"]==0


def test_positive_mass_scaling_is_not_a_massless_limit():
    one,two=kernel.moment_bound(1),kernel.moment_bound(2)
    assert one["primitive_C0_upper"]==two["primitive_C0_upper"]==sp.Rational(1,2)
    assert one["inverse_squared_frequency_spectral_moment"]==4*two["inverse_squared_frequency_spectral_moment"]
    with pytest.raises(ValueError):
        kernel.mass(0)


def test_scalar_range_inverse_does_not_invert_full_two_source_block():
    assert cut.data()["full_two_source_inverse"] is False
    B=sp.Symbol("nonzero_block",nonzero=True)
    assert sp.diag(0,B)*sp.diag(0,1/B)==sp.diag(0,1)
    assert sp.diag(0,1)!=sp.eye(2)


def test_native_input_guards_are_checked_after_caches_are_populated():
    spectral.pair("T")
    spectral.contact("T")
    cut.fraction(sp.Rational(1,2))
    kernel.moment_bound(1)
    assert audit.controls()["rejected_inputs"]==78


def test_exact_report_payload_roundtrip_and_inexact_rejection():
    data=original.serialize({"block":spectral.data(),"cut":cut.data(),"kernel":kernel.data()})
    assert json.loads(json.dumps(data))==data
    for value in (1.0,sp.Float(1),sp.oo,sp.nan):
        with pytest.raises((TypeError,ValueError)):
            original.serialize(value)
