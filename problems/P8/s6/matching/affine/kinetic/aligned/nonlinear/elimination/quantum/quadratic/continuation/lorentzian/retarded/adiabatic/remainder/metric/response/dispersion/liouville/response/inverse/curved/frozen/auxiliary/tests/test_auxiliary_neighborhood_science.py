"""Exact domain, source, boundary and actual Hamiltonian controls."""
import pytest
import sympy as sp
from p8_auxiliary_neighborhood import bounds, center, model, verify


def test_all_exact_anchors_and_continuous_gates():
    verify.affine.certify_residuals(verify.residuals())
    assert all(value is True for value in bounds.gates().values())


def test_all_nine_spatial_canonical_invariants_retained():
    g=model.generic()
    assert len(model.VARIABLES)==9
    for field in model.VARIABLES:
        assert sp.diff(g["Hamiltonian_over_hat_volume"],field)!=0
    assert len(g["pieces"])==8


def test_electric_magnetic_normalizations_are_opposite_in_zeta():
    g=model.generic()
    zeta=sp.Symbol("zeta",positive=True)
    E,B=sp.symbols("E B")
    assert g["pieces"]["electric"].subs(model.electric,E/zeta)==E/(2*zeta*g["eo"])
    assert g["pieces"]["magnetic"].subs(model.magnetic,zeta*B)==zeta*B/(4*g["eo"])


def test_full_scalar_complex_majorant_uses_all_33_monomials():
    data=bounds.scalar_potential_bound()
    assert data["monomials"]==33
    assert data["upper"]==sp.Rational(
        101704307704782433883533238264601170673664,
        422661243303360596017265796661376953125)
    assert data["upper"]<1000


def test_complex_lapse_image_and_coefficient_branch_discs():
    data=bounds.coefficients()
    assert data["s_minus_one_absolute_upper"]==sp.Rational(1,99)
    assert data["x_plus_one_absolute_upper"]==sp.Rational(199,9801)
    assert data["ratio_minus_one_absolute_upper"]<sp.Rational(3,100)
    assert data["ratio_minus_one_absolute_upper"]/2<sp.Rational(1,50)


def test_complex_retuned_mass_denominators_are_bounded():
    data=bounds.coefficients()
    assert data["gamma_t_denominator_absolute_lower"]>6
    assert data["gamma_t_absolute_upper"]<2
    assert data["gamma_s_absolute_upper"]<2
    assert data["inverse_gamma_s_absolute_upper"]<2


def test_actual_source_ode_complex_bound():
    data=bounds.coefficients()
    assert data["Q_ODE_coefficient_absolute_upper"]<sp.Rational(3,2)
    assert data["Q_ODE_forcing_over_distance_upper"]<7
    assert data["Q_over_distance_squared_upper"]==sp.Rational(140,37)<4


def test_original_boundary_primitive_and_phi_bounds():
    data=bounds.coefficients()
    assert data["boundary_primitive_absolute_upper"]==sp.Rational(2,99)
    assert data["boundary_primitive_phi_absolute_upper"]==sp.Rational(100,99)
    assert data["boundary_primitive_s_absolute_upper"]<2


def test_complete_Hamiltonian_modulus_bound():
    data=bounds.coefficients()
    assert sum(data["full_Hamiltonian_piece_upper_before_N"].values())==sp.Rational(6361,2)
    assert data["full_Hamiltonian_absolute_upper_from_pieces"]==6361
    assert data["full_Hamiltonian_absolute_upper_from_pieces"]<10000
    assert data["full_post_boundary_potential_absolute_upper"]<3000


def test_cauchy_derivatives_keep_factorials_and_nine_directions():
    data=bounds.contraction()
    assert data["H_NY_absolute_upper"]==4000000
    assert data["H_NNN_absolute_upper"]==480000000000
    assert data["H_NNY_absolute_upper"]==1600000000
    assert data["lapse_force_at_center_upper"]==9*4000000*bounds.JET_RADIUS


def test_newton_rate_and_actual_disc_image():
    data=bounds.contraction()
    assert data["Newton_contraction_upper"]==sp.Rational(3000000000009,31250000000000)
    assert data["Newton_contraction_upper"]<sp.Rational(1,10)
    assert data["Newton_closed_disc_image_radius_upper"]<bounds.ROOT_RADIUS


def test_actual_lapse_displacement_and_secondary_pivot():
    data=bounds.contraction()
    assert data["lapse_displacement_over_input_radius_upper"]==800000000
    assert data["lapse_displacement_upper"]==sp.Rational(1,1250000000000000)
    assert sp.Rational(1,20)-data["lapse_pivot_change_upper"]>sp.Rational(1,25)


def test_temporal_reconstruction_uses_quadratic_source():
    data=bounds.contraction()
    assert data["temporal_vector_over_input_radius_upper"]==sp.Rational(502168640000243,125000000000000)
    assert data["temporal_vector_over_input_radius_upper"]<5
    g=model.generic()
    assert sp.factor(g["temporal_reconstruction"]-g["d"]*g["trace_reconstruction"]-g["c"]+g["gt"]*model.j/g["U"])==0


def test_closed_bounce_boundary_primitive_is_not_dropped():
    data=center.data()
    actual=data["actual_boundary_primitive_phi_at_bounce"]
    assert actual.subs(center.r,1)==0
    assert actual.subs(center.r,4)!=0
    assert sp.diff(actual,center.r,2).subs(center.r,1)==-18


def test_full_closed_bounce_lapse_derivatives():
    assert center.data()["lapse_background_derivatives"]==(
        sp.Rational(801,100),0,-sp.Rational(749377,250000),
        sp.Rational(1077891,31250),-sp.Rational(31494387,200000))


def test_bounce_invariant_decomposition_replays_all_terms():
    data=center.data()
    assert data["full_closed_bounce_Hamiltonian"]==data["background"]+data["linear_invariant_part"]+data["quadratic_invariant_part"]
    assert sp.Poly(data["linear_invariant_part"],*model.VARIABLES).total_degree()==1
    assert sp.Poly(data["quadratic_invariant_part"],*model.VARIABLES).total_degree()==2
    assert center.series()["lapse"][0]!=0


def test_bounce_stationary_first_order_satisfies_actual_constraint():
    data=center.data()
    first=center.series()["lapse"][0]
    assert sp.factor(data["lapse_background_derivatives"][2]*first+
                     sp.diff(data["linear_invariant_part"],model.N).subs(model.N,1))==0


def test_multivariate_Taylor_counts_and_remainder():
    assert bounds.taylor_majorant(3)["monomial_count"]==165
    assert bounds.taylor_majorant(4)["monomial_count"]==495
    data=bounds.contraction()
    assert data["quartic_invariant_Taylor_remainder_input_radius"]==sp.Rational(1,10**27)
    assert data["quartic_invariant_Taylor_remainder_upper"]==sp.Rational(3861,299300000000)
    assert data["quartic_invariant_Taylor_remainder_upper"]<sp.Rational(1,10**7)


def test_actual_radius_function_and_zero_source():
    assert bounds.response_radius(0)["lapse_displacement_upper"]==0
    assert bounds.response_radius(bounds.JET_RADIUS)["lapse_displacement_upper"]==8*sp.Rational(1,10**16)
    assert bounds.response_radius(bounds.JET_RADIUS/1000)["temporal_vector_absolute_upper"]==5*sp.Rational(1,10**27)


def test_invalid_types_ranges_and_float_serializer_rejected():
    assert verify.controls()["rejected_inputs"]==28
    for value in (True,0.0,sp.Float(0),"0"):
        with pytest.raises(TypeError):
            bounds.response_radius(value)
    with pytest.raises(ValueError):
        bounds.response_radius(2*bounds.JET_RADIUS)


def test_no_cutoff_or_nonlocal_constraint_is_inferred():
    result=bounds.response_radius(bounds.JET_RADIUS)
    assert "not a frequency or quantum-response estimate" in result["boundaries"]
    assert bounds.gates()["classical_tree_plus_margin_not_the_nonlocal_quantum_constraint"]
