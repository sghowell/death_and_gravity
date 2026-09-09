"""Scientific checks without the recursive parent report build."""
import json

import pytest
import sympy as sp
from p8_auxiliary_neighborhood import model
from p8_coupled_energy import (
    audit,
    bounds,
    energy,
    other_modes,
    scalars,
    tree,
    verify,
    window,
)
from p8_coupled_vertices import majorant


def test_all_exact_scalar_phase_and_other_mode_identities():
    rows=audit.residuals()
    assert len(rows)==32
    assert sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in rows.values())==71
    for value in rows.values():
        assert all(x==0 for x in value) if isinstance(value,sp.MatrixBase) else value==0


def test_all_continuous_and_finite_tree_gates():
    assert len(audit.gates())==51
    assert all(value is True for value in audit.gates().values())


def test_all_actual_time_jets_have_continuous_denominator_certificates():
    data=bounds.backgrounds()
    assert len(data["actual_jets"])==15
    for symbol,value in data["actual_jets"].items():
        numerator,denominator=sp.fraction(sp.cancel(value))
        num,den=sp.Poly(numerator,model.u),sp.Poly(denominator,model.u)
        assert den.nth(0)>0
        assert all(coefficient>=0 and not degree[0]%2 for degree,coefficient in den.terms())
        expected=sum(abs(c)*sp.Rational(1,2)**degree[0] for degree,c in num.terms())/den.nth(0)
        assert data["positive_even_denominator_proofs"][symbol]["absolute_half_interval_upper"]==expected


@pytest.mark.parametrize("name",("unitary","gamma"))
def test_actual_scalar_principal_gradient_and_margin_are_distinct(name):
    data=scalars.data(name)
    mapping=scalars.actual_background_jets()
    actual=scalars.checks()[name+"_actual_principal_gradient"]
    assert actual==sp.zeros(2)
    point=sp.Rational(1,2) if name=="unitary" else 0
    alpha0=data["alpha"].subs(scalars.z,0).subs(mapping,simultaneous=True).subs(model.u,point)
    gradient=data["principal_gradient"].subs(mapping,simultaneous=True).subs(model.u,point)
    assert sp.factor(alpha0[0,0]-gradient[0,0])>0
    assert sp.factor((alpha0-gradient).det())==0


@pytest.mark.parametrize("name,point",(
    ("unitary",sp.Rational(-1,2)),("unitary",sp.Rational(-1,4)),
    ("unitary",sp.Rational(1,4)),("unitary",sp.Rational(1,2)),
    ("gamma",sp.Rational(-1,4)),("gamma",0),("gamma",sp.Rational(1,4))))
def test_actual_finite_frequency_positive_matrices_at_exact_chart_anchors(name,point):
    mapping={key:sp.factor(value.subs(model.u,point)) for key,value in scalars.actual_background_jets().items()}
    mapping[scalars.z]=sp.Rational(1,10**12)
    for key,lower in (("alpha",sp.Rational(1,1000)),("principal_gradient",sp.Rational(1,1000)),
                      ("potential",sp.Rational(10**12,2000))):
        matrix=scalars.data(name)[key].subs(mapping,simultaneous=True).applyfunc(sp.factor)
        shifted=matrix-lower*sp.eye(2)
        assert shifted==shifted.T
        assert shifted[0,0]>0 and sp.factor(shifted.det())>0


@pytest.mark.parametrize("name",("unitary","gamma"))
def test_high_q_remainder_mixing_and_energy_budget(name):
    data=scalars.data(name)
    row=bounds.chart_bounds(name)
    assert data["beta_leading"]==data["beta_leading"].T
    assert data["antisymmetric_mixing"]==-data["antisymmetric_mixing"].T
    assert row["sufficient_q_lower"]>=10**8
    assert row["sufficient_q_lower"]>=2000*row["coefficient_bounds"]["potential_remainder"]["operator_upper"]
    assert row["energy_logarithmic_growth_upper"]<10**9
    assert row["kinetic_lower"]==sp.Rational(1,1000)
    assert row["full_potential_lower_over_q"]==sp.Rational(1,2000)


def test_gamma_finite_q_rank_one_increment_is_not_deleted():
    row=scalars.data("gamma")["alpha"]
    increment=scalars.clean(row-row.subs(scalars.z,0))
    assert increment!=sp.zeros(2)
    assert sp.factor(increment.det())==0
    assert scalars.checks()["gamma_finite_q_positive_rank_one_kinetic_increment"]==sp.zeros(2)


def test_fixed_comoving_derivative_includes_momentum_redshift():
    assert scalars.derivative(1/scalars.z)==-2*scalars.H/scalars.z
    assert scalars.derivative(scalars.H)==scalars.FIRST[0]


def test_full_antisymmetric_mixing_is_nonzero_on_actual_background():
    matrix=scalars.data("gamma")["antisymmetric_mixing"]
    actual=matrix.subs(scalars.actual_background_jets(),simultaneous=True).subs(
        {model.u:sp.Rational(1,5),scalars.z:sp.Rational(1,10**12)}).applyfunc(sp.factor)
    assert actual!=sp.zeros(2)


def test_canonical_initial_data_and_wrong_full_beta_shift():
    rows=energy.checks()
    assert rows["canonical_initial_column_Wronskian"]==sp.zeros(2)
    assert rows["canonical_initial_column_isotropy"]==sp.zeros(2)
    assert rows["initial_velocity_retains_antisymmetric_mixing"]==sp.zeros(2)
    assert energy.initial_mixing_negative_control()==sp.Matrix([[0,1],[-1,0]])


def test_two_chart_cover_and_every_explicit_column_margin():
    assert all(value>0 for value in bounds.elementary_margins().values())
    margins=window.data()["margins"]
    assert margins["energy_growth_exponent_at_most_half"]==0
    assert all(value>0 for name,value in margins.items() if name!="energy_growth_exponent_at_most_half")


def test_named_window_has_many_radians_and_high_q():
    d=window.data()
    assert d["full_window"]==sp.Rational(1,2*10**9)
    assert window.LOWER_K*window.HALF_WINDOW==250
    assert d["Hamiltonian_derivative_bound"]==8*10**13
    assert d["York_inverse_transfer_bound"]==sp.Rational(2,10**12)
    assert d["common_raw_phase_seed_after_center_Fourier_conversion"]==10**42


def test_wrong_chart_and_zero_frequency_denominators_fail_closed():
    s=scalars
    D=s.lam**2-(s.Je+s.ell**2*s.lam**2/2)*s.z
    assert bounds.coefficient_bound(1/D,"gamma")["absolute_upper"]==32
    with pytest.raises(ValueError):
        bounds.coefficient_bound(1/D,"unitary")
    with pytest.raises(ValueError):
        bounds.coefficient_bound(1/s.z,"gamma")


def test_exact_native_input_validation_before_cache():
    scalars.data("gamma")
    bounds.chart_bounds("gamma")
    assert audit.controls()["rejected_inputs"]==76


def test_actual_tensor_pump_and_wrong_formula_negative_control():
    u=model.u
    H=model.coefficients()["background"]["H"]
    actual=sp.factor(sp.Rational(3,2)*sp.diff(H,u)+sp.Rational(9,4)*H**2)
    assert sp.factor(actual-(6+30*u**2)/(1+u**2)**2)==0
    assert sp.factor((actual-(6+36*u**2)/(1+u**2)**2).subs(u,sp.Rational(1,2)))!=0
    assert other_modes.tensors()["global_I_energy_inflation_upper"]==3**10<10**5


def test_original_vector_mixing_and_all_polarization_weights():
    d=other_modes.vectors()
    assert d["mass"]==1000 and d["zeta"]==sp.Rational(1,10**6)
    assert d["actual_selected_mixing_coefficient"]==1813229+sp.Rational(5347035781757616,1000**4)
    assert 0<d["all_momentum_mixing_upper"]<1
    assert d["selected_state_unchanged"] is True
    assert d["reference_not_reselected_as_the_exact_state"] is True
    assert all(value==0 for value in other_modes.checks().values())


def test_both_remaining_mode_sectors_fit_raw_seed():
    assert other_modes.vectors()["raw_seed_margin"]>0
    assert other_modes.tensors()["raw_seed_margin"]>0
    assert tree.data()["propagating_mode_columns"]==7


def test_full_phase_homothety_for_momentum_dependent_vertices():
    assert len(tree.scale_checks())==3
    assert all(value==0 for value in tree.scale_checks().values())


def test_fixed_center_volume_is_not_lost_in_column_conversion():
    d=tree.data()
    raw=majorant.bound(window.RAW_SEED,8*window.UPPER_K,sp.Rational(2,window.LOWER_K))
    assert d["fixed_center_background_volume_upper"]==8
    assert d["cubic_mode_kernel_upper"]==8*720*sp.Rational(raw["cubic_labelled_kernel_upper"])
    assert d["quartic_mode_kernel_upper"]==8*720*sp.Rational(raw["quartic_labelled_kernel_upper"])


def test_exchange_species_orientation_time_order_and_Schur_constants():
    d=tree.data()
    assert d["Wick_majorant"]==720
    assert d["connected_exchange_contraction_majorant"]==84
    mu=7*(4*window.UPPER_K+1)**3
    assert d["momentum_species_measure_upper"]==mu
    assert d["fixed_total_momentum_Schur_upper"]==8*mu**3
    assert (8*mu**3)**2>(7*mu)**2


def test_full_contact_plus_exchange_not_contact_only():
    d=tree.data()
    T=window.data()["full_window"]
    Schur=d["fixed_total_momentum_Schur_upper"]
    assert d["quartic_connected_tree_numerator"]==Schur*(
        T*d["quartic_mode_kernel_upper"]+84*T*T*d["cubic_mode_kernel_upper"]**2)
    assert d["quartic_connected_tree_numerator"]>Schur*T*d["quartic_mode_kernel_upper"]


def test_named_sufficient_scale_not_an_old_example_transfer():
    d=tree.data()
    assert d["sufficient_named_M_tau"]==10**353
    assert d["cubic_block_bound_at_named_scale"]<=sp.Rational(1,1000)
    assert d["quartic_connected_tree_bound_at_named_scale"]<=sp.Rational(1,1000)
    previous=10**352
    assert d["cubic_transition_numerator"]/previous>tree.TARGET or (
        d["quartic_connected_tree_numerator"]/previous**2>tree.TARGET)
    assert d["not_the_previous_M_tau_10_to_24_example"] is True
    assert d["not_a_light_only_heavy_elimination_domain"] is True
    assert d["not_a_Wilsonian_or_all_orders_cutoff"] is True


def test_exact_bound_and_matrix_serialization_roundtrip():
    data=verify.serialize({"tree":tree.data(),"backgrounds":bounds.backgrounds(),
                           "scalar":scalars.data("gamma"),"other":other_modes.vectors()})
    assert json.loads(json.dumps(data))==data
    for value in (1.0,sp.Float(1),sp.oo,sp.nan):
        with pytest.raises((ValueError,TypeError)):
            verify.serialize(value)
