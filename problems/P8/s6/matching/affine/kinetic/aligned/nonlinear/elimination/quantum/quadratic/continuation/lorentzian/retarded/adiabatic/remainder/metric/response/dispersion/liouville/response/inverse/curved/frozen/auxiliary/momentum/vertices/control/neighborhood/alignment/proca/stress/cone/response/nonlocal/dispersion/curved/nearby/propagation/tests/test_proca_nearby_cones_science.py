"""Independent full geometric, time-dependent principal and whole-box tests."""
import json

import pytest
import sympy as sp
from p8_affine import verify as original
from p8_proca_nearby_cones import audit, bounds, geometry, principal, rational


def test_all_exact_identities_and_scalar_counts():
    rows=audit.residuals()
    assert len(rows)==61
    assert sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in rows.values())==70
    assert all(v==(sp.zeros(*v.shape) if isinstance(v,sp.MatrixBase) else 0) for v in rows.values())


def test_all_whole_box_and_solution_bridge_gates():
    assert len(audit.gates())==21
    assert all(value is True for value in audit.gates().values())


@pytest.mark.parametrize("left,right",[(left,right) for i,left in enumerate(geometry.CHANNELS)
                                       for right in geometry.CHANNELS[i:]])
def test_every_full_geometric_scalar_phase_pair(left,right):
    row=geometry.pair(left,right)
    assert row["residual"]==0
    assert sp.factor(row["actual"]-row["expected"])==0


def test_omitting_moving_volume_trace_boundary_changes_offcenter_quadratic():
    correction=9*geometry.H*geometry.p0
    assert correction!=0
    row=geometry.pair("v","v")
    assert sp.factor(row["actual"]-correction-row["expected"])==-sp.factor(correction)


def test_freezing_the_potential_is_not_the_Euler_gradient():
    data,form=principal.data(),principal.formula()
    difference=sp.factor(form["Euler_gradient"][0,0]-data["potential_over_q"][0,0])
    assert sp.factor(difference-form["B_dot"]-principal.H*form["B"])==0
    assert difference!=0


def test_time_dependent_momentum_chart_retains_canonical_boundary():
    assert principal.data()["Hamiltonian"].coeff(principal.H)==-principal.b*principal.Pb
    assert principal.formula()["division_by_Hubble_or_gamma_crossing_coefficient"] is False


def test_finite_q_Legendre_zero_is_not_hidden_by_high_q_statement():
    determinant=principal.data()["finite_q_Legendre_determinant"]
    sample={principal.beta:0,principal.r:-1,principal.h:-1,principal.m:1,principal.rn:1}
    assert determinant.subs(sample).subs(principal.q,sp.Rational(1,4))==0
    assert principal.formula()["all_finite_q_Legendre_regular"] is False


def test_principal_limit_rejects_a_divergent_declared_limit():
    with pytest.raises(ValueError,match="diverges"):
        principal.infinity_limit(principal.q)
    assert principal.infinity_limit(1/principal.q)==0


def test_physical_characteristic_basis_has_luminal_matter_and_clock_Schur():
    p=principal.physical()
    assert p["physical_matter_speed_squared"]==1
    assert p["diagonal_kinetic"][0,1]==p["diagonal_gradient"][0,1]==0
    assert sp.factor(p["diagonal_kinetic"][0,0]+principal.h/(4*p["matter_square_N"]**2))==0


@pytest.mark.parametrize("name",("d0","V0","I0","b0","A","M","B","C"))
def test_actual_shifted_coefficient_is_rational_and_has_no_primitive(name):
    value=rational.coefficients()[name]
    assert value.is_rational_function(rational.u,rational.N)
    assert not value.has(rational.model.I,rational.model.Iphi)


def test_original_primitive_endpoint_and_time_derivative_are_retained():
    checks=rational.checks()
    assert checks["actual_primitive_endpoint_lapse_derivative"]==0
    assert checks["actual_potential_retains_original_primitive_time_derivative"]==0
    assert rational.coefficients()["I0"]!=0


def test_constraint_density_elimination_is_actual_algebra():
    d=rational.system()
    z=rational.zf
    assert d["A"]*z*z+d["B"]*z+d["M"]*d["constraint_matter_square"]+d["C"]==0


def test_exact_central_speed_control_and_original_clock_are_distinct():
    value=rational.system()["clock_physical_speed_squared"]
    clock=value.subs([(rational.uf,0),(rational.Nf,1),(rational.zf,0)]).as_expr()
    shifted=value.subs([(rational.uf,0),(rational.Nf,rational.QQ(1000001,1000000)),(rational.zf,0)]).as_expr()
    assert clock==sp.Rational(749375,749377)<1
    assert shifted>1
    assert rational.checks()["actual_central_clock_speed_reproduces_full_Euler_result"]==0


def test_actual_lapse_and_flow_fixed_phase_chain_identities():
    assert len(rational.chain_checks())==5
    assert all(value==0 for value in rational.chain_checks().values())


def test_full_speed_parity_is_not_only_a_central_jet():
    speed=rational.system()["clock_physical_speed_squared"]
    for poly in (speed.numer,speed.denom):
        assert all((powers[0]+powers[2])%2==0 for powers in poly)


def test_exact_polynomial_signatures_pin_all_nonzero_terms():
    signatures=rational.signatures()
    speed=signatures["clock_physical_speed_squared"]
    assert speed["numerator"]["monomials"]==1457
    assert speed["denominator"]["monomials"]==1387
    poly=rational.system()["clock_physical_speed_squared"].numer
    assert rational.polynomial_signature(poly+1)!=rational.polynomial_signature(poly)
    assert len(speed["numerator"]["exact_sorted_monomial_rational_coefficient_sha256"])==64


@pytest.mark.parametrize("name",("clock_physical_speed_squared_excess",
    "on_constraint_rescaled_lapse_Hessian","constraint_matter_square","M"))
def test_every_exact_box_has_a_separated_denominator_and_nontrivial_width(name):
    box=bounds.enclosures()[name]
    assert box["lower"]<box["upper"]
    assert all(isinstance(box[key],sp.Rational) for key in ("lower","upper"))
    assert box["denominator"]["lower"]*box["denominator"]["upper"]>0
    assert box["numerator"]["shifted_monomials"]>0


def test_clock_speed_not_just_speed_squared_is_uniformly_separated():
    box=bounds.enclosures()["clock_physical_speed_squared_excess"]
    assert box["lower"]>(1+bounds.CLOCK_SPEED_EXCESS_LOWER)**2-1
    assert box["upper"]<(1+bounds.CLOCK_SPEED_EXCESS_UPPER)**2-1
    assert box["lower"]>sp.Rational(5,10**6)
    assert box["upper"]<sp.Rational(16,10**6)


def test_actual_phase_and_odd_boundary_primitive_fit_entire_box():
    b=bounds.solution_bridge()
    assert b["actual_phase_displacement_upper"]==sp.Rational(1,10**6)
    assert b["actual_trace_linear_absolute_upper"]==sp.Rational(100000,9999999999)
    assert b["rescaled_trace_absolute_upper"]<bounds.TRACE_RADIUS
    assert b["actual_lapse_lower"]>1
    assert b["actual_eomega_upper"]**2>b["actual_lapse_upper"]


def test_all_scalar_kinetic_and_tensor_lower_bounds_are_positive():
    b=bounds.solution_bridge()
    assert b["clock_diagonal_kinetic_lower"]>7
    assert b["matter_diagonal_kinetic_lower"]>sp.Rational(99,100)
    assert b["tensor_principal_coefficient_lower"]>sp.Rational(99,100)
    assert b["actual_scalar_and_tensor_gradients_positive"] is True


def test_same_emission_characteristics_have_exact_positive_separation():
    value=bounds.characteristic_separation()
    assert value["classical_comoving_characteristic_separation_lower"]==sp.Rational(99,2*10**14)
    assert value["classical_comoving_characteristic_separation_upper"]==sp.Rational(202,125*10**12)
    assert value["same_emission_event_and_same_final_clock_slice"] is True
    assert value["finite_frequency_wavepacket_error_or_cutoff_established"] is False
    assert value["not_a_UV_exclusion"] is True


def test_shorter_characteristic_interval_has_exact_linear_scaling():
    whole=bounds.characteristic_separation()
    half=bounds.characteristic_separation(bounds.TIME_WINDOW/2)
    assert 2*half["classical_comoving_characteristic_separation_lower"]==whole["classical_comoving_characteristic_separation_lower"]
    assert 2*half["classical_comoving_characteristic_separation_upper"]==whole["classical_comoving_characteristic_separation_upper"]


def test_native_guards_reject_all_bad_calls_after_valid_cache_population():
    bounds.characteristic_separation()
    geometry.pair("v","p")
    assert audit.controls()["rejected_inputs"]==49


def test_symbolic_and_exact_bound_report_data_roundtrips_without_float():
    data=original.serialize({"coefficients":rational.coefficients(),"signatures":rational.signatures(),
        "enclosures":bounds.enclosures(),"bridge":bounds.solution_bridge(),"principal":principal.formula()})
    assert json.loads(json.dumps(data))==data
    for value in (1.0,sp.Float(1),sp.oo,sp.nan):
        with pytest.raises((TypeError,ValueError)):
            original.serialize(value)
