"""Native literal action, fresh solution, full boxes and central continuum tests."""
import json

import pytest
import sympy as sp
from p8_affine import verify as certificate
from p8_proca_retuned_margin import audit, central, domain, flow, model, original
from sympy.polys.domains import QQ
from sympy.polys.fields import field


def zero(value):
    return all(v==0 for v in value) if isinstance(value,sp.MatrixBase) else value==0


def test_exact_named_and_scalar_identity_counts():
    rows=audit.residuals()
    assert len(rows)==62
    assert sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in rows.values())==62
    assert all(zero(value) for value in rows.values())


@pytest.mark.parametrize("name",tuple(audit.residuals()))
def test_each_exact_identity(name):
    assert zero(audit.residuals()[name])


def test_every_actual_solution_and_continuum_gate():
    rows=audit.gates()
    assert len(rows)==54
    assert all(value is True for value in rows.values())


def test_new_total_margin_keeps_frozen_old_margin_and_adds_only_difference():
    assert model.OLD_MARGIN==sp.Rational(1,10**6)
    assert model.NEW_MARGIN==sp.Rational(1,200)
    assert model.INCREMENT==sp.Rational(4999,10**6)
    before=model.old.coefficients()
    after=model.coefficients()
    assert before is not after
    assert before["V0"]!=after["V0"]
    assert before["C"]!=after["C"]


def test_literal_density_keeps_volume_lapse_and_hamiltonian_sign():
    a=model.action()
    e=model.old.coefficients()["e"]
    X=a["physical_clock_X"]
    expected=model.N*e**3*a["added_physical_lower_scalar"].subs(X,model.N**-2)
    assert sp.factor(a["added_hat_volume_Lagrangian_density"]-expected)==0
    assert a["added_hat_volume_Hamiltonian_density"]==-a["added_hat_volume_Lagrangian_density"]
    assert a["clock_fixed_phase_lapse_Hessian_change"]==-8*model.INCREMENT/(1+model.u**2)**6


def test_new_rational_potential_has_full_conformal_fourth_power():
    a=model.action()
    c=model.coefficients()
    assert a["added_rationalized_potential"]==sp.factor(model.INCREMENT*(model.N**2-1)**2/(model.N*c["D"]*c["h"]))


def test_new_lapse_and_trace_velocity_both_differ_from_old_action():
    new,old=model.system(),model.old.system()
    assert new["constraint_matter_square"]!=old["constraint_matter_square"]
    assert new["on_constraint_rescaled_lapse_Hessian"]!=old["on_constraint_rescaled_lapse_Hessian"]
    assert new["lapse_flow"]!=-old["constraint_time_forcing"]/old["on_constraint_rescaled_lapse_Hessian"]


@pytest.mark.parametrize("denominator",(model.FIELD.ring.zero,model.FIELD.ring.gens[0],
    model.FIELD.ring.gens[1]-domain.box_arithmetic._qq(domain.N0)))
def test_uncancelled_pair_box_rejects_whole_box_poles(denominator):
    with pytest.raises(ValueError,match="denominator"):
        domain.pair_box(model.FIELD.ring.one,denominator)


def test_exact_box_rejects_foreign_coordinate_rings():
    other,x=field((sp.Symbol("foreign_coordinate"),),QQ)
    with pytest.raises(ValueError,match="field"):
        domain.rational_box(x)
    with pytest.raises(ValueError,match="ring"):
        domain.pair_box(other.ring.one,other.ring.one)


def test_uncancelled_fraction_enclosure_handles_negative_denominator():
    ring=model.FIELD.ring
    row=domain.pair_box(-2*ring.one,-3*ring.one)
    assert row["lower"]==row["upper"]==sp.Rational(2,3)


def test_new_real_and_complex_bounds_cover_entire_changed_system():
    rows=domain.enclosures()
    assert set(rows)=={"clock_physical_speed_squared","clock_physical_subluminal_margin",
        "on_constraint_rescaled_lapse_Hessian","constraint_matter_square","M",
        "lapse_flow","rescaled_trace_flow","hat_Hubble","conformal_log_flow"}
    assert rows["clock_physical_speed_squared"]["lower"]>sp.Rational(493,500)
    assert rows["clock_physical_speed_squared"]["upper"]<sp.Rational(247,250)
    assert domain.complex_annulus(rows["constraint_matter_square"])["lower"]>sp.Rational(9,1000)


def test_physical_acceleration_contains_all_three_chain_terms_and_proper_lapse():
    row=domain.acceleration()
    assert len(row["three_coordinate_time_chain_term_enclosures"])==3
    assert set(domain.physical_hubble_partials())=={"u","N","z"}
    assert row["physical_proper_time_acceleration_lower"]>sp.Rational(39997,10000)
    assert row["physical_proper_time_acceleration_upper"]<sp.Rational(4001,1000)
    assert row["three_coordinate_time_chain_term_enclosures"][1]!=(0,0)


def test_Picard_uses_new_velocities_and_weighted_phase_radii():
    d=flow.picard()
    assert d["whole_outer_complex_velocity_bounds"]==(sp.Rational(1,10000),9)
    assert d["Picard_contraction_upper"]==sp.Rational(3,25)
    assert d["actual_solution_displacement_bounds"]==(sp.Rational(1,10**11),sp.Rational(9,10**7))
    assert all(v>0 for v in d["Cauchy_phase_Lipschitz_entry_upper"])
    assert d["complex_conformal_D_minus_one_upper"]<sp.Rational(1,2)


def test_new_initial_charge_and_physical_scale_are_recomputed():
    d=flow.initial()
    assert d["new_positive_conserved_matter_charge_squared"]>0
    assert d["new_positive_conserved_matter_charge_squared"]!=d["old_charge_squared_at_same_lapse"]
    assert d["old_charge_residual_in_new_rescaled_constraint"]!=0
    assert d["central_physical_scale"]==sp.sqrt(domain.N0)
    assert d["central_physical_Hubble"]==d["central_lapse_velocity"]==0


def test_new_physical_principal_forms_use_real_new_lapse_and_conformal_bounds():
    d=flow.health()
    assert d["new_clock_kinetic_Schur_lower"]>11
    assert d["new_matter_kinetic_lower"]>sp.Rational(99,100)
    assert d["new_tensor_principal_lower"]>sp.Rational(99,100)
    assert d["new_matter_tensor_and_ordinary_Proca_principal_squared_speeds"]==(1,1,1)


BAD=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
    sp.Symbol("unfixed"),sp.sqrt(2),[])


@pytest.mark.parametrize("bad",BAD+(0,-1,2*domain.T))
def test_new_duration_domain_rejects_bad_inputs_after_valid_use(bad):
    assert flow.picard()["Picard_contraction_upper"]<1
    for call in (flow.duration,flow.picard):
        with pytest.raises((TypeError,ValueError)):
            call(bad)


@pytest.mark.parametrize("bad",BAD)
def test_central_diagnostic_rejects_bad_types_in_both_arguments(bad):
    for pair in ((bad,1),(model.NEW_MARGIN,bad)):
        with pytest.raises((TypeError,ValueError)):
            central.diagnostic(*pair)


@pytest.mark.parametrize("pair",((0,1),(-1,1),(1,1),(model.NEW_MARGIN,0),
    (model.NEW_MARGIN,-1),(model.NEW_MARGIN,2),(model.NEW_MARGIN,sp.Rational(1,2))))
def test_central_diagnostic_rejects_outside_margin_and_clock_tube(pair):
    with pytest.raises(ValueError):
        central.diagnostic(*pair)


def test_weaker_margin_is_an_actual_positive_matter_fast_bounce_control():
    weak=central.diagnostic(sp.Rational(1,10000))
    assert weak["positive_matter_charge_possible"] is True
    assert weak["positive_clock_kinetic_pivot"] is True
    assert weak["strict_classical_bounce_acceleration"] is True
    assert weak["clock_physical_squared_speed"]>1
    assert weak["physical_Hubble_proper_time_acceleration"]>3


def test_selected_margin_repairs_control_lapse_but_changes_its_charge():
    weak=central.diagnostic(sp.Rational(1,10000))
    strong=central.diagnostic()
    assert 0<strong["clock_physical_squared_speed"]<1
    assert strong["rescaled_matter_square"]>0
    assert strong["rescaled_matter_square"]!=weak["rescaled_matter_square"]


def test_negative_matter_square_is_flagged_not_assigned_a_real_state():
    d=central.diagnostic(lapse=sp.Rational(1001,1000))
    assert d["rescaled_matter_square"]<0
    assert d["positive_matter_charge_possible"] is False


def test_central_continuum_proof_uses_exact_monotonic_polynomials():
    d=central.theorem()
    assert d["matter_polynomial_derivative_lower_on_full_tube"]==sp.Rational(3223,9)
    assert d["matter_polynomial_positive_at_branch_cap"]==sp.Rational(2619,12500)
    assert d["positive_clock_squared_speed_lower"]==sp.Rational(49753,93129)
    assert d["positive_clock_squared_cone_margin_lower"]==sp.Rational(129954,485046875)
    assert d["positive_clock_squared_cone_margin_lower"]>sp.Rational(1,4000)


@pytest.mark.parametrize("name",("actual_acceleration_numerator_Bernstein_certificate",
    "actual_acceleration_denominator_Bernstein_certificate"))
def test_actual_Bernstein_coefficients_are_exact_positive_and_reconstruct(name):
    d=central.theorem()[name]
    assert len(d["coefficients"])==d["degree"]+1
    assert all(isinstance(v,sp.Rational) and v>0 for v in d["coefficients"])
    assert d["basis_reconstruction_residual"]==0


def test_new_central_acceleration_bound_is_actual_uniform_rational_lower():
    d=central.theorem()
    assert d["positive_actual_central_acceleration_lower"]==sp.Rational(207386958411035601,60630859375000000)
    assert d["positive_actual_central_acceleration_lower"]>3


def test_original_global_clock_and_relative_deformation_budget_keep_total_and_increment():
    d=original.data()
    assert d["positive_all_time_clock_squared_speed_lower"]==sp.Rational(1199,1215)
    assert d["relative_light_kinetic_increment_over_frozen_margin_upper"]==4*model.INCREMENT/(sp.Rational(1199,800)+4*model.OLD_MARGIN)
    assert d["added_scalar_on_real_clock_tube_upper"]==model.INCREMENT/100
    assert d["no_old_loop_potential_profile_or_quantum_state_transferred"] is True


def test_all_invalid_duration_and_central_domain_inputs_are_counted():
    d=audit.controls()
    assert d["groups"]["new_actual_solution_duration"]["rejected_inputs"]==34
    assert d["groups"]["central_family"]["rejected_inputs"]==35
    assert d["rejected_inputs"]==69


def test_all_new_data_serialize_without_float_rounding():
    blocks=(model.action(),flow.picard(),flow.initial(),flow.health(),original.data(),central.theorem())
    serialized=json.loads(json.dumps(certificate.serialize(blocks)))
    assert len(serialized)==6
    assert serialized[0]["new_total_margin"]=="1/200"
