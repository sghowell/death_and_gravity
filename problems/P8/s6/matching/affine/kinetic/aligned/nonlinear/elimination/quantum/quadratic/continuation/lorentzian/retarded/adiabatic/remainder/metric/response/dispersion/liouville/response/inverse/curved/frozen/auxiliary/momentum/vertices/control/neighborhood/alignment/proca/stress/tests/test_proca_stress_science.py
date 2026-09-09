"""Actual ordinary Proca stress and explicitly fixed new-profile tests."""
import json

import pytest
import sympy as sp
from p8_affine import verify as original
from p8_clock_tadpole import window
from p8_constant_proca import quantum
from p8_proca_stress import audit, estimates, profiles, readouts, subtraction
from p8_vector_regularity import estimates as previous
from p8_vector_regularity import spectral
from p8_vector_state import wkb


def test_all_exact_ordinary_stress_residuals():
    rows=audit.residuals()
    assert len(rows)==170
    assert sum(len(v) if isinstance(v,sp.MatrixBase) else 1 for v in rows.values())==186
    for value in rows.values():
        assert all(x==0 for x in value) if isinstance(value,sp.MatrixBase) else value==0


def test_all_continuous_and_scope_gates():
    assert len(audit.gates())==46
    assert all(value is True for value in audit.gates().values())


@pytest.mark.parametrize("kind",("transverse","longitudinal"))
def test_actual_new_unit_energy_weights_and_canonical_readout(kind):
    assert readouts.weights(kind)["A"]==readouts.weights(kind)["B"]==1
    assert readouts.rows(kind,0)==sp.Matrix([[readouts.omega2,0,1]])
    assert quantum.matrices()[kind]["canonical_energy_change"]!=sp.zeros(2)


@pytest.mark.parametrize("kind",("transverse","longitudinal"))
def test_exact_mode_Ward_identity_with_unchanged_pressure(kind):
    H=wkb.background()["H"]
    assert spectral.clean(readouts.rows(kind,1)+3*H*(
        readouts.rows(kind,0)+spectral.rows(kind+"_pressure",0)))==sp.zeros(1,3)
    assert quantum.matrices()[kind]["canonical_pressure_change"]==sp.zeros(2)


@pytest.mark.parametrize("order",range(6))
def test_both_actual_derivative_rows_have_continuous_frequency_bounds(order):
    for kind in ("transverse","longitudinal"):
        row=readouts.row_bounds(kind,order)
        assert row["all_frequency_powers_at_most_derivative_order_plus_one"] is True
        assert all(v==0 for v in row["box_reconstructions"])
        assert row["normalized_reference_product_envelope"]>0


def test_actual_ordinary_adiabatic_coefficients():
    for kind in ("transverse","longitudinal"):
        c=readouts.weights(kind)["c1"]
        d=wkb.frequency(kind)
        ad=subtraction.adiabatic_rows(kind,0)
        assert ad[0]==2
        assert sp.factor(ad[1]-c*c)==0
        assert sp.factor(ad[2]-d["P2"]**2-c*d["B2"]+c*c*d["P2"])==0


@pytest.mark.parametrize("order",range(6))
def test_all_new_subtraction_tails_are_integrable(order):
    for kind in ("transverse","longitudinal"):
        result=subtraction.reference_tail(kind,order)
        assert set(result["low_coefficient_residuals"])==set(range(5))
        assert all(value==0 for value in result["low_coefficient_residuals"].values())
        assert result["integrable_tail_certified"] is True
        assert result["no_lower_Laurent_power"] is True
        assert result["all_majorant_coefficients_nonnegative"] is True


def test_new_polarization_specific_integer_tail_bounds():
    expected={
        "transverse":[40976759,330164332,5423162383,407150071277,20454860860415,699227829249482],
        "longitudinal":[39243120,343842778,6026341152,423718215792,25184980639455,1131096194226166]}
    for kind,values in expected.items():
        assert [subtraction.reference_tail(kind,j)["reference_bracket_tail_integer_upper"]
                for j in range(6)]==values


def test_lower_reference_failure_is_actual_new_fourth_derivative():
    result=subtraction.reference_tail("longitudinal",4,2)
    assert result["integrable_tail_certified"] is False
    assert result["low_coefficient_residuals"][4].subs({readouts.u:0,readouts.z:1})==23808


def test_actual_new_local_coefficients_differentiated_before_bounding():
    for order in range(6):
        d=estimates.local_derivatives(order)
        for n,row in quantum.local_coefficients().items():
            assert sp.factor(d["actual_new_local_energy_derivatives"][n]
                             -sp.diff(row["energy"],wkb.u,order))==0
        assert all(v==0 for v in d["reconstructions"].values())


def test_old_nonminimal_local_energy_is_not_silently_reused():
    old=previous.local_derivatives("energy",0)["coefficients"]
    new=estimates.local_derivatives(0)["actual_new_local_energy_derivatives"]
    assert any(sp.factor(old[n]-new[n])!=0 for n in range(3))


def test_actual_fixed_state_radial_envelope_and_pressure_transfer():
    d=estimates.physical_bounds(estimates.SCALE)
    old=previous.physical_bounds(estimates.SCALE,estimates.MASS)
    assert d["unchanged_initial_state_radial_envelope"]==previous.radial_envelope(estimates.MASS)
    assert d["identical_pressure_derivative_bounds"]==old["normalized_derivative_bounds"]["pressure"]
    assert d["fixed_canonical_mass_time_product"]==1000
    assert d["M_tau"]==10**400


@pytest.mark.parametrize("order",range(6))
def test_all_new_energy_and_pressure_derivatives_and_changes(order):
    d=estimates.physical_bounds(estimates.SCALE)
    energy=d["new_energy_derivative_bounds"][order]
    assert energy["total"]==sum(value for key,value in energy.items() if key!="total")
    assert 0<energy["total"]<sp.Rational(1,10**770)
    assert 0<d["identical_pressure_derivative_bounds"][order]["total"]<sp.Rational(1,10**770)
    assert energy["total"]<=d["new_minus_old_energy_derivative_upper"][order]<sp.Rational(1,10**770)


def test_exact_scale_homogeneity_for_fixed_mass_and_state():
    a,b=estimates.physical_bounds(10**400),estimates.physical_bounds(2*10**400)
    for j in range(6):
        assert b["new_energy_derivative_bounds"][j]["total"]==a["new_energy_derivative_bounds"][j]["total"]/4
        assert b["new_minus_old_energy_derivative_upper"][j]==a["new_minus_old_energy_derivative_upper"][j]/4


def test_continuum_integrands_use_two_transverse_and_one_longitudinal():
    d=profiles.integrands()
    rows=d["new_energy_modes"]
    assert sp.expand(d["new_energy_combined_finite_integrand"]
        -2*(rows["transverse"]["physical_new_energy"]-rows["transverse"]["actual_fourth_order_subtraction"])
        -(rows["longitudinal"]["physical_new_energy"]-rows["longitudinal"]["actual_fourth_order_subtraction"]))==0
    assert d["finite_energy_local_coefficients"]==quantum.local_coefficients()


def test_literal_fixed_profile_replacement_and_clock_norm_jets():
    d=profiles.action()
    assert sp.expand(d["new_clock_tube_profile"]-d["old_clock_tube_profile"]
                     -(d["new_fixed_energy"]-d["old_fixed_energy"])*(profiles.x+1)/2)==0
    assert sp.diff(d["new_clock_tube_profile"],profiles.x,2)==0
    assert d["explicit_new_minus_old_profile"]!=0


def test_actual_one_point_and_clock_equation_identities():
    rows=profiles.checks()
    for name in ("new_fixed_profile_cancels_actual_energy",
                 "new_fixed_profile_cancels_actual_pressure",
                 "new_fixed_profile_clock_equation_is_new_conserved_Ward_balance"):
        assert rows[name]==0


def test_old_profile_not_misidentified_as_new_energy_cancellation():
    d=profiles.action()
    old=d["old_clock_tube_profile"]
    energy=(2*profiles.x*sp.diff(old,profiles.x)-old).subs(profiles.x,-1)
    assert sp.expand(energy+d["new_fixed_energy"])==d["new_fixed_energy"]-d["old_fixed_energy"]
    assert sp.expand(energy+d["new_fixed_energy"])!=0


def test_global_window_preserves_open_vacuum_norm_neighborhood():
    d=profiles.action()
    assert window.cutoff(sp.Integer(-1))==1
    assert d["new_global_profile"].subs(profiles.x,0)==0
    for j in range(1,6):
        assert sp.diff(d["new_global_profile"],profiles.x,j).subs(profiles.x,0)==0


def test_all_twenty_one_global_mixed_bounds_and_changes():
    d=profiles.profile_bounds(estimates.SCALE)
    for name in ("new_profile_global_mixed_derivative_bounds",
                 "profile_change_global_mixed_derivative_bounds"):
        assert len(d[name])==21
        assert all(0<value<sp.Rational(1,10**770) for value in d[name].values())
    assert d["fixed_profiles_not_recomputed_under_metric_or_state_variation"] is True


def test_actual_new_profile_quadratic_contact_and_force_changes():
    rows=profiles.checks()
    assert all(rows[name]==0 for name in rows if name.startswith("actual_replacement_"))
    d=profiles.profile_bounds(estimates.SCALE)
    assert d["changed_scale_square_coefficient_exactly_zero"] is True
    assert d["changed_mixed_lapse_scale_coefficient_upper"]==3*d["changed_lapse_square_coefficient_upper"]
    assert d["changed_lapse_square_coefficient_upper"]<sp.Rational(1,10**770)


def test_invalid_inputs_rejected_after_cache_population():
    estimates.physical_bounds(10**400)
    profiles.profile_bounds(10**400)
    readouts.rows("transverse",0)
    subtraction.reference_tail("longitudinal",0)
    assert audit.controls()["rejected_inputs"]==94


def test_exact_json_roundtrip_and_inexact_value_rejection():
    d=original.serialize({"bounds":estimates.physical_bounds(estimates.SCALE),
        "profile":profiles.action(),"integrands":profiles.integrands()})
    assert json.loads(json.dumps(d))==d
    for value in (1.0,sp.Float(1),sp.oo,sp.nan):
        with pytest.raises((TypeError,ValueError)):
            original.serialize(value)
