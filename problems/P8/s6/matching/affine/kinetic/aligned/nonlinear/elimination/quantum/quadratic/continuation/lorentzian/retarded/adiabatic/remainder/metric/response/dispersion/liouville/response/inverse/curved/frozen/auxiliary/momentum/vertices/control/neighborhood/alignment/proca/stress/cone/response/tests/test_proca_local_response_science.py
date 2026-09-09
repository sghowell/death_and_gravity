"""Actual ordinary Proca finite local response, chart contacts and exact bounds."""
import json

import pytest
import sympy as sp
from p8_affine import verify as original
from p8_proca_local_response import audit, bounds, chart, local


def test_all_new_local_response_exact_residuals():
    rows=audit.residuals()
    assert len(rows)==246
    assert sum(len(value) if isinstance(value,sp.MatrixBase) else 1 for value in rows.values())==250
    for value in rows.values():
        assert all(x==0 for x in value) if isinstance(value,sp.MatrixBase) else value==0


def test_all_local_response_gates():
    assert len(audit.gates())==15
    assert all(value is True for value in audit.gates().values())


@pytest.mark.parametrize("order",range(3))
def test_actual_mass_jets_are_zero_before_dimensional_matching(order):
    rows=local.checks()
    for output in ("N","Z"):
        assert rows[output+f"_actual_constant_mass_pole_{order}"]==0
        assert rows[output+f"_independent_ordinary_finite_heat_action_{order}"]==0


def test_no_old_mass_jet_or_unfixed_background_survives():
    fields=set(local.n+local.v+(local.u,))
    assert all(value==0 for value in local.ZERO.values())
    for order in range(3):
        for output in ("N","Z"):
            assert local.operator(output,order).free_symbols<=fields


@pytest.mark.parametrize("order",range(3))
def test_physical_readout_normalization_cancels_constant_spatial_rescaling(order):
    for output in ("energy","pressure"):
        assert local.physical_operator(output,order).coeff(local.v[0])==0


def test_constant_density_physical_stress_response_is_zero():
    assert local.physical_operator("energy",0)==0
    assert local.physical_operator("pressure",0)==0
    assert local.operator("N",0)!=0
    assert local.operator("Z",0)!=0


def test_actual_clock_map_keeps_coefficient_time_jets():
    mapping=chart.mapping()
    assert mapping["omega_N"]==1/(2*(1+local.u**2)**3)
    assert mapping["omega_NN"]==-sp.Rational(3,2)/chart.h+1/chart.h**2
    expected=chart.vhat[1]+mapping["omega_N"]*local.n[1]+sp.diff(mapping["omega_N"],local.u)*local.n[0]
    assert sp.factor(mapping["linear_physical_logscale_jet"][1]-expected)==0
    assert sp.diff(mapping["omega_N"],local.u)!=0


@pytest.mark.parametrize("order",range(3))
def test_independent_nonlinear_covariant_metric_map_has_same_Euler_currents(order):
    rows=chart.checks()
    assert rows[f"N_independent_nonlinear_metric_map_Euler_{order}"]==0
    assert rows[f"V_independent_nonlinear_metric_map_Euler_{order}"]==0


def test_metric_pullback_keeps_nonzero_second_map_contact():
    d=chart.density(0)
    linear_map=dict(zip(local.v,chart.mapping()["linear_physical_logscale_jet"]))
    linear=local.density(0).subs(linear_map,simultaneous=True)
    contact=local.clean(d["metric_linear_pullback_and_second_map_contact"]-linear,extra=(chart.vhat,))
    assert contact!=0
    assert sp.factor(contact-3*sp.Rational(5,2)*chart.mapping()["omega_NN"]*local.n[0]**2/2)==0


def test_local_profile_is_fixed_and_zero_order_cancellation_is_exact():
    d=chart.density(0)
    assert d["actual_fixed_local_profile_quadratic_density"]!=0
    assert d["metric_linear_pullback_and_second_map_contact"]!=0
    assert d["combined_local_quadratic_density"]==0
    assert chart.operator("N",0)==0
    assert chart.operator("V",0)==0


def test_fourth_local_matrix_has_rank_one_and_the_actual_null_direction():
    d=chart.top()
    actual=d["actual_fourth_derivative_coefficient"]
    assert actual==d["expected_rank_one_coefficient"]
    assert actual.det()==0
    assert actual[1,1]==-4
    assert sp.simplify(actual*d["kernel_vector"])==sp.zeros(2,1)
    assert sp.simplify(actual.subs(local.u,0))==sp.Matrix([[-1,-2],[-2,-4]])


def test_old_invertible_fourth_matrix_is_not_assumed():
    matrix=chart.top()["actual_fourth_derivative_coefficient"]
    with pytest.raises(sp.matrices.exceptions.NonInvertibleMatrixError):
        matrix.inv()


def test_weighted_hessian_self_adjointness_in_both_charts():
    rows=audit.self_adjoint_checks()
    assert len(rows)==90
    assert all(value==0 for value in rows.values())


def test_weighted_adjoint_is_not_flat_time_adjoint():
    assert audit.adjoint((sp.Integer(0),sp.Integer(1)))==(-3*local.H,-1)
    assert audit.adjoint((sp.Integer(1),))==(1,)


@pytest.mark.parametrize("group",("physical_finite_local_response","clock_chart_finite_local_plus_its_fixed_profile"))
def test_all_coefficient_enclosures_are_continuous_exact_reconstructions(group):
    for orders in bounds.envelopes()[group].values():
        for data in orders.values():
            assert len(data["box_reconstructions"])==10
            assert all(value==0 for value in data["box_reconstructions"].values())
            assert all(value>=0 for row in data["coefficient_absolute_envelopes"].values() for value in row.values())


@pytest.mark.parametrize("group",("physical_finite_local_response","clock_chart_finite_local_plus_its_fixed_profile"))
def test_declared_C4_to_C0_bounds_at_named_scale(group):
    data=bounds.at_scale(bounds.SCALE)
    assert data["M_tau"]==10**400
    assert data["fixed_mass_time_product"]==1000
    for row in data["response_bounds"][group].values():
        assert 0<row["C4_to_C0_upper"]<sp.Rational(1,10**790)
        assert row["adiabatic_half_order_bounds"][0]==0
    assert data["full_nonlocal_or_mixed_loop_bound"] is False


def test_exact_scale_law_and_no_hidden_mass_rescaling():
    a,b=bounds.at_scale(bounds.SCALE),bounds.at_scale(2*bounds.SCALE)
    assert a["fixed_mass_time_product"]==b["fixed_mass_time_product"]==1000
    for group,outputs in a["response_bounds"].items():
        for output,row in outputs.items():
            assert row["C4_to_C0_upper"]==4*b["response_bounds"][group][output]["C4_to_C0_upper"]


@pytest.mark.parametrize("bad",(local.n[0]**2,sp.Symbol("unknown")*local.n[0],local.n[5]))
def test_nonlinear_unknown_or_higher_derivative_sources_are_rejected(bad):
    with pytest.raises(ValueError):
        bounds.coefficient_bounds(bad,{"N":local.n,"Z":local.v})


def test_validation_precedes_populated_caches():
    local.density(1)
    chart.direct_metric_density(1)
    chart.operator("N",1)
    bounds.at_scale(bounds.SCALE)
    assert audit.controls()["rejected_inputs"]==100


def test_exact_serialization_and_float_rejection():
    value=original.serialize({"local_bounds":bounds.at_scale(bounds.SCALE),"top":chart.top()})
    assert json.loads(json.dumps(value))==value
    for bad in (1.0,sp.Float(1),sp.nan,sp.oo):
        with pytest.raises((TypeError,ValueError)):
            original.serialize(bad)

