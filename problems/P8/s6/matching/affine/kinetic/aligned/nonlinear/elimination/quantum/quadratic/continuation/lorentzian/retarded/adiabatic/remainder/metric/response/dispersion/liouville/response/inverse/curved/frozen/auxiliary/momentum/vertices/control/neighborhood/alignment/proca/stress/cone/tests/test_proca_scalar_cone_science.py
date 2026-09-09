"""Independent central time jets, physical scalar cones and finite-domain limits."""
import json

import pytest
import sympy as sp
from p8_affine import verify as original
from p8_proca_scalar_cone import audit, background, domain, principal, spatial


def test_all_exact_central_scalar_residuals():
    rows=audit.residuals()
    assert len(rows)==46
    assert sum(len(v) if isinstance(v,sp.MatrixBase) else 1 for v in rows.values())==91
    for value in rows.values():
        assert all(x==0 for x in value) if isinstance(value,sp.MatrixBase) else value==0


def test_all_domain_and_scope_gates():
    assert len(audit.gates())==32
    assert all(v is True for v in audit.gates().values())


def test_actual_homogeneous_canonical_equations():
    rows=background.checks()
    assert rows["actual_homogeneous_scale_Hamilton_equation"]==0
    assert rows["actual_homogeneous_trace_Hamilton_equation"]==0


def test_actual_trace_time_derivative_keeps_boundary_primitive():
    d=background.data()
    assert d["actual_trace_linear_time_derivative"]!=0
    assert d["actual_boundary_primitive_phi"]!=0
    assert background.checks()["actual_bounce_trace_momentum_derivative"]==0


def test_actual_hat_Hubble_is_not_canonical_background_coefficient():
    d=background.data()
    difference=sp.factor(d["hat_Hubble_time_derivative"]+d["trace_momentum_time_derivative"]/2)
    assert difference!=0
    assert difference.subs(background.N,1)==0


def test_actual_time_parity_and_zero_central_lapse_velocity():
    rows=background.checks()
    assert all(value==0 for name,value in rows.items() if "even_" in name or "odd_" in name)
    assert background.data()["lapse_time_derivative"]==0
    assert background.data()["matter_density_time_derivative"]==0


def test_partial_lapse_derivative_not_differentiation_along_constraint_family():
    d=background.data()
    N=background.N
    wrong=sp.diff(N*(d["trace_momentum_time_derivative"]
                      -d["actual_trace_linear_time_derivative"])/(2*d["trace_a"]),N)
    assert sp.factor(wrong-d["trace_force_time_derivative"])!=0
    assert principal.checks()["wrong_constraint_family_derivative_negative_control"]==0


@pytest.mark.parametrize("left",spatial.CHANNELS)
def test_actual_spatial_phase_pairs_and_constraints(left):
    for right in spatial.CHANNELS:
        result=spatial.pair(left,right)
        assert all(result["spatial_checks"].values())
        assert result["fixed_phase_derivative_precedes_background_equations"] is True
        reverse=spatial.pair(right,left)["kernel_first_time_jet"].subs(spatial.k,-spatial.k)
        assert spatial.clean(result["kernel_first_time_jet"]-reverse)==0


def test_moving_trace_boundary_is_kept_in_previous_matrix_comparison():
    rows=principal.checks()
    assert rows["previous_central_spatial_scalar_matrix_plus_actual_boundary"]==sp.zeros(4)
    actual=spatial.pair("curvature","curvature")["moving_boundary_kernel"]
    assert actual.subs(spatial.eta,0)==3*spatial.at_lapse(background.data()["trace_momentum_time_derivative"])


def test_full_canonical_gamma_map():
    assert principal.checks()["actual_scalar_gamma_map_symplectic"]==sp.zeros(4)
    assert principal.blocks()["gamma_canonical_map"][1,0]==-2*principal.q


def test_kinetic_limit_is_derived_from_actual_phase_matrix():
    assert principal.checks()["actual_high_q_kinetic_from_full_spatial_time_jet"]==sp.zeros(2)
    assert principal.blocks()["high_q_mixed_over_q"]==sp.zeros(2)


def test_clock_gradient_requires_actual_mixed_time_derivative():
    b=principal.blocks()
    assert b["high_q_potential_over_q"][0,0]==0
    assert b["high_q_mixed_time_derivative_over_q"][0,0].subs(spatial.t,1)==12
    assert b["Euler_first_gradient"][0,0].subs(spatial.t,1)==12


def test_freeze_before_Euler_gives_false_negative_gradient_determinant():
    b=principal.blocks()
    point={spatial.t:1,spatial.P:sp.Rational(1,10)}
    wrong=b["high_q_potential_over_q"].subs(point)
    right=b["Euler_first_gradient"].subs(point)
    assert wrong.det()==-sp.Rational(1,100)
    assert right.det()==sp.Rational(1199,100)


def test_characteristic_factorization_uses_the_physical_matter_metric():
    f=principal.formula()
    assert f["actual_lapse"]==background.N
    assert f["physical_spatial_metric_over_hat_metric"]==background.N
    assert f["physical_squared_momentum"]==principal.q/background.N
    assert principal.checks()["exact_physical_characteristic_factorization"]==0


def test_exact_luminal_matter_direction_and_global_exceptional_identity():
    f=principal.formula()
    assert f["matter_physical_speed_squared"]==1
    assert f["physical_cone_comparison_matrix"][0,1]==0
    assert f["physical_cone_comparison_matrix"][1,1]==0
    assert background.checks()["global_exceptional_luminal_matter_relation"]==0


def test_original_clock_margin_is_recovered():
    assert principal.formula()["clock_physical_speed_squared"].subs(background.N,1)==sp.Rational(749375,749377)


@pytest.mark.parametrize("radius",(sp.Rational(1,10**28),sp.Rational(1,10**10),domain.RADIUS))
def test_continuous_positive_subluminal_domain(radius):
    d=domain.domain(radius)
    assert all(d["gates"].values())
    boxes=d["continuous_rational_enclosures"]
    assert boxes["physical_clock_cone_margin"]["lower"]>sp.Rational(1,10**6)
    assert boxes["clock_gradient_Schur_over_sqrt_N"]["lower"]>11
    assert d["not_the_old_nine_invariant_interaction_polydisc"] is True


@pytest.mark.parametrize("lapse",(1-domain.RADIUS,1,1+domain.RADIUS))
def test_actual_central_constraint_data_are_positive_and_subluminal(lapse):
    p=domain.point(lapse)
    assert p["positive_matter_momentum_squared"]>0
    assert p["fixed_phase_lapse_Hessian"]<0
    assert p["clock_gradient_Schur"]>0
    assert 0<p["clock_physical_speed_squared"]<1


def test_larger_interval_is_not_falsely_declared_subluminal():
    d=domain.domain(domain.OUTER_RADIUS)
    assert d["gates"]["strictly_subluminal_clock_on_this_whole_interval"] is False
    assert d["gates"]["positive_matter_root"] is True
    assert d["gates"]["positive_scalar_gradient_Schur"] is True


def test_unique_simple_cone_crossing_bracket():
    b=domain.boundary()
    assert b["root_bracket"]==(1+sp.Rational(2,10**7),1+sp.Rational(21,10**8))
    assert b["polynomial_at_lower"]<0<b["polynomial_at_upper"]
    assert b["positive_polynomial_derivative_enclosure"]["lower"]>0
    assert b["denominator_on_outer_interval"]["upper"]<0


def test_actual_superluminal_constraint_datum_is_not_a_ghost_or_gradient_failure():
    b=domain.boundary()
    p=b["actual_superluminal_fixed_phase_constraint_datum"]
    assert p["lapse"]==1+sp.Rational(1,10**6)
    assert b["control_is_ghost_and_gradient_positive"] is True
    assert b["control_is_strictly_superluminal"] is True
    assert p["clock_physical_speed_squared"]==sp.Rational(
        5995038990009995014500031625011,5994976989968995020500043125015)


def test_continuous_enclosure_detects_an_uncontrolled_pole():
    with pytest.raises(ValueError,match="denominator"):
        domain.enclosure(1/(background.N-1),domain.RADIUS)


def test_inexact_and_outside_domain_controls_after_cache_population():
    domain.domain(domain.RADIUS)
    domain.point(1)
    spatial.pair("curvature","scalar_p")
    assert audit.controls()["rejected_inputs"]==68


def test_exact_json_roundtrip_and_float_rejection():
    d=original.serialize({"background":background.data(),"principal":principal.formula(),
        "domain":domain.domain(),"boundary":domain.boundary()})
    assert json.loads(json.dumps(d))==d
    for value in (1.0,sp.Float(1),sp.oo,sp.nan):
        with pytest.raises((TypeError,ValueError)):
            original.serialize(value)

