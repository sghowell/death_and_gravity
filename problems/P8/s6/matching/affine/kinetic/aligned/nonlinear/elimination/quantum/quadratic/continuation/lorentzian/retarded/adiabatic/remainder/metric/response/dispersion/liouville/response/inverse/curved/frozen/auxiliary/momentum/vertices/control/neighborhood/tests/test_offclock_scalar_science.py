"""Scientific checks without the recursive parent certificate build."""
import json

import pytest
import sympy as sp
from p8_affine import verify as affine
from p8_auxiliary_neighborhood import model
from p8_offclock_scalar import audit, datum, euler, principal, spatial


def test_all_exact_coupled_offclock_residuals():
    rows=audit.residuals()
    assert len(rows)==45
    assert sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in rows.values())==97
    for value in rows.values():
        assert all(x==0 for x in value) if isinstance(value,sp.MatrixBase) else value==0


def test_all_actual_domain_and_scope_gates():
    assert len(audit.gates())==21
    assert all(value is True for value in audit.gates().values())


def test_constraint_family_from_fixed_phase_lapse_derivative():
    d=datum.data()
    H=d["homogeneous_Hamiltonian"]
    P2=d["matter_momentum_squared"]
    assert sp.factor(sp.diff(H,model.N).subs(datum.P**2,P2))==0
    assert sp.factor(sp.diff(H,model.N,2).subs(datum.P**2,P2)-d["fixed_phase_lapse_Hessian"])==0


def test_constraint_tangent_is_not_the_auxiliary_pivot():
    d=datum.data()
    along=sp.diff(sp.diff(d["homogeneous_Hamiltonian"],model.N).subs(
        datum.P**2,d["matter_momentum_squared"]),model.N)
    assert sp.factor(along)==0
    assert d["fixed_phase_lapse_Hessian"].subs(model.N,1)<0


@pytest.mark.parametrize("point",(1-datum.RADIUS,sp.Integer(1),1+datum.RADIUS))
def test_exact_family_points_stay_in_continuous_bounds(point):
    d=datum.data()
    P2=d["matter_momentum_squared"].subs(model.N,point)
    assert sp.Rational(1,200)<P2
    assert (sp.Rational(1,10)-sp.Rational(1,10**24))**2<P2
    assert P2<(sp.Rational(1,10)+sp.Rational(1,10**24))**2
    assert d["fixed_phase_lapse_Hessian"].subs(model.N,point)<0
    assert d["joint_temporal_auxiliary_pivot"].subs(model.N,point)<0


def test_continuous_lapse_polynomial_remainder_proofs():
    proofs=datum.domain()["continuous_rational_proofs"]
    assert len(proofs)==6
    for row in proofs.values():
        assert row["denominator_center"]-row["denominator_tail_upper"]>0
        assert row["lower"]<=row["upper"]
        assert row["radius"]==sp.Rational(1,10**28)


def test_actual_temporal_pivot_includes_trace_legendre_coupling():
    d=datum.data()
    assert datum.checks()["full_joint_temporal_pivot_not_isolated_mass"]==0
    assert datum.checks()["joint_mixed_lapse_temporal_pivot_zero_at_actual_datum"]==0
    correction=sp.factor(d["gamma_t"]-d["longitudinal_kinetic_times_sqrt_N"])
    assert sp.factor(correction-sp.Rational(3,2)*model.N**2*d["source_delta"]**2)==0
    assert correction.subs(model.N,1-datum.RADIUS)>0


def test_original_clock_datum_is_recovered_exactly():
    d=datum.data()
    assert d["matter_momentum_squared"].subs(model.N,1)==sp.Rational(1,100)
    assert d["fixed_phase_lapse_Hessian"].subs(model.N,1)==-sp.Rational(749377,250000)
    assert d["longitudinal_kinetic_times_sqrt_N"].subs(model.N,1)==1
    assert d["canonical_gamma_vector_mixing"].subs(model.N,1)==0


def test_no_independent_lapse_choice_or_physical_metric_redefinition():
    d=datum.data()
    assert d["physical_spatial_metric_over_hat_metric"]==model.N
    assert d["physical_lapse"]==model.N
    assert datum.domain()["positive_matter_momentum_deviation_upper"]<sp.Rational(1,10**24)


def test_all_spatial_constraints_for_every_phase_pair():
    count=0
    for a,left in enumerate(spatial.CHANNELS):
        for right in spatial.CHANNELS[a:]:
            row=spatial.pair(left,right)
            assert all(row["spatial_checks"].values())
            assert row["matter_density_relation_imposed_after_N_derivative"] is True
            count+=len(row["spatial_checks"])
    assert count==63


def test_independent_full_Fourier_reversal_and_label_exchange():
    for a,left in enumerate(spatial.CHANNELS):
        for right in spatial.CHANNELS[a+1:]:
            assert spatial.root_relation(spatial.pair(right,left)["kernel"]-
                spatial.pair(left,right)["kernel"].subs(spatial.k,-spatial.k))==0


def test_York_cancellation_does_not_delete_scalar_vector_coupling():
    assert spatial.pair("scalar_p","scalar_p")["kernel"]==0
    coupling=spatial.pair("scalar_p","longitudinal_p")["kernel"]
    assert coupling!=0
    assert sp.factor(coupling+spatial.at_lapse(datum.data()["canonical_gamma_vector_mixing"])/2)==0
    assert coupling.subs(spatial.t,1)==0


def test_lapse_force_does_not_cancel_longitudinal_cross_term():
    for pair in (("scalar_p","scalar_p"),("longitudinal_p","longitudinal_p")):
        row=spatial.pair(*pair)
        assert row["lapse_force_left"]==row["lapse_force_right"]==0
    row=spatial.pair("scalar_p","longitudinal_p")
    assert row["kernel"]==row["unreduced_kernel"]


def test_gamma_map_is_exactly_canonical_without_Theta_division():
    check=principal.checks()["full_spatial_gamma_canonical_map"]
    assert check==sp.zeros(6)
    T=principal.blocks()["canonical_gamma_map"]
    assert T[0,3]==1/(2*principal.q)
    assert T[1,0]==-2*principal.q


def test_actual_mixing_is_only_the_expected_q_column_at_slice():
    B=principal.blocks()["momentum_coordinate_block"]
    expected=sp.zeros(3)
    expected[2,0]=principal.q*spatial.at_lapse(datum.data()["canonical_gamma_vector_mixing"])
    assert (B-expected).applyfunc(principal.clean)==sp.zeros(3)


@pytest.mark.parametrize("point",(1-datum.RADIUS,1+datum.RADIUS))
def test_full_finite_q_kinetic_matrix_positive_at_exact_anchors(point):
    A=principal.blocks()["momentum_Hessian"]
    for size in (1,2,3):
        minor=principal.clean(A[:size,:size].det())
        value=minor.subs({spatial.t:sp.sqrt(point),principal.q:10**80})
        assert sp.factor(value)>0


def test_positive_limit_is_coupled_and_not_independent_oscillators():
    A=principal.blocks()["momentum_Hessian_limit"]
    assert A[0,1]!=0
    assert A[0,2]==A[1,2]==0
    assert principal.checks()["leading_light_kinetic_determinant"]==0
    assert principal.checks()["actual_three_scalar_high_q_kinetic_Hessian"]==sp.zeros(3)


def test_full_longitudinal_mass_and_electric_weights():
    assert principal.checks()["actual_longitudinal_gradient_coefficient"]==0
    assert principal.checks()["actual_longitudinal_finite_q_kinetic_correction"]==0
    assert principal.checks()["actual_matter_gradient_leading_coefficient"]==0


def test_bare_coordinate_square_is_absent_but_Legendre_square_is_not():
    C=principal.blocks()["coordinate_Hessian"]
    gamma=euler.actual_orders()["Lagrangian_coordinate_matrix"]
    assert C[0,0]==0
    assert gamma[0,0]!=0
    assert euler.checks()["actual_Lagrangian_q_squared_negative_entry"]==0


def test_arbitrary_allowed_time_jets_cannot_cancel_leading_determinant():
    d=euler.generic_determinant()
    assert d["residual"]==0
    assert d["leading_determinant_over_q_fourth"]==d["expected_leading_coefficient"]
    assert len(d["matrix_with_all_allowed_time_jet_orders"].free_symbols)>12


@pytest.mark.parametrize("point",(1-datum.RADIUS,1+datum.RADIUS))
def test_actual_Euler_leading_coefficient_negative_off_clock(point):
    coefficient=euler.result()["Euler_constant_determinant_leading_coefficient"]
    assert coefficient.subs(model.N,point)<0
    assert coefficient.subs(model.N,1)==0


def test_false_frozen_Hamiltonian_clock_sign_is_detected():
    assert euler.result()["clock_negative_control_without_Euler_time_jets"]==-sp.Rational(1,100)
    assert euler.result()["actual_clock_Euler_leading_determinant"]==sp.Rational(1199,100)
    assert euler.checks()["actual_Euler_first_clock_sign_is_positive_not_false_frozen_H_sign"]==0


def test_fixed_metric_shortcut_would_remove_the_actual_coefficient():
    coefficient=euler.result()["Euler_constant_determinant_leading_coefficient"]
    assert coefficient!=0
    assert datum.data()["canonical_gamma_vector_mixing"]!=0
    assert audit.gates()["no_numerical_cutoff_or_UV_exclusion_from_q_infinity"] is True


def test_domain_validation_precedes_pair_cache():
    spatial.pair("curvature","curvature")
    assert audit.controls()["rejected_inputs"]==36


def test_zero_denominator_and_inexact_lapse_boxes_rejected():
    for value in (1/(model.N-1),1/(model.N-1+datum.RADIUS/2),1.0,sp.Float(1),True,sp.sqrt(2)):
        with pytest.raises((TypeError,ValueError)):
            datum.enclosure(value)


def test_exact_matrices_domain_and_diagnostic_serialization_roundtrip():
    data=affine.serialize({"datum":datum.data(),"domain":datum.domain(),"blocks":principal.blocks(),
                           "orders":euler.actual_orders(),"result":euler.result()})
    assert json.loads(json.dumps(data))==data
    for value in (1.0,sp.Float(1),sp.oo,sp.nan):
        with pytest.raises((TypeError,ValueError)):
            affine.serialize(value)
