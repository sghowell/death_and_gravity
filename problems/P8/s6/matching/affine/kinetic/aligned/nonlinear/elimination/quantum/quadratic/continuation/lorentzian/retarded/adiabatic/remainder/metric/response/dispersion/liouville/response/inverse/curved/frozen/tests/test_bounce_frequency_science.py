"""Independent exact checks for the frozen diagnostic and its strict boundary."""
import pytest
import sympy as sp
from p8_bounce_frequency import bracket, verify
from p8_bounce_frequency import coefficients as coeff
from p8_prepared_volterra import tree
from p8_vector_metric_response import tadpole


def test_all_exact_residuals_and_continuous_gates():
    verify.affine.certify_residuals(verify.residuals())
    assert all(value is True for value in bracket.gates().values())


def test_actual_bounce_jets_are_not_flattened():
    point=coeff.background_jets()
    assert point[tree.theta]==0
    assert point[sp.diff(tree.theta,tree.u)]==3
    assert point[sp.diff(tree.delta,tree.u,2)]==-3
    assert point[tree.Je]==sp.Rational(1199,800)+4*sp.Rational(1,10**6)


def test_tree_regular_coefficients_from_physical_euler_operator():
    E0,E1,E2=coeff.tree_rows()
    assert E0==sp.Matrix([[sp.Rational(749377,250000),-sp.Rational(3,200)],
                          [-sp.Rational(3,200),-sp.Rational(9,100)]])
    assert E1==sp.zeros(2)
    assert E2==sp.diag(0,6)


def test_dropping_theta_prime_is_a_nonzero_freezing_defect():
    physical=tree.data()["physical_metric_Euler_currents"]
    C=sp.Matrix([[sp.diff(row,field) for field in (tree.n,tree.zeta)] for row in physical])
    wrong=coeff.background_jets().copy()
    wrong[sp.diff(tree.theta,tree.u)]=0
    defect=coeff.T.T*C.xreplace(wrong)*coeff.T-coeff.tree_rows()[0]
    assert defect.applyfunc(sp.factor)==sp.Matrix([[0,0],[18,0]])


def test_all_local_orders_and_nonsymmetric_lower_rows_retained():
    assert coeff.local_row(1)==coeff.local_row(3)==sp.zeros(2)
    assert coeff.local_row(4)==coeff.finite()
    assert sp.factor(coeff.local_row(2)[0,1]-coeff.local_row(2)[1,0])==-sp.Rational(1712,27)


def test_finite_lower_local_absolute_maxima():
    data=bracket.data()
    assert data["local_zero_entrywise_absolute_max"]==22500239999424
    assert data["local_second_entrywise_absolute_max"]==19999856


def test_fixed_tadpole_uses_euler_hessian_not_stress_jacobian():
    data=tadpole.physical_vertices()
    actual=coeff.T.T*data["density_hessian"]*coeff.T
    assert (actual-coeff.tadpole_matrix()).applyfunc(sp.expand)==sp.zeros(2)
    wrong=coeff.T.T*data["physical_stress_jacobian"]*coeff.T
    assert (wrong-actual).applyfunc(sp.expand)!=sp.zeros(2)


def test_fixed_profile_box_is_physical_not_unscaled():
    box=coeff.profile_box()
    original=coeff.pinned(coeff.response.REPORT,coeff.RESPONSE_SHA)
    stored=original["complete_background_cancelled_response_scale_example"]["existing_tadpole_physical_C0_to_C0_upper_bound"]
    assert sum(box.values())==sp.Rational(stored)
    assert 5760*coeff.L**2*sum(box.values())==sp.Rational(84242272974349549,1296)


def test_positive_pair_gram_and_integer_envelopes():
    data=bracket.data()
    assert data["pair_Gram_entrywise_raw_upper"]==sp.Matrix([
        [sp.Rational(24946,2187),sp.Rational(1340,81)],
        [sp.Rational(1340,81),24]])
    assert all(x>=0 for x in bracket.GRAM_UPPER-data["pair_Gram_entrywise_raw_upper"])
    assert coeff.gram().subs(coeff.spectral.z,sp.Rational(1,2)).is_positive_definite


def test_exact_log_bound_does_not_use_leading_log_approximation():
    data=bracket.data()
    assert data["log_argument_upper"]==25000000000000000001
    assert sp.Integer(100)**18/sp.factorial(18)>data["log_argument_upper"]
    assert bracket.BLOCK_UPPER==sp.Matrix([[302,428],[428,604]])


def test_pi_bounds_are_exact_positive_integral_bounds():
    data=bracket.data()
    assert data["pi_lower_partial_sum"]==sp.Rational(135904,45045)>3
    assert data["pi_upper"]**2<10
    x=sp.Symbol("x",real=True)
    assert sp.integrate(x**4*(1-x)**4/(1+x*x),(x,0,1))==sp.Rational(22,7)-sp.pi


def test_lower_local_and_fixed_tadpole_uniform_error():
    data=bracket.data()
    assert data["normalized_tadpole_entrywise_upper"]<10**15
    assert data["full_lower_order_entrywise_error_upper"]==sp.Rational(
        16666666666666666667,96000000000000000000000000000000000)
    assert data["full_lower_order_entrywise_error_upper"]<sp.Rational(1,10**12)


def test_lower_endpoint_determinant_strictly_positive():
    data=bracket.endpoint("lower")
    assert data["NN_lower"]>2
    assert data["ZZ_lower"]>5*10**24
    assert data["both_off_diagonal_absolute_upper"]<1
    assert data["determinant_comparison"]>0


def test_upper_endpoint_determinant_strictly_negative():
    data=bracket.endpoint("upper")
    assert data["NN_upper"]<-12
    assert data["ZZ_lower"]>5*10**26
    assert data["both_off_diagonal_absolute_upper"]<10000
    assert data["determinant_comparison"]<0


def test_entire_bracket_ZZ_cofactor_stays_positive():
    assert bracket.data()["ZZ_positive_lower_on_entire_band"]>0
    A,B,C,D=sp.symbols("A B C D")
    assert sp.Matrix([[A,B],[C,D]]).adjugate()[0,0]==D


def test_tree_alone_has_no_zero_in_this_band():
    E0,_,E2=coeff.tree_rows()
    s=sp.Symbol("s",positive=True)
    determinant=sp.expand((E0+s*s*E2).det())
    assert determinant.subs(s,bracket.LOW)>0
    assert sp.diff(determinant,s).is_positive


def test_validation_runs_before_cache_and_rejects_conversions():
    coeff.local_row(1)
    bracket.endpoint("lower")
    assert verify.controls()["rejected_inputs"]==28
    for value in (True,sp.Integer(1),1.0):
        with pytest.raises(ValueError):
            coeff.local_row(value)


def test_pinned_ancestor_report_cannot_silently_drift(monkeypatch):
    monkeypatch.setattr(coeff.affine,"sha",lambda path:"wrong source hash")
    with pytest.raises(ValueError,match="changed"):
        coeff.pinned(coeff.local.REPORT,coeff.LOCAL_SHA)


def test_stationary_scope_does_not_claim_actual_instability_or_cutoff():
    description=bracket.description()
    assert "actual nonstationary pole" in description["not_claimed"]
    assert "justified interacting cutoff" in description["not_claimed"]
    assert "every" in description["uniformity"].lower()
    assert bracket.LOW/coeff.L==sp.Rational(1,10**12)
    assert bracket.HIGH/coeff.L==sp.Rational(1,10**11)
