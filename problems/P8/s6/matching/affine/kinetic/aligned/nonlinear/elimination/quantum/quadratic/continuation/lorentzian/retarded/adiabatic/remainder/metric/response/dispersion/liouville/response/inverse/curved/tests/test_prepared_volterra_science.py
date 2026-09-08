"""Scientific UV, finite-contact, primitive and literal-Euler controls."""
import pytest
import sympy as sp
from p8_prepared_volterra import primitive, proofs, tree, verify, vertices
from p8_vector_metric_dispersion import principal, spectral


def zero(values):
    for value in values.values():
        if isinstance(value,sp.MatrixBase):
            assert all(entry==0 for entry in value)
        else:
            assert value==0


def test_actual_physical_bilocal_and_dimensional_coefficients():
    zero(vertices.checks())


def test_all_momentum_dimensional_fourth_Taylor_contacts_without_freezing_curvature():
    values=vertices.adiabatic_principal_checks()
    assert len(values)==8
    zero(values)


def test_all_time_primitive_and_multiplier_commutator_identities():
    zero(primitive.checks())


def test_literal_matter_charge_and_both_classical_metric_Eulers():
    zero(tree.checks())


def test_complete_identity_entry_and_audit_counts():
    values=proofs.residuals()
    zero(values)
    assert len(values)==47
    assert sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in values.values())==72
    assert len(proofs.checks())==18
    assert all(value is True for value in proofs.checks().values())


def test_dimensional_pairs_keep_the_longitudinal_finite_jet():
    T,L=vertices.high("T")["pair"],vertices.high("L")["pair"]
    assert T.subs(vertices.D,3)==sp.zeros(2,1)
    assert sp.diff(T,vertices.D)==sp.Matrix([0,2])
    assert sp.diff((vertices.D-1)*T*T.T,vertices.D).subs(vertices.D,3)==sp.zeros(2)
    assert sp.diff(L,vertices.D)==sp.Matrix([0,2])


def test_scale_cancellation_has_no_dimensional_log_a_remainder():
    value=vertices.a**(-vertices.D-2)*(vertices.a*vertices.p)**(vertices.D+1)*vertices.a
    assert sp.simplify(value-vertices.p**(vertices.D+1))==0
    assert sp.simplify(sp.diff(value-vertices.p**(vertices.D+1),vertices.D))==0


def test_finite_fourth_Taylor_contact_survives_in_the_nonlogarithmic_channel():
    M=spectral.actual(spectral.matrix())
    w=principal.change()[:,0]
    contact=sp.factor((w.T*M*w)[0]/32)
    assert contact.subs({spectral.z:sp.Rational(1,2),spectral.h:1})==sp.Rational(1,17496)
    assert (w.T*principal.pole()*w)[0].simplify()==0


def test_potentials_do_not_enter_the_two_highest_bilocal_powers():
    poly=vertices.polynomial()
    for j in (3,4):
        assert not poly.nth(j).has(vertices.Ut,vertices.Us)
    assert poly.nth(2).has(vertices.Ut,vertices.Us)


def test_connection_subleading_fixture_is_nonzero_and_antisymmetric():
    A=vertices.high("L")["A"].subs({vertices.D:3,vertices.h:1})
    b=vertices.high("L")["pair"].subs({vertices.D:3,vertices.h:1})
    value=2*(b*A.T-A*b.T)
    assert value==sp.Matrix([[0,sp.Rational(248,27)],[-sp.Rational(248,27),0]])


def test_curved_lag_geometry_is_not_frozen_away_from_diagonal():
    tau=sp.Symbol("lag",real=True)
    fixture=(1+tau)**(-4)*(1-tau/2)**(-5)
    assert fixture.subs(tau,0)==1
    assert sp.diff(fixture,tau).subs(tau,0)==-sp.Rational(3,2)


def test_local_coefficient_commutator_is_not_optional():
    c=sp.Function("local_coefficient")(primitive.s)
    actual=primitive.local_kernel(2).subs(c,primitive.s).doit()
    naive=(primitive.t-primitive.s)*primitive.s
    assert sp.factor(actual-naive)==-(primitive.t-primitive.s)**2


def test_constant_fourth_local_coefficient_requires_the_instantaneous_term():
    c=sp.Function("local_coefficient")(primitive.s)
    regular=primitive.local_kernel(4).subs(c,1).doit()
    assert regular==0
    assert sp.diff(primitive.lag**3/6,primitive.lag,3)==1


def test_subleading_primitive_is_logarithmic_but_leading_one_is_inverse_lag():
    assert primitive.primitive_shape(3)==-sp.log(primitive.lag)/6
    assert primitive.primitive_shape(4)==1/(24*primitive.lag)


def test_short_interval_log_majorant_has_integral_tending_to_zero():
    T=sp.Symbol("positive_T",positive=True)
    bound=T*(2-sp.log(T))
    assert sp.limit(bound,T,0,dir="+")==0
    assert bound.subs(T,sp.Rational(1,2))>0


def test_wrong_matter_reconstruction_sign_violates_the_conserved_charge():
    data=tree.data()
    wrong=data["zero_prepared_charge"].subs(sp.diff(tree.matter,tree.u),3*tree.ell*tree.v+tree.w*tree.n)
    assert sp.simplify(wrong.subs({tree.ell:1,tree.v:1,tree.n:0}))==6
    assert sp.simplify(data["zero_prepared_charge"].subs(sp.diff(tree.matter,tree.u),data["matter_reconstruction_rate"]))==0


def test_physical_tree_highest_matrix_is_regular_and_not_inverted():
    currents=tree.data()["physical_metric_Euler_currents"]
    matrix=sp.Matrix([[sp.diff(row,sp.diff(field,tree.u,2)) for field in (tree.n,tree.zeta)] for row in currents])
    assert matrix.subs(tree.delta,sp.Rational(1,2))==sp.Matrix([[sp.Rational(3,2),-3],[-3,6]])
    assert sp.factor(matrix.det())==0


def test_strict_native_interfaces_remain_strict_after_cache_population():
    vertices.high("L")
    primitive.local_kernel(4)
    primitive.primitive_shape(2)
    assert verify.controls()["rejected_inputs"]==40
    for value in (True,sp.Integer(2),2.0,-1,5,"2",None):
        with pytest.raises((TypeError,ValueError)):
            primitive.local_kernel(value)


def test_independent_prepared_polynomial_IBP_integral():
    t,s=primitive.t,primitive.s
    c=sp.Function("local_coefficient")(s)
    original=sp.integrate((t-s)**3*(1+s)*12*s**2/6,(s,0,t))
    transformed=sp.integrate(primitive.local_kernel(2).subs(c,1+s).doit()*s**4,(s,0,t))
    assert sp.simplify(original-transformed)==0
