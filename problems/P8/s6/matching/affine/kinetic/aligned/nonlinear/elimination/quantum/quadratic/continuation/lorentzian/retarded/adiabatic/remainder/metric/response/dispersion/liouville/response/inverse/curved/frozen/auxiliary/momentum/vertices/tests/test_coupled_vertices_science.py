"""Scientific replay without the recursive parent certificate build."""
import json
from fractions import Fraction

import pytest
import sympy as sp
from p8_auxiliary_neighborhood import model
from p8_coupled_momentum import constraints
from p8_coupled_vertices import (
    anchors,
    audit,
    coefficients,
    majorant,
    physical,
    taylor,
    verify,
)


def test_all_physical_and_symbolic_residuals():
    rows=audit.residuals()
    assert len(rows)==154
    assert sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in rows.values())==455
    for value in rows.values():
        assert all(x==0 for x in value) if isinstance(value,sp.MatrixBase) else value==0


def test_all_actual_audit_gates():
    assert len(audit.gates())==17
    assert all(value is True for value in audit.gates().values())


def test_complete_lapse_coefficient_inventory():
    rows=coefficients.jets()
    assert len(rows)==14
    assert sum(map(len,rows.values()))==70
    assert sum(sum(p)==1 for p in rows)==9
    assert sum(sum(p)==2 for p in rows)==4


def test_original_I_phi_fourth_derivative_is_not_discarded():
    d=coefficients.boundary_jets()
    row=taylor.derivatives(d["Iphi"])
    assert row[0]==0
    assert row[4].subs(model.u,0)!=0
    assert row[4]==sp.factor(sp.diff(taylor.derivatives(d["I"])[4],model.u))


def test_original_Q_ODE_and_fixed_basepoint():
    checks=coefficients.checks()
    assert len(checks)==90
    assert checks["original_Q_fixed_basepoint"]==0
    assert checks["original_Q_second_lapse_derivative"]==0
    assert all(checks[f"original_Q_ODE_lapse_jet_{n}"]==0 for n in range(4))


def test_all_time_quadratic_algebraic_bridge():
    checks=physical.quadratic_algebraic_bridge()
    assert len(checks)==20
    assert all(value==0 for value in checks.values())


@pytest.mark.parametrize("point",(sp.Rational(-1,2),0,sp.Rational(1,2)))
def test_reduced_quadratic_scalar_vector_tensor_physical_samples(point):
    for pair in (("curvature","curvature"),("scalar_p","matter"),("Wy","Wy"),("Px","Px"),("tensor2","tensor2")):
        assert physical.quadratic_pair(*pair,point=point)==0


@pytest.mark.parametrize("n,point",audit.POINTS)
def test_complete_mixed_spatial_and_lapse_constraints(n,point):
    data=physical.fixture(n,point)
    assert all(data["momentum"]["checks"].values())
    assert len(data["momentum"]["checks"])==3*(n-1)
    assert len(data["reduction"]["checks"])==n+1
    assert all(data["reduction"]["checks"].values())
    assert data["full_labelled_kernel"]!=0


def test_independent_transverse_electric_quartic_and_factorial():
    data=anchors.pure_electric()
    Je=sp.Rational(1199,800)+4*model.MARGIN
    expected=sp.factor(sp.factorial(4)/(64*Je*sp.Rational(1,10**6)**2))
    assert data["kernel"]==expected==sp.Rational(187500000000000000,749377)
    assert data["residual"]==0


def test_omitted_lapse_response_is_detected():
    assert anchors.pure_electric()["omitting_lapse_response_quartic_defect"]<0


def test_Fourier_reality_is_full_momentum_reversal():
    assert sp.im(physical.fixture(3)["full_labelled_kernel"])!=0
    for n in (3,4):
        data=anchors.transformed_fixture(n,reverse=True)
        assert data["actual"]==sp.conjugate(physical.fixture(n)["full_labelled_kernel"])
        assert data["residual"]==0


def test_all_label_permutation_keeps_fields_with_their_wavevectors():
    for n in (3,4):
        data=anchors.transformed_fixture(n,permutation=tuple(reversed(range(n))))
        assert data["actual"]==physical.fixture(n)["full_labelled_kernel"]
        assert data["residual"]==0


def test_continuous_denominator_and_independent_Cauchy_certificates():
    bound=majorant.coefficient_bounds()
    for powers,values in coefficients.jets().items():
        for n,value in enumerate(values):
            _num,den=sp.fraction(sp.cancel(value))
            poly=sp.Poly(den,model.u)
            assert poly.nth(0)>0
            assert all(c>=0 and not d[0]%2 for d,c in poly.terms())
            proof=bound["rows"][powers]["proof_data"][n]
            assert proof["denominator_positive_even_constant_lower"]==poly.nth(0)
            assert proof["independent_S6_75_Cauchy_bound"]==10**4*sp.factorial(n)*100**n


def test_actual_coefficients_below_declared_bounds_at_endpoints_and_center():
    bound=majorant.coefficient_bounds()
    for powers,values in coefficients.jets().items():
        for n,value in enumerate(values):
            for point in (sp.Rational(-1,2),0,sp.Rational(1,2)):
                assert abs(value.subs(model.u,point))<=bound["rows"][powers]["derivative_bounds"][n]


def test_every_mixed_fixture_below_continuous_raw_majorant():
    bounds=majorant.bound(1,4,1)["physical_Hamiltonian_degree_majorant"].coefficients
    for n,point in audit.POINTS:
        value=physical.fixture(n,point)["full_labelled_kernel"]
        assert abs(sp.re(value))+abs(sp.im(value))<bounds[n]


def test_raw_field_degree_scaling_is_exact():
    first=majorant.bound(1,4,1)["physical_Hamiltonian_degree_majorant"].coefficients
    second=majorant.bound(2,4,1)["physical_Hamiltonian_degree_majorant"].coefficients
    assert all(second[n]==2**n*first[n] for n in range(5))


def test_full_vector_source_changes_higher_York_bounds():
    full=majorant.bound(1,4,1)["York_correction_majorants"]
    omitted=majorant.bound(1,4,1,include_vector=False)["York_correction_majorants"]
    assert full[0]==omitted[0]
    assert full[1]>omitted[1]
    assert full[2]>omitted[2]


def test_hard_wave_example_is_raw_and_not_a_cutoff():
    data=majorant.bound(1,10**9,Fraction(1,10**8))
    assert data["cubic_labelled_kernel_upper"]>0
    assert data["quartic_labelled_kernel_upper"]>0
    assert "no M*tau factors or free mode functions" in data["normalization"]


def test_positive_majorant_is_monotone_in_derivative_and_inverse_transfer():
    a=majorant.bound(1,4,1)["physical_Hamiltonian_degree_majorant"].coefficients
    b=majorant.bound(1,8,2)["physical_Hamiltonian_degree_majorant"].coefficients
    assert all(x<=y for x,y in zip(a,b))
    assert a[3]<b[3] and a[4]<b[4]


def test_taylor_fixed_basepoint_and_native_domain():
    ctx,e=taylor.context(model.u)
    value=taylor.power(1+e,sp.Rational(1,2))
    integral=taylor.integrate_zero(value,e)
    assert integral.coefficient(0)==0
    assert taylor.derivatives(value)==(1,sp.Rational(1,2),sp.Rational(-1,4),sp.Rational(3,8),sp.Rational(-15,16))
    for bad in (True,False,4.0,sp.Integer(4),0,5,None):
        with pytest.raises((ValueError,TypeError)):
            taylor.context(model.u,bad)
    with pytest.raises(ValueError):
        taylor.power(ctx.jet(),sp.Rational(1,2))


def test_public_validation_before_cache():
    physical.fixture(2)
    anchors.transformed_fixture(3)
    assert audit.controls()["rejected_inputs"]==98


def test_rational_time_and_symbolic_background_are_distinct_valid_domains():
    assert physical.background()["H"].has(model.u)
    assert physical.background(Fraction(1,2))["H"]==sp.Rational(8,5)
    with pytest.raises(TypeError):
        physical.background(0.5)
    with pytest.raises(ValueError):
        physical.background(sp.Rational(3,4))


def test_exact_raw_bound_and_jet_serialization_round_trip():
    data=verify.serialize({"bound":majorant.serialize_bound(majorant.bound(1,4,1)),
                           "jets":coefficients.jets(),"anchors":anchors.pure_electric()})
    assert json.loads(json.dumps(data))==data
    for value in (1.0,sp.Float(1),sp.oo,sp.nan):
        with pytest.raises((ValueError,TypeError)):
            verify.serialize(value)


def test_nonzero_transfer_requirement_not_silently_relaxed():
    with pytest.raises(constraints.ZeroMomentumConstraint):
        constraints.context(((1,0,0),(-1,0,0),(0,1,0),(0,-1,0)))
