"""Exact algebra, whole-interval arithmetic and finite-width normalization."""
from fractions import Fraction

import pytest
import sympy as sp
from flint import ctx
from p8_proca_finite_width import audit, bounds, compact, intervals, time_jets


@pytest.mark.parametrize("name",list(audit.residuals()))
def test_exact_identity(name):
    value=audit.residuals()[name]
    entries=list(value) if isinstance(value,sp.MatrixBase) else [value]
    assert all(item==0 for item in entries),name

@pytest.mark.parametrize("name",list(audit.gates()))
def test_continuous_proof_gate(name):
    assert audit.gates()[name] is True,name

def test_rejected_invalid_domains():
    assert audit.controls()["rejected_inputs"]==38

@pytest.mark.parametrize("degree",(2,4,6,8))
def test_actual_full_triangular_Riccati_recursion(degree):
    data=time_jets.build(degree)
    assert len(data["at_center"])==degree+1
    assert data["at_center"][degree]!=sp.zeros(2)

@pytest.mark.parametrize("degree",(2,3,5))
def test_exact_time_jet_inverse_and_positive_root(degree):
    T=time_jets.TimeJets(degree)
    u=intervals.model.u
    value=T.rational(1+u*u)
    inverse=T.inv(value)
    assert T.mul(value,inverse)==T.one
    root=T.sqrt(value,1)
    assert T.mul(root,root)==value

@pytest.mark.parametrize("point",(0,sp.Rational(1,200),sp.Rational(-1,100),sp.Rational(1,100)))
def test_whole_interval_rational_jets_contain_independent_derivatives(point):
    with ctx.workprec(160):
        T=intervals.BallJets(5)
        u=intervals.model.u
        expr=(1+2*u+3*u**3)/(1+u*u)
        jet=T.rational(expr)
        for n in range(6):
            exact=sp.diff(expr,u,n).subs(u,point)/sp.factorial(n)
            assert jet[n].contains(intervals.scalar(exact))

def test_interval_precision_is_scoped_and_integer_enclosures_reproducible():
    before=ctx.prec
    with ctx.workprec(71):
        actual=intervals.build(radius="1/100")
        assert ctx.prec==71
        assert actual["R_derivative_upper_integers"]==intervals.build()["R_derivative_upper_integers"]
    assert ctx.prec==before

@pytest.mark.parametrize("order",(0,1,2,3))
def test_IBP_constant_inverse_reduces_to_ordinary_product_rule(order):
    L=[2,3,5,7]
    out=bounds.ibp_coefficients(L,[2,0,0,0],order)
    assert out==[2**order*int(sp.binomial(order,j))*L[order-j] for j in range(order+1)]

def test_new_state_is_not_a_finite_truncation_at_zero_momentum():
    R=time_jets.build(2)["at_center"]
    bad=-(R[0]+R[2])/sp.I
    assert bad[1,1]==sp.Rational(-201,200)
    assert compact.center_bounds()["positive_graph_threshold"]==16

def test_whole_interval_defect_is_nonzero_and_kept():
    data=bounds.data()
    assert data["scaled_Riccati_defect_order_seven_constant"]==91215897
    assert data["normalized_packet_mode_error_constant"]>0

def test_finite_width_fourier_sampler_coefficients_are_explicit():
    data=bounds.data()
    assert data["time_IBP_Sobolev_L1_coefficients"]==[619407,63666,3078,81]
    assert data["spatial_IBP_Sobolev_L1_coefficients"]==[405837,41877,2025,54]
    assert data["inverse_derivative_bounds"]==[3,18,351,11673]

def test_actual_error_split_factor_two_in_both_constants():
    data=bounds.data()
    assert data["approximate_high_band_IBP_squared_coefficient"]==sp.exp(sp.Rational(1,25))/(160*sp.pi**3)
    Ce=data["normalized_packet_mode_error_constant"]
    assert sp.simplify(data["high_band_error_L2_squared_coefficient"]-257*Ce**2/(16*sp.pi**2*32**8))==0

def test_complete_low_band_includes_zero_and_initial_mixture():
    data=compact.transfer_bounds(Fraction(1,100))
    assert data["energy_lower"]==Fraction(1,14)
    assert data["energy_upper"]==13
    assert data["norm_log_rate"]==392

def test_exact_root_isolation_covers_cancelling_eighth_order_coefficient():
    data=compact.center_bounds()
    left,right=data["root_isolation"]
    assert left*left<17985<right*right
    assert data["coefficient_norm_upper_integers"]=={2:31,4:505,6:8752,8:574765}
