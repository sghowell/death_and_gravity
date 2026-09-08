"""Independent matrix, sign, endpoint and domain controls."""
import pytest
import sympy as sp
from p8_vector_dispersion_inverse import block, cut, kernel, proofs, verify
from p8_vector_metric_dispersion import principal, spectral


def zero(values):
    for value in values.values():
        if isinstance(value,sp.MatrixBase):
            assert all(entry==0 for entry in value)
        else:
            assert value==0


def test_first_sheet_and_congruence_identities():
    zero(block.checks())


def test_generic_inverse_positive_Gram_and_cut_sign():
    zero(cut.generic_checks())
    zero(cut.checks())


def test_endpoint_and_causal_normalization_identities():
    zero(kernel.checks())


def test_complete_science_counts_and_audit_gates():
    values=proofs.residuals()
    zero(values)
    assert len(values)==40
    assert sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in values.values())==100
    assert len(proofs.checks())==15
    assert all(value is True for value in proofs.checks().values())


def test_strict_physical_spectral_positivity_at_interior_fraction():
    M=block.data()["spectral_matrix"].subs({block.h:1,block.z:sp.Rational(1,2)})
    assert M[0,0]>0
    assert M.det()==sp.Rational(32,6561)
    assert M==2*spectral.actual(spectral.pair("T")).subs({block.h:1,block.z:sp.Rational(1,2)})*spectral.actual(spectral.pair("T")).subs({block.h:1,block.z:sp.Rational(1,2)}).T+spectral.actual(spectral.pair("L")).subs({block.h:1,block.z:sp.Rational(1,2)})*spectral.actual(spectral.pair("L")).subs({block.h:1,block.z:sp.Rational(1,2)}).T


def test_subthreshold_endpoint_is_indefinite_not_positive_definite():
    item=block.data()
    T=item["threshold_matrix"].subs(block.h,1)
    Q=item["chart"].subs(block.h,1)
    assert (Q.T*T*Q)[0,0]>0
    assert (Q.T*T*Q)[1,1]==-sp.Rational(16,15)
    assert T.det()==-sp.Rational(526352,1476225)


def test_nonreal_dispersion_sign_has_no_scalar_sign_ambiguity():
    p=1+sp.I
    y=sp.Rational(1,2)
    weight=-p/(4*(4+(1-y*y)*p))
    assert sp.im(weight)<0
    assert sp.im(weight)==-1/abs(4+(1-y*y)*p)**2


def test_cut_inverse_sign_reversal_would_violate_positive_density():
    D=sp.diag(1,-1)
    W=sp.eye(2)
    correct=(D-sp.I*sp.pi*W).inv().applyfunc(lambda entry:sp.simplify(sp.im(entry)/sp.pi))
    wrong=(D+sp.I*sp.pi*W).inv().applyfunc(lambda entry:sp.simplify(sp.im(entry)/sp.pi))
    assert correct==sp.eye(2)/(1+sp.pi**2)
    assert wrong==-correct


def test_threshold_density_vanishes_with_nonzero_rank_one_coefficient():
    value=kernel.data()["density_over_sqrt_z_at_threshold"].subs(block.h,1)
    assert value.det()==0
    assert value.trace()>0
    assert value==block.clean(block.data()["threshold_matrix"].subs(block.h,1).inv()*block.data()["spectral_matrix"].subs({block.h:1,block.z:0})*block.data()["threshold_matrix"].subs(block.h,1).inv()/8)


def test_static_spectral_moment_is_positive_and_cannot_be_dropped():
    item=block.data()
    moment=item["static_spectral_moment"].subs(block.h,1)
    assert moment[0,0]>0
    assert moment.det()==sp.Rational(177147,99536384)
    assert item["inverse_infinity"]!=item["local_matrix"].inv()
    assert block.clean(item["inverse_infinity"]-item["local_matrix"].inv()-item["static_spectral_moment"])==sp.zeros(2)


def test_instantaneous_matrix_is_nonzero_and_cannot_have_small_interval_norm_zero():
    value=block.data()["inverse_infinity"].subs(block.h,1)
    assert value.det()==0
    assert value.trace()==sp.Rational(102255,15616)>6


def test_large_frequency_density_has_the_recorded_rank_one_coefficient():
    value=kernel.data()["density_times_log_tau_over_mass2_squared_at_infinity"]
    assert value.det()==0
    assert value[1,1]==sp.Rational(156025,119072)>0


def test_leading_log_root_is_not_in_its_large_frequency_regime():
    ell=-sp.Rational(9137,58560)
    assert ell<0
    assert sp.factor((principal.asymptotic()-2*principal.pole()*ell).det())==0
    assert principal.real_axis()["positive_diagonal_lower"]>0
    assert principal.real_axis()["determinant_upper"]<0


def test_exact_primitive_trace_bounds_at_both_interval_endpoints():
    one=kernel.moment_bound(1)
    end=kernel.moment_bound(sp.Rational(125,64))
    assert one["primitive_spectral_norm_upper"]==sp.Rational(106392093,24884096)
    assert end["primitive_spectral_norm_upper"]==sp.Rational(1394087584389,101925257216)
    assert end["instantaneous_spectral_norm_upper"]==sp.Rational(1553463015,63963136)
    assert end["primitive_spectral_norm_upper"]>one["primitive_spectral_norm_upper"]


def test_positive_exact_domain_rejects_bool_float_zero_and_symbol_after_population():
    block.scale(1)
    kernel.moment_bound(1)
    assert verify.controls()["rejected_inputs"]==38
    for value in (True,sp.true,1.0,0,"1",None,[1],sp.Symbol("h")):
        with pytest.raises((TypeError,ValueError)):
            block.scale(value)


def test_variable_h_is_congruence_not_scalar_rescaling_or_pulled_through_convolution():
    item=block.data()
    S=block.scale(sp.Rational(125,64))
    assert item["local_matrix"].subs(block.h,sp.Rational(125,64))==S*item["local_matrix"].subs(block.h,1)*S
    assert item["local_matrix"].subs(block.h,sp.Rational(125,64))!=item["local_matrix"].subs(block.h,1)/(sp.Rational(125,64)**2)
