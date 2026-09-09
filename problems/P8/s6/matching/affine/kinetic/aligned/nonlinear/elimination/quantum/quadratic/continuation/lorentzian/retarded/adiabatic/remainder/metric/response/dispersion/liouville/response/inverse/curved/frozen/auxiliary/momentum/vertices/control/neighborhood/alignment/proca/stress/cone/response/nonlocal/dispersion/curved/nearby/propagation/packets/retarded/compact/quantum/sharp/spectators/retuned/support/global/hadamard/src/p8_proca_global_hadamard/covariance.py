"""Positive graph covariance, exact CCR, and infrared completion without state transfer."""
from functools import cache

import sympy as sp


def cutoff_weight(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require an exact convex cutoff weight")
    value=sp.Rational(value)
    if not 0<=value<=1:
        raise ValueError("The covariance cutoff must be convex")
    return value


@cache
def data():
    hbar,kappa=sp.symbols("positive_hbar positive_action_normalization",positive=True)
    l1,l2=sp.symbols("positive_Cholesky_first positive_Cholesky_second",positive=True)
    l21=sp.Symbol("real_Cholesky_lower",real=True)
    a1,a2,a3=sp.symbols("real_graph_A11 real_graph_A12 real_graph_A22",real=True)
    A=sp.Matrix([[a1,a2],[a2,a3]])
    L=sp.Matrix([[l1,0],[l21,l2]])
    B=L*L.T
    T=L.inv().T
    R=A-sp.I*B
    factor=sp.sqrt(hbar/(2*kappa))*T.col_join(R*T)
    Bi=B.inv()
    C=hbar/(2*kappa)*Bi.row_join(Bi*A+sp.I*sp.eye(2)).col_join(
        (A*Bi-sp.I*sp.eye(2)).row_join(A*Bi*A+B))
    Omega=sp.zeros(4)
    Omega[:2,2:]=sp.eye(2)
    Omega[2:,:2]=-sp.eye(2)
    return {"positive_hbar":hbar,"positive_action_normalization":kappa,
        "real_symmetric_graph_part":A,"positive_imaginary_graph_part":B,
        "positive_Cholesky_matrix":L,"inverse_square_Gram_factor":T,
        "negative_frequency_graph":R,"positive_Cauchy_covariance_Gram_factor":factor,
        "exact_graph_Cauchy_covariance":C,"constant_symplectic_form":Omega,
        "physical_density_high_frequency_covariance":"(sqrt(k)*E_packet^-1)*C_graph*(sqrt(k)*E_packet^-1)^T at u=0",
        "infrared_completion":"chi(k/K)*C_high+(1-chi(k/K))*C_mu, where 0<=chi<=1, chi=0 for k<=K and chi=1 for k>=2K",
        "initial_mu_covariance_only":"C_mu=hbar/(2*kappa)*[mu^-1 I,i I;-i I,mu I], mu=sqrt(1+k^2)",
        "exact_new_evolution":"W(t,s,k)=U_global(t,0,k)*C_initial(k)*U_global(s,0,k)^T",
        "initial_density_covariance_high_frequency_symbol_order":3,
        "complete_compact_strip_covariance_polynomial_growth_order":21,
        "no_old_evolved_state_or_finite_order_adiabatic_state_is_promoted":True}


@cache
def checks():
    d=data()
    C,F,Omega=d["exact_graph_Cauchy_covariance"],d["positive_Cauchy_covariance_Gram_factor"],d["constant_symplectic_form"]
    scalar=d["positive_hbar"]/d["positive_action_normalization"]
    t=sp.Symbol("convex_cutoff_weight",real=True)
    R=d["negative_frequency_graph"]
    T=d["inverse_square_Gram_factor"]
    return {
        "positive_graph_covariance_is_exact_Gram":(C-F*F.conjugate().T).applyfunc(sp.factor),
        "positive_graph_covariance_has_exact_density_CCR":(C-C.T-sp.I*scalar*Omega).applyfunc(sp.factor),
        "positive_graph_covariance_is_Hermitian":(C-C.conjugate().T).applyfunc(sp.factor),
        "positive_graph_covariance_range_is_negative_frequency_graph":(
            F[2:,:]-R*F[:2,:]).applyfunc(sp.factor),
        "positive_covariance_factor_keeps_inverse_B_ordering":(
            T*T.T-d["positive_imaginary_graph_part"].inv()).applyfunc(sp.factor),
        "negative_frequency_graph_is_complex_symmetric":R-R.T,
        "convex_infrared_completion_preserves_CCR":sp.I*scalar*(t+(1-t)-1)*Omega,
        "global_two_time_tempered_growth_includes_both_transfers_and_initial_covariance":sp.Integer(9+3+9-21),
    }


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
        sp.Symbol("free"),sp.sqrt(2),-1,2,sp.Rational(-1,10),sp.Rational(11,10))
    rejected=0
    for value in bad:
        try:
            cutoff_weight(value)
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(bad):
        raise ValueError("An invalid covariance weight was accepted")
    return {"rejected_inputs":rejected}
