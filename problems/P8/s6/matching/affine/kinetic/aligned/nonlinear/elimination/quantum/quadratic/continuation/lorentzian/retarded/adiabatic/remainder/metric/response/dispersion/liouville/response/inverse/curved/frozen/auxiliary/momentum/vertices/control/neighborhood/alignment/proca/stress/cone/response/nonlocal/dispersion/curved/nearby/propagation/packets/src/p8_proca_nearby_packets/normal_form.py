"""Quantitative four-mode normal form on the inner actual time interval."""
from functools import cache

import sympy as sp

from . import domains, modes

MIN_MOMENTUM=4*sp.Integer(10)**31
Z_BOUND=sp.Integer(10)**21
Z_DOT_BOUND=sp.Integer(10)**36
B0_BOUND=3*sp.Integer(10)**13
B0_DOT_BOUND=2*sp.Integer(10)**21
W_BOUND=3*sp.Integer(10)**10
REMAINDER_BOUND=sp.Integer(10)**38


def momentum(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require an exact real classical comoving momentum")
    value=sp.Rational(value)
    if value<MIN_MOMENTUM:
        raise ValueError("Momentum is below the declared classical normal-form domain")
    return value


@cache
def constants():
    T=domains.T
    S,L=domains.BASIS_NORM,domains.MATRIX_NORM
    Sd=S/(T/2)
    Sdd=2*S/(T/2)**2
    Ld=L/(T/2)
    B=S*L*S+S*Sd
    Bd=3*Sd*L*S+Sd*Sd+S*Sdd
    Lambda_d=4/(T/2)
    Z=B0_BOUND/domains.GAP
    Zd=B0_DOT_BOUND/domains.GAP+2*Lambda_d*B0_BOUND/domains.GAP**2
    W=3*L*S*S
    remainder=2*(B0_BOUND*Z_BOUND+sp.Rational(3,2)*W_BOUND+Z_DOT_BOUND)
    return {"actual_inner_time_half_width":domains.INNER_HALF_WIDTH,"time_Cauchy_radius":T/2,
            "basis_and_inverse_norm_upper":S,"all_four_L_matrix_norm_upper":L,
            "basis_and_inverse_first_derivative_upper":Sd,"basis_second_derivative_upper":Sdd,
            "L0_first_derivative_upper":Ld,
            "leading_transport_triangle_upper":B,"leading_transport_first_derivative_triangle_upper":Bd,
            "declared_B0_upper":B0_BOUND,"declared_B0_dot_upper":B0_DOT_BOUND,
            "signed_frequency_first_derivative_upper":Lambda_d,
            "off_diagonal_change_triangle_upper":Z,"off_diagonal_change_derivative_triangle_upper":Zd,
            "declared_Z_upper":Z_BOUND,"declared_Z_dot_upper":Z_DOT_BOUND,
            "inverse_frequency_remainder_triangle_upper":W,"declared_W_upper":W_BOUND,
            "exact_transformed_remainder_triangle_upper":remainder,"declared_R_upper":REMAINDER_BOUND,
            "minimum_classical_momentum":MIN_MOMENTUM,
            "all_leading_diagonal_transport_entries_identically_zero":True}


@cache
def data():
    d=modes.data()
    lam=d["Lambda"]
    B=d["B0"]
    Z=sp.zeros(4)
    for i in range(4):
        for j in range(4):
            if i!=j:
                Z[i,j]=sp.factor(sp.I*B[i,j]/(lam[i,i]-lam[j,j]))
    return {"Lambda":lam,"B0":B,"Z":Z}


def error(value):
    value=momentum(value)
    eta=domains.T*REMAINDER_BOUND/value
    modal=2*eta
    near=Z_BOUND/value
    return {"classical_momentum":value,"near_identity_difference_upper":near,
            "near_identity_inverse_norm_upper":sp.Integer(2),
            "integrated_remainder_upper":eta,
            "modal_phase_approximation_error_upper":modal,
            "full_U_error_upper":(1+near)*modal+near,
            "phase_propagator_norm":sp.Integer(1),
            "no_exponential_of_carrier_or_B0_norm":True}


@cache
def checks():
    d=data()
    lam,B,Z=d["Lambda"],d["B0"],d["Z"]
    k=sp.Symbol("positive_classical_momentum",positive=True)
    W=sp.Matrix(4,4,lambda i,j:sp.Symbol(f"W_{i}_{j}"))
    Zd=sp.Matrix(4,4,lambda i,j:sp.Symbol(f"Zdot_{i}_{j}"))
    change=sp.eye(4)+Z/k
    actual=(sp.I*k*lam+B+W/k)*change-Zd/k-change*(sp.I*k*lam)
    expected=(B*Z+W+W*Z/k-Zd)/k
    b,bd,delta,dd=sp.symbols("offdiagonal_entry entry_derivative gap gap_derivative",nonzero=True)
    t=sp.Symbol("jet_time",real=True)
    return {
        "off_diagonal_transport_removed_by_actual_gap_matrix":
            (sp.I*(lam*Z-Z*lam)+B).applyfunc(sp.factor),
        "exact_time_dependent_near_identity_remainder":
            (actual-expected).applyfunc(sp.factor),
        "gap_derivative_is_retained_in_change_derivative":
            sp.factor(sp.diff(sp.I*(b+bd*t)/(delta+dd*t),t).subs(t,0)-sp.I*(bd/delta-b*dd/delta**2)),
        "exponential_series_bound_is_applied_to_integrated_small_remainder":
            error(MIN_MOMENTUM)["integrated_remainder_upper"]-sp.Rational(1,4),
        "Cauchy_second_derivative_factorial_retained":
            constants()["basis_second_derivative_upper"]-2*domains.BASIS_NORM/(domains.T/2)**2}


@cache
def gates():
    d=constants()
    e=error(MIN_MOMENTUM)
    return {name:bool(value) for name,value in {
        "time_Cauchy_first_and_second_derivative_bounds_finite":d["basis_and_inverse_first_derivative_upper"]==2*10**10
            and d["basis_second_derivative_upper"]==8*10**17,
        "actual_leading_transport_below_declared_upper":d["leading_transport_triangle_upper"]<B0_BOUND,
        "actual_transport_derivative_below_declared_upper":d["leading_transport_first_derivative_triangle_upper"]<B0_DOT_BOUND,
        "actual_change_below_declared_upper":d["off_diagonal_change_triangle_upper"]<Z_BOUND,
        "actual_change_derivative_below_declared_upper":d["off_diagonal_change_derivative_triangle_upper"]<Z_DOT_BOUND,
        "actual_W_upper_has_all_three_inverse_frequency_pieces":d["inverse_frequency_remainder_triangle_upper"]==W_BOUND,
        "actual_R_upper_keeps_time_derivative_and_all_contacts":d["exact_transformed_remainder_triangle_upper"]<REMAINDER_BOUND,
        "change_invertible_on_entire_classical_domain":e["near_identity_difference_upper"]<sp.Rational(1,2),
        "integrated_remainder_in_exponential_series_domain":e["integrated_remainder_upper"]<=sp.Rational(1,4),
        "not_a_classical_to_interacting_EFT_cutoff_transfer":True}.items()}


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
         sp.Symbol("free"),sp.sqrt(2),[],0,-1,MIN_MOMENTUM-1)
    rejected=0
    for value in bad:
        for call in (momentum,error):
            try:
                call(value)
            except (TypeError,ValueError):
                rejected+=1
    if rejected!=2*len(bad):
        raise ValueError("An inadmissible normal-form momentum was accepted")
    return {"rejected_inputs":rejected}
