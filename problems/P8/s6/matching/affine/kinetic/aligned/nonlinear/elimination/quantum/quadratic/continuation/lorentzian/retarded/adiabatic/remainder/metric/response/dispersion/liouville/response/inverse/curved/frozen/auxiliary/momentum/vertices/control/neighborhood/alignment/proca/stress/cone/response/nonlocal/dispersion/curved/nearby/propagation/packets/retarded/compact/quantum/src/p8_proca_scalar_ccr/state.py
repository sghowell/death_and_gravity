"""An explicit positive tempered Gaussian covariance; no Hadamard claim."""
from functools import cache

import sympy as sp

from . import canonical

mu,hbar,kappa=sp.symbols("positive_initial_frequency Planck_constant overall_action_normalization",positive=True)


def frequency_weight(momentum_squared):
    """Exact smooth all-momentum positive weight; no dispersion or Hadamard assertion."""
    if isinstance(momentum_squared,bool) or not isinstance(momentum_squared,(int,sp.Expr)):
        raise TypeError("The squared momentum must be an exact nonnegative finite expression")
    value=sp.sympify(momentum_squared)
    if (value.has(sp.Float) or value.is_finite is not True
            or value.is_real is not True or value.is_nonnegative is not True):
        raise ValueError("The squared momentum must be an exact nonnegative finite expression")
    return sp.sqrt(1+value)


@cache
def data():
    I=sp.eye(2)
    gram=I/sp.sqrt(mu)
    gram=gram.col_join(-sp.I*sp.sqrt(mu)*I)*sp.sqrt(hbar/(2*kappa))
    covariance=hbar/(2*kappa)*((I/mu).row_join(sp.I*I).col_join((-sp.I*I).row_join(mu*I)))
    return {"initial_positive_frequency":"mu(k)=sqrt(1+|k|^2)",
            "initial_reference_time":-sp.Rational(1,2*10**7),
            "ordered_Cauchy_covariance":covariance,
            "positive_Cauchy_covariance_Gram_factor":gram,
            "symmetric_Cauchy_covariance":hbar/(2*kappa)*sp.diag(1/mu,1/mu,mu,mu),
            "commutator_normalization":hbar/kappa,
            "new_scalar_Gaussian_choice_not_a_transfer_of_old_Proca_state":True,
            "tempered_linear_fields_and_Wick_n_point_state_only":True,
            "Hadamard_stress_tensor_and_interacting_quantum_matching_not_asserted":True}


@cache
def checks():
    d=data()
    C,B=d["ordered_Cauchy_covariance"],d["positive_Cauchy_covariance_Gram_factor"]
    Omega=canonical.data()["constant_symplectic_form"]
    Ut=sp.Matrix(4,4,lambda i,j:sp.Symbol(f"target_transfer_{i}_{j}",real=True))
    Us=sp.Matrix(4,4,lambda i,j:sp.Symbol(f"source_transfer_{i}_{j}",real=True))
    transport_difference=Ut*C*Us.T-(Us*C*Ut.T).T
    rows={
        "explicit_all_momentum_initial_frequency_has_no_zero_mode_pole":
            frequency_weight(sp.Symbol("squared_momentum",nonnegative=True))**2
            -1-sp.Symbol("squared_momentum",nonnegative=True),
        "initial_frequency_is_exactly_one_at_zero_momentum":frequency_weight(0)-1,
        "initial_ordered_covariance_is_a_positive_Gram_matrix":C-B*B.conjugate().T,
        "initial_covariance_has_exact_canonical_antisymmetric_part":C-C.T-sp.I*hbar*Omega/kappa,
        "initial_ordered_covariance_is_Hermitian":C-C.conjugate().T,
        "symmetric_covariance_keeps_both_inverse_frequency_weights":(C+C.T)/2-d["symmetric_Cauchy_covariance"],
        "two_time_commutator_depends_only_on_CCR_not_symmetric_state_covariance":
            transport_difference-sp.I*hbar*Ut*Omega*Us.T/kappa}
    return {name:value.applyfunc(sp.factor) if isinstance(value,sp.MatrixBase) else sp.factor(value)
            for name,value in rows.items()}


def controls():
    bad=[True,False,"1",1.0,sp.Float(1),None,[],{},-1,-sp.Rational(1,2),
         sp.oo,-sp.oo,sp.zoo,sp.nan,sp.I,sp.Symbol("unconstrained"),
         sp.Symbol("only_real",real=True)]
    rejected=0
    for value in bad:
        try:
            frequency_weight(value)
        except (TypeError,ValueError):
            rejected+=1
        else:
            raise ValueError("An invalid initial Gaussian frequency input was accepted")
    return {"rejected_inputs":rejected,"zero_momentum_is_included":frequency_weight(0)==1,
            "exact_positive_Gaussian_weight_not_a_physical_dispersion_relation":True}
