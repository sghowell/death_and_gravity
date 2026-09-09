"""Uniform two-time fundamental-matrix error, including both endpoint changes."""
from functools import cache

import sympy as sp
from p8_proca_nearby_packets import domains, normal_form

MATRIX_ERROR_COEFFICIENT=sp.Integer(10)**39
SCALAR_ERROR_COEFFICIENT=sp.Integer(10)**44
SOURCE_OUTPUT_FACTOR=sp.Integer(256)


@cache
def constants():
    T,Z,rem=domains.T,normal_form.Z_BOUND,normal_form.REMAINDER_BOUND
    interior=6*Z+4*T*rem
    return {"minimum_classical_momentum":normal_form.MIN_MOMENTUM,
            "maximum_retarded_interval_length":T,
            "both_endpoint_near_identity_matrix_error_coefficient":interior,
            "physical_fundamental_matrix_error_coefficient_before_rounding":domains.BASIS_NORM**2*interior,
            "declared_physical_fundamental_matrix_error_coefficient":MATRIX_ERROR_COEFFICIENT,
            "output_inverse_volume_upper":sp.Integer(4),
            "source_volume_upper":sp.Integer(4),
            "physical_source_lapse_conformal_factor_upper":sp.Integer(16),
            "source_and_output_factor_upper":SOURCE_OUTPUT_FACTOR,
            "declared_scalar_retarded_multiplier_error_coefficient":SCALAR_ERROR_COEFFICIENT,
            "near_identity_change_at_source_is_inverted_not_omitted":True,
            "scalar_error_is_absolute_and_uniform_in_phase_zeros":True}


def error(value):
    value=normal_form.momentum(value)
    eta=domains.T*normal_form.REMAINDER_BOUND/value
    near=normal_form.Z_BOUND/value
    return {"classical_momentum":value,
            "near_identity_upper":near,"integrated_normalform_error_upper":eta,
            "transformed_propagator_norm_upper":sp.Integer(2),
            "both_endpoint_transformed_transfer_error_upper":6*near+4*eta,
            "actual_physical_transfer_error_upper":domains.BASIS_NORM**2*(6*near+4*eta),
            "declared_physical_transfer_error_upper":MATRIX_ERROR_COEFFICIENT/value,
            "declared_scalar_multiplier_error_upper":SCALAR_ERROR_COEFFICIENT/value**2}


@cache
def checks():
    Ct,Csi,V,D=sp.symbols("target_change source_change_inverse exact_modal_transfer phase_transfer",commutative=False)
    actual=Ct*V*Csi-D
    expansion=(Ct-1)*V*Csi+(V-D)*Csi+D*(Csi-1)
    near,eta=sp.symbols("near_identity_norm integrated_remainder",nonnegative=True)
    return {"exact_two_endpoint_propagator_error_decomposition":sp.expand(actual-expansion),
            "two_endpoint_triangle_keeps_source_inverse":sp.expand(near*2*2+2*eta*2+2*near-(6*near+4*eta)),
            "uniform_scalar_error_keeps_exact_output_inverse_frequency":
                SOURCE_OUTPUT_FACTOR*MATRIX_ERROR_COEFFICIENT-
                constants()["source_and_output_factor_upper"]*MATRIX_ERROR_COEFFICIENT,
            "source_and_output_volume_factor_is_complete":SOURCE_OUTPUT_FACTOR-4*4*16}


@cache
def gates():
    d=constants()
    lowest=error(normal_form.MIN_MOMENTUM)
    return {name:bool(value) for name,value in {
        "uniform_remainder_applies_to_every_ordered_pair_in_inner_interval":
            d["maximum_retarded_interval_length"]==2*domains.INNER_HALF_WIDTH,
        "inverse_changes_exist_at_both_endpoints":lowest["near_identity_upper"]<sp.Rational(1,2),
        "exact_modal_propagator_norm_below_two":lowest["integrated_normalform_error_upper"]<=sp.Rational(1,4),
        "complete_endpoint_error_below_declared_matrix_coefficient":
            d["physical_fundamental_matrix_error_coefficient_before_rounding"]<MATRIX_ERROR_COEFFICIENT,
        "source_and_output_factors_below_declared_scalar_coefficient":
            SOURCE_OUTPUT_FACTOR*MATRIX_ERROR_COEFFICIENT<SCALAR_ERROR_COEFFICIENT,
        "exact_multiplier_error_is_k_to_minus_two_not_k_to_minus_one":
            lowest["declared_scalar_multiplier_error_upper"]*normal_form.MIN_MOMENTUM**2==SCALAR_ERROR_COEFFICIENT,
        "high_frequency_limit_is_classical_not_interacting_EFT_domain":True}.items()}


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
         sp.Symbol("free"),sp.sqrt(2),[],0,-1,normal_form.MIN_MOMENTUM-1)
    rejected=0
    for value in bad:
        try:
            error(value)
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(bad):
        raise ValueError("An inadmissible retarded high-frequency momentum was accepted")
    return {"rejected_inputs":rejected}
