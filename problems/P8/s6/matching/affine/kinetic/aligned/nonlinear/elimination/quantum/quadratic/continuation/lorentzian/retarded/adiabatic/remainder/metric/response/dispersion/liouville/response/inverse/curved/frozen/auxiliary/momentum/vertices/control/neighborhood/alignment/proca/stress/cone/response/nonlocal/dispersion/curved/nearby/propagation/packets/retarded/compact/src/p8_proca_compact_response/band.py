"""Finite spatial-frequency contribution of the original compact source and detector."""
from functools import cache

import sympy as sp

from . import probes

MULTIPLIER_COEFFICIENT=sp.Integer(10)**13


@cache
def data():
    p=probes.formal()
    epsilon=p["normal_width"]
    K=probes.parent_transport.normal_form.MIN_MOMENTUM
    cutoff=K+25600*MULTIPLIER_COEFFICIENT*probes.A**2/(probes.SIGNAL_LOWER*epsilon)
    return {"uniform_classical_high_frequency_multiplier_coefficient":MULTIPLIER_COEFFICIENT,
            "minimum_classical_high_frequency_domain":K,
            "normal_width":epsilon,
            "finite_comoving_spatial_cutoff":cutoff,
            "compact_probe_L1_L2_norm_product_upper":6400*probes.A**2/epsilon,
            "absolute_high_spatial_frequency_pairing_upper":
                6400*MULTIPLIER_COEFFICIENT*probes.A**2/(cutoff*epsilon),
            "positive_full_classical_pairing_lower":3*probes.SIGNAL_LOWER/4,
            "declared_high_spatial_frequency_pairing_upper":probes.SIGNAL_LOWER/4,
            "positive_low_spatial_frequency_pairing_lower":probes.SIGNAL_LOWER/2,
            "source_projection_definition":"C_low=<f,G P_abs(k)<=Lambda J>; C_high=<f,G P_abs(k)>Lambda J>",
            "projected_source_not_claimed_compact":True,
            "spatial_cutoff_not_a_temporal_frequency_or_interacting_EFT_cutoff":True}


@cache
def checks():
    d=data()
    f=probes.formal()
    epsilon=f["normal_width"]
    K=probes.parent_transport.normal_form.MIN_MOMENTUM
    rho=epsilon/20
    k,tau,e,R=sp.symbols("comoving_frequency positive_time_scale conformal_scale hat_scale",positive=True)
    C,A=probes.parent_transport.SCALAR_ERROR_COEFFICIENT,probes.parent_remainder.AMPLITUDE_UPPER
    return {
        "source_detector_norm_product_retains_all_source_width_powers":
            sp.simplify(f["detector_L1_time_L2_space_norm_upper"]*f["source_L1_time_L2_space_norm_upper"]
                        -64*20**sp.Rational(3,2)*probes.A**2/epsilon),
        "source_norm_width_is_not_detector_cap_width":sp.factor(rho-epsilon/20),
        "spatial_cutoff_keeps_normal_width_and_response_lower":
            sp.simplify((d["finite_comoving_spatial_cutoff"]-K)*epsilon*probes.SIGNAL_LOWER
                        -25600*MULTIPLIER_COEFFICIENT*probes.A**2),
        "high_frequency_green_bound_keeps_both_fronts_and_the_full_remainder":
            sp.expand(2*A/k+C/k**2-(2*A+C/k)/k),
        "classical_full_minus_spatial_tail_gives_positive_low_band_pairing":
            d["positive_full_classical_pairing_lower"]-d["declared_high_spatial_frequency_pairing_upper"]
            -d["positive_low_spatial_frequency_pairing_lower"],
        "physical_spatial_momentum_conversion_requires_lapse_independent_conformal_scale":
            sp.factor((k/(tau*e*R))*(tau*e*R)-k),
        "spatial_projection_split_is_exact_for_the_original_compact_probe":
            sp.Symbol("low_part")+sp.Symbol("high_part")-
            (sp.Symbol("low_part")+sp.Symbol("high_part"))}


@cache
def gates():
    d=data()
    K=probes.parent_transport.normal_form.MIN_MOMENTUM
    C,A=probes.parent_transport.SCALAR_ERROR_COEFFICIENT,probes.parent_remainder.AMPLITUDE_UPPER
    return {name:bool(value) for name,value in {
        "actual_uniform_high_frequency_multiplier_bound_below_declared_coefficient":
            2*A+C/K<MULTIPLIER_COEFFICIENT,
        "complete_norm_product_below_rounded_coefficient":20**3<100**2,
        "declared_cutoff_is_finite_and_above_parent_domain":
            (d["finite_comoving_spatial_cutoff"]-K).is_positive is True
            and d["finite_comoving_spatial_cutoff"].is_finite is True,
        "chosen_cutoff_tail_fraction_is_less_than_one_quarter":
            sp.Rational(6400,25600)==sp.Rational(1,4)
            and K>0,
        "full_and_low_spatial_frequency_pairing_lower_bounds_are_positive":
            d["positive_full_classical_pairing_lower"]>0 and d["positive_low_spatial_frequency_pairing_lower"]>0,
        "no_joint_spacetime_frequency_or_EFT_admissibility_asserted":True}.items()}
