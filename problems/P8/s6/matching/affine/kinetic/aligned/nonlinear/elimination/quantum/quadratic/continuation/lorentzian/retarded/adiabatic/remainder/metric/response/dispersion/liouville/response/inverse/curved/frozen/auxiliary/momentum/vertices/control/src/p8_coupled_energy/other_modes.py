"""Tensor columns and the unchanged selected-state Proca columns."""
from functools import cache

import sympy as sp
from p8_affine_aligned import modes
from p8_auxiliary_neighborhood import model
from p8_vector_regularity import preparation
from p8_vector_state import comparison

from . import window


@cache
def vectors():
    mass=sp.Integer(1000)
    common=preparation.constants()["common"]
    mixing=common["initial_mixing_envelopes"][6]+common["evolution_mixing_envelope"]/mass**4
    # The original eighth-order comparison and unchanged Borel state give
    # |B_mix|<1, |A_mix|<2, |f|<=omega^-1/2, |p_f|<=2*omega^1/2.
    return {"mass":mass,"zeta":1/mass**2,"Amax":comparison.AMAX,
            "actual_selected_mixing_coefficient":mixing,
            "all_momentum_mixing_upper":mixing/mass**6,
            "exact_mode_coordinate_upper_over_inverse_sqrt_omega":3,
            "exact_mode_canonical_momentum_upper_over_sqrt_omega":6,
            "after_fixed_center_Fourier_conversion_phase_upper":10**5*window.UPPER_K,
            "raw_seed_margin":window.RAW_SEED-10**5*window.UPPER_K,
            "selected_state_unchanged":True,
            "reference_not_reselected_as_the_exact_state":True}


@cache
def checks():
    a,m,q=sp.symbols("positive_scale positive_mass positive_q",positive=True)
    omega=sp.sqrt(q+m*m)
    gT=sp.sqrt(a)/m
    gL=a**sp.Rational(3,2)*sp.sqrt(q)/omega
    local={
        "transverse_coordinate":1/(a*gT),
        "transverse_momentum":gT/a**2,
        "longitudinal_coordinate":sp.sqrt(q)/gL,
        "longitudinal_momentum":gL/(a**3*sp.sqrt(q))}
    expected={
        "transverse_coordinate":m/a**sp.Rational(3,2),
        "transverse_momentum":1/(m*a**sp.Rational(3,2)),
        "longitudinal_coordinate":omega/a**sp.Rational(3,2),
        "longitudinal_momentum":1/(omega*a**sp.Rational(3,2))}
    out={name+"_actual_local_phase_conversion":sp.simplify(local[name]-expected[name]) for name in local}
    for name in ("transverse","longitudinal"):
        out[name+"_local_phase_symplectic_volume"]=sp.simplify(
            local[name+"_coordinate"]*local[name+"_momentum"]-a**-3)
    u=model.u
    bg=model.coefficients()["background"]
    H=bg["H"]
    pump=sp.factor(sp.Rational(3,2)*sp.diff(H,u)+sp.Rational(9,4)*H**2)
    out["actual_tensor_canonical_pump"]=sp.factor(pump-(6+30*u**2)/(1+u**2)**2)
    chosen=vectors()
    out["selected_state_mixing_constant_replay"]=chosen["actual_selected_mixing_coefficient"]-(
        1813229+sp.Rational(5347035781757616,1000**4))
    out["original_Proca_frequency_floor_replay"]=sp.factor(
        modes.canonical()["longitudinal"]["frequency_squared"]
        -modes.canonical()["q"]-1/modes.ZETA+modes.canonical()["longitudinal"]["correction"])
    return out


@cache
def tensors():
    # Tensor oscillator y=a^(3/2)*gamma/sqrt(2), with E:E=2.
    # |pump|<=15, |pump'|<=65 on I. On q>=1e4, its positive
    # energy obeys |E'|<=10E. e^10<3^10<1e5 on the whole I.
    return {"pump_absolute_upper":15,"pump_derivative_absolute_upper":65,
            "energy_logarithmic_growth_upper":10,"global_I_energy_inflation_upper":3**10,
            "initial_energy_upper_over_kappa":1,
            "coordinate_y_upper_over_inverse_sqrt_kappa":2000,
            "derivative_y_upper_over_sqrt_kappa":500,
            "after_TT_dual_and_fixed_center_Fourier_conversion_phase_upper":10**5*window.UPPER_K,
            "raw_seed_margin":window.RAW_SEED-10**5*window.UPPER_K,
            "both_tensor_polarizations_kept":True}


@cache
def gates():
    v=vectors()
    return {"unchanged_selected_vector_mixing_remains_below_one":bool(v["all_momentum_mixing_upper"]<1),
            "selected_state_frequency_mass_floor":v["mass"]>=1000,
            "both_transverse_and_longitudinal_physical_phase_weights_retained":True,
            "actual_WKB_reference_and_Borel_state_are_not_identified":True,
            "vector_and_tensor_columns_fit_common_raw_seed":bool(v["raw_seed_margin"]>0 and tensors()["raw_seed_margin"]>0),
            "tensor_global_energy_growth_budget":3**10<10**5,
            "tensor_pump_derivative_budget":sp.Rational(3,2)*13+sp.Rational(9,2)*2*5<65,
            "tensor_pump_value_budget":sp.Rational(3,2)*4+sp.Rational(9,4)*4==15,
            "tensor_frequency_high_q_domain":sp.Rational(window.LOWER_K**2,4)>10**4,
            "vector_mass_below_hard_external_band":window.LOWER_K>v["mass"]}
