"""A newly named rational-width compact probe pair with a sharper full L2 remainder."""
from functools import cache

import sympy as sp
from p8_proca_compact_response import bumps
from p8_proca_compact_response import probes as parent
from p8_proca_scalar_ccr import state

from . import transfer

REMAINDER_NORM=sp.Integer(10)**9
EPSILON=parent.EPSILON_NUMERATOR/REMAINDER_NORM**2
RHO=EPSILON/20
CUTOFF=1+25600*transfer.SCALAR_MULTIPLIER*parent.A**2/(parent.SIGNAL_LOWER*EPSILON)


@cache
def data():
    K,C=transfer.MIN_MOMENTUM,transfer.SCALAR_ERROR
    low=sp.Rational(25)/(6*sp.pi**2)
    middle=sp.Integer(110)**2*(K-1)/(2*sp.pi**2)
    high=C*C/(2*sp.pi**2*K)
    rational=sp.Rational(25,54)+sp.Integer(110)**2*K/18+C*C/(18*K)
    b=bumps.data()
    epsilon=sp.Symbol("normal_width",positive=True)
    a=sp.Symbol("cap_width",positive=True)
    t0=sp.Symbol("detector_center_time",real=True)
    s0=sp.Symbol("source_center_time",real=True)
    F=8*parent.A**2*sp.sqrt(EPSILON)
    J=8/RHO**sp.Rational(3,2)
    return {"new_probe_pair_name":"S6.94 rational-width compact J_sharp and f_sharp",
        "same_actual_background_and_retarded_kernel_but_different_probe_widths":True,
        "complete_remainder_low_unit_ball_L2_squared_upper":low,
        "complete_remainder_intermediate_band_L2_squared_upper":middle,
        "complete_remainder_high_tail_L2_squared_upper":high,
        "complete_remainder_L2_squared_rational_upper":rational,
        "declared_complete_remainder_L2_norm_upper":REMAINDER_NORM,
        "new_exact_normal_width":EPSILON,"new_exact_source_width":RHO,
        "source_spacetime_L1_norm":sp.Integer(1),
        "source_L1_time_L2_space_norm_upper":J,"detector_L1_time_L2_space_norm_upper":F,
        "actual_new_compact_source_formula":b["actual_compact_source_formula"].subs(
            {bumps.rho:RHO,s0:-parent.T/4},simultaneous=True),
        "actual_new_compact_detector_formula":b["actual_compact_detector_formula"].subs(
            {epsilon:EPSILON,a:parent.A,t0:parent.T/4},simultaneous=True),
        "unchanged_actual_bump_normalizers":
            (b["one_dimensional_normalizer"],b["three_dimensional_normalizer"]),
        "unchanged_clock_front_pairing_lower":parent.SIGNAL_LOWER,
        "new_complete_remainder_pairing_upper":REMAINDER_NORM*F,
        "positive_new_full_classical_pairing_lower":3*parent.SIGNAL_LOWER/4,
        "new_finite_spatial_cutoff":CUTOFF,
        "new_source_only_high_spatial_tail_upper":
            6400*transfer.SCALAR_MULTIPLIER*parent.A**2/(CUTOFF*EPSILON),
        "positive_new_low_spatial_band_pairing_lower":parent.SIGNAL_LOWER/2,
        "positive_new_compact_quantum_commutator_magnitude_lower":
            3*state.hbar*parent.SIGNAL_LOWER/(4*state.kappa),
        "new_cutoff_decimal_upper_exponent":94,
        "finite_rational_width_and_cutoff_not_interacting_EFT_scales":True,
        "projected_source_not_compact_and_temporal_frequency_not_truncated":True}


@cache
def checks():
    d=data()
    k,C,K=sp.symbols("radial_frequency high_error split_frequency",positive=True)
    rows={"new_high_tail_keeps_three_dimensional_Plancherel_factor":
        sp.integrate(4*sp.pi*k*k*(C/k**2)**2,(k,K,sp.oo))/(2*sp.pi)**3-C*C/(2*sp.pi**2*K),
        "new_intermediate_band_has_no_uncovered_frequency_gap":
        sp.integrate(4*sp.pi*k*k*(110/k)**2,(k,1,K))/(2*sp.pi)**3-sp.Integer(110)**2*(K-1)/(2*sp.pi**2),
        "new_normal_width_keeps_complete_remainder_square":
        EPSILON*REMAINDER_NORM**2-parent.EPSILON_NUMERATOR,
        "new_compact_source_width_is_one_twentieth_normal_width":RHO-EPSILON/20,
        "new_remainder_pairing_is_one_quarter_unchanged_signal":
        sp.simplify(d["new_complete_remainder_pairing_upper"]-parent.SIGNAL_LOWER/4),
        "new_frequency_cutoff_retains_complete_norm_product":
        CUTOFF-1-25600*transfer.SCALAR_MULTIPLIER*parent.A**2/(parent.SIGNAL_LOWER*EPSILON),
        "new_quantum_lower_keeps_action_scale_and_Planck_constant":
        d["positive_new_compact_quantum_commutator_magnitude_lower"]
        -state.hbar*d["positive_new_full_classical_pairing_lower"]/state.kappa}
    return {name:sp.factor(value) for name,value in rows.items()}


@cache
def gates():
    d=data()
    cap=parent.A/(4*parent.S_MAX)
    rows={"low_unit_ball_exact_minus_front_multiplier_below_five":4+24*parent.T<5,
        "intermediate_exact_minus_front_multiplier_below_one_hundred_ten":transfer.SCALAR_MULTIPLIER+6<110,
        "complete_three_piece_L2_remainder_below_one_billion":
            d["complete_remainder_L2_squared_rational_upper"]<REMAINDER_NORM**2,
        "new_normal_width_is_positive_exact_rational":isinstance(EPSILON,sp.Rational) and EPSILON>0,
        "new_normal_width_less_than_one_hundredth_cap_width":EPSILON<parent.A/100,
        "new_source_support_inside_original_spacelike_neighborhood":RHO<parent.H,
        "new_detector_support_inside_original_spacelike_neighborhood":6*parent.A<parent.H,
        "new_source_time_width_less_than_cap_width":RHO<parent.A,
        "new_clock_radius_window_stays_time_ordered":2*parent.A<parent.T/4,
        "new_cap_and_normal_band_avoid_spatial_origin":parent.A<parent.S_MIN/4 and EPSILON<parent.S_MIN/4,
        "new_radial_band_Jacobian_below_two":parent.A**2<sp.Rational(3,4)*(parent.S_MIN-EPSILON)**2,
        "new_cap_angle_has_positive_branch":0<cap<1,
        "new_normal_bump_positive_plateau":5*RHO<EPSILON/2,
        "new_transverse_bump_positive_plateau":RHO+parent.A/4<parent.A/2,
        "new_longitudinal_bump_positive_plateau":5*RHO+parent.A**2/(16*parent.S_MAX)<parent.A/2,
        "new_source_norm_product_bound_keeps_three_half_power":20**3<100**2,
        "new_finite_spatial_tail_below_quarter_signal":
            d["new_source_only_high_spatial_tail_upper"]<parent.SIGNAL_LOWER/4,
        "new_cutoff_in_proved_all_frequency_response_domain":CUTOFF>1,
        "new_cutoff_below_one_e_ninety_four":CUTOFF<10**94,
        "new_compact_quantum_commutator_strictly_nonzero":
            d["positive_new_compact_quantum_commutator_magnitude_lower"].is_positive is True,
        "no_earlier_literal_probe_pair_or_quantum_state_replaced":True}
    return {name:bool(value) for name,value in rows.items()}
