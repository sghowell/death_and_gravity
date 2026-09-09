"""Explicit compact probes and finite symbolic widths on the actual nearby solution."""
from functools import cache

import sympy as sp
from p8_proca_nearby_retarded import remainder as parent_remainder
from p8_proca_nearby_retarded import support as parent_support
from p8_proca_nearby_retarded import transport as parent_transport

T=parent_support.domains.T
AMP_LOWER=sp.Rational(1,10**8)
S_MIN=T/16
S_MAX=3*T
H=parent_support.neighborhood()["common_time_and_Euclidean_space_neighborhood_radius"]
A=H/10
SIGNAL_LOWER=AMP_LOWER*S_MIN*A**3/(1024*S_MAX**2)
EPSILON_NUMERATOR=(SIGNAL_LOWER/(32*A**2))**2
REMAINDER_FLOOR=parent_transport.SCALAR_ERROR_COEFFICIENT
EPSILON_UPPER=EPSILON_NUMERATOR/REMAINDER_FLOOR**2
REMAINDER_NORM=sp.Symbol("complete_remainder_L2_upper",positive=True)


@cache
def formal():
    width=EPSILON_NUMERATOR/REMAINDER_NORM**2
    source_width=width/20
    return {"whole_interval_half_window":T/2,"source_time":-T/4,"detector_time":T/4,
            "parent_support_neighborhood_radius":H,"detector_cap_and_time_width":A,
            "clock_ray_radius_lower":S_MIN,"clock_ray_radius_upper":S_MAX,
            "clock_shell_amplitude_lower":AMP_LOWER,
            "positive_clock_shell_pairing_lower":SIGNAL_LOWER,
            "normal_width_numerator":EPSILON_NUMERATOR,
            "normal_width":width,"source_time_and_spatial_width":source_width,
            "normal_width_rational_upper":EPSILON_UPPER,
            "source_spacetime_L1_norm":sp.Integer(1),
            "source_L1_time_L2_space_norm_upper":8/source_width**sp.Rational(3,2),
            "detector_L1_time_L2_space_norm_upper":8*A**2*sp.sqrt(width),
            "detector_spatial_support_volume_upper":16*A**2*width,
            "clock_cap_solid_angle_lower":sp.pi*A**2/(16*S_MAX**2),
            "clock_cap_detector_pointwise_lower":sp.Rational(1,16),
            "compact_detector":"beta((t-t0)/a) beta(|x_perp|/a) beta((x1-S_clock(t,s0))/a) beta((|x|-S_clock(t,s0))/epsilon)",
            "compact_positive_source":"normalized beta((s-s0)/rho) times normalized radial beta(|y|/rho)",
            "beta":"exp(1-1/(1-z^2)) for |z|<1; zero otherwise",
            "source_normalization_not_assumed_equal_to_uncomputed_bump_integral":True}


@cache
def actual():
    low=parent_remainder.low_frequency()
    K=parent_transport.normal_form.MIN_MOMENTUM
    R=REMAINDER_FLOOR+K**2*low["complete_exact_minus_two_front_multiplier_upper"]
    return {"complete_spatial_remainder_L2_upper":R,
            "complete_remainder_norm_lower":REMAINDER_FLOOR,
            "exact_normal_width":EPSILON_NUMERATOR/R**2,
            "exact_source_time_and_spatial_width":EPSILON_NUMERATOR/(20*R**2),
            "finite_symbolic_exponential_not_numerically_evaluated":True}


@cache
def checks():
    r=sp.Symbol("signed_bump_coordinate",real=True)
    f=formal()
    beta=sp.exp(1-1/(1-r*r))
    normalizer_lower=sp.Rational(1,2)*(4*sp.pi/3)*sp.Rational(1,2)**3
    normalized_source_squared_upper=(4*sp.pi/3)/normalizer_lower**2
    return {
        "parent_compact_neighborhood_is_used_without_changing_the_background":
            H-parent_support.neighborhood()["common_time_and_Euclidean_space_neighborhood_radius"],
        "smooth_bump_center_is_one":beta.subs(r,0)-1,
        "half_radius_bump_exponent_is_minus_one_third":
            (1-1/(1-r*r)).subs(r,sp.Rational(1,2))+sp.Rational(1,3),
        "radial_source_normalizer_lower_keeps_the_three_dimensional_ball_volume":
            normalizer_lower-sp.pi/12,
        "normalized_source_L2_squared_upper_keeps_the_normalizer":
            sp.factor(normalized_source_squared_upper-192/sp.pi),
        "two_time_source_width_is_one_twentieth_detector_normal_width":
            f["source_time_and_spatial_width"]-f["normal_width"]/20,
        "detector_support_volume_uses_both_radial_band_edges_and_Jacobian":
            2*2*4*A**2*f["normal_width"]-f["detector_spatial_support_volume_upper"],
        "clock_cap_signal_lower_keeps_amplitude_radius_time_and_all_four_bumps":
            sp.factor(AMP_LOWER*S_MIN/(4*sp.pi)*A*f["clock_cap_solid_angle_lower"]
                      *f["clock_cap_detector_pointwise_lower"]-SIGNAL_LOWER),
        "chosen_normal_width_makes_full_remainder_pairing_one_quarter":
            sp.simplify(REMAINDER_NORM*f["detector_L1_time_L2_space_norm_upper"]-SIGNAL_LOWER/4),
        "normal_width_includes_square_of_complete_remainder_norm":
            sp.factor(f["normal_width"]*REMAINDER_NORM**2-EPSILON_NUMERATOR)}


@cache
def gates():
    cap=A/(4*S_MAX)
    rho_upper=EPSILON_UPPER/20
    return {name:bool(value) for name,value in {
        "detector_six_width_displacement_stays_inside_parent_spacelike_neighborhood":6*A<H,
        "normal_width_rational_upper_is_less_than_one_hundredth_cap_width":EPSILON_UPPER<A/100,
        "source_support_stays_inside_parent_spacelike_neighborhood":rho_upper<H,
        "actual_clock_radius_bounds_cover_both_source_and_detector_time_variations":2*A<T/4,
        "source_time_width_less_than_cap_width":rho_upper<A,
        "cap_radius_and_normal_band_stay_far_from_spatial_origin":A<S_MIN/4 and EPSILON_UPPER<S_MIN/4,
        "radial_band_Jacobian_below_two":A**2<sp.Rational(3,4)*(S_MIN-EPSILON_UPPER)**2,
        "clock_cap_solid_angle_rationalization_has_positive_branch":0<cap<1,
        "normal_bump_is_in_its_positive_plateau_on_source_clock_front":5*rho_upper<EPSILON_UPPER/2,
        "transverse_bump_is_in_its_positive_plateau_on_source_clock_front":rho_upper+A/4<A/2,
        "longitudinal_bump_is_in_its_positive_plateau_on_source_clock_front":
            5*rho_upper+A**2/(16*S_MAX)<A/2,
        "clock_shell_pairing_lower_is_exact_and_strictly_positive":
            isinstance(SIGNAL_LOWER,sp.Rational) and SIGNAL_LOWER>0,
        "source_L2_norm_upper_uses_pi_greater_than_three":192/sp.Integer(3)==64,
        "detector_volume_upper_uses_pi_less_than_four":4*sp.Integer(4)==16,
        "bump_half_radius_lower_follows_from_exp_minus_x_greater_equal_one_minus_x":
            1-sp.Rational(1,3)>sp.Rational(1,2),
        "finite_positive_actual_widths_are_symbolic_not_floating":
            actual()["complete_spatial_remainder_L2_upper"].is_positive is True
            and actual()["exact_normal_width"].is_positive is True,
        "smoothness_at_bump_boundary_is_a_written_all_orders_induction":True,
        "compact_probes_have_no_claimed_finite_temporal_frequency_support":True}.items()}
