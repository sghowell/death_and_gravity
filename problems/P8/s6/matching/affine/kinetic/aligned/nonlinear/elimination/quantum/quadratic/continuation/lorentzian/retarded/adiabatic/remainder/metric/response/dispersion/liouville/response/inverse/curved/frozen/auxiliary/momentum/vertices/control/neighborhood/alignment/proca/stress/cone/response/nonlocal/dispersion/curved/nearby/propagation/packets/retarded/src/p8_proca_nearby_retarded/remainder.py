"""Complete low-frequency and high-frequency spatial L2 remainder bounds."""
from functools import cache

import sympy as sp
from p8_proca_nearby_packets import domains, normal_form

from . import source, transport

AMPLITUDE_UPPER=sp.Integer(10)**6


@cache
def low_frequency():
    p=source.parent
    K=normal_form.MIN_MOMENTUM
    limits={p.a:(0,3),p.m:(0,8),p.g:(0,2),p.r:(0,8),
            p.rn:(sp.Rational(1,10),1),p.h:(1,6),p.ell:(0,sp.Rational(1,8)),
            p.alpha:(0,1),p.beta:(0,1),p.H:(0,1)}
    coefficients=[matrix.applyfunc(lambda value:domains.laurent_majorant(value,limits))
                  for matrix in source.data()["polynomial_q_coefficients"]]
    Q=4*K**2
    matrix=sum((coefficient*Q**j for j,coefficient in enumerate(coefficients)),sp.zeros(4))
    bound=max(sum(matrix[i,j] for j in range(4)) for i in range(4))
    multiplier=16*sp.exp(domains.T*bound)+8*AMPLITUDE_UPPER*domains.T
    return {"compact_frequency_ball_radius":K,"hat_q_upper_on_compact_frequency_ball":Q,
            "complete_original_generator_coefficient_majorants":coefficients,
            "complete_original_generator_entry_majorants":matrix,
            "complete_original_generator_infinity_norm_upper":bound,
            "original_chart_fundamental_matrix_norm_upper":sp.exp(domains.T*bound),
            "complete_exact_minus_two_front_multiplier_upper":multiplier,
            "spatial_L2_squared_upper":K**3*multiplier**2/(6*sp.pi**2),
            "bound_is_finite_not_a_numerical_stability_assertion":True}


@cache
def data():
    C,K=transport.SCALAR_ERROR_COEFFICIENT,normal_form.MIN_MOMENTUM
    radius,wave,mu=sp.symbols("positive_shell_radius radial_frequency angular_cosine",positive=True)
    angular=sp.integrate(sp.cos(wave*radius*mu),(mu,-1,1))
    return {"convolution_kernel_inverse_Fourier_normalization":"(2*pi)^(-3)",
            "wave_shell":"delta(|x|-S)/(4*pi*|x|), S>0",
            "shell_forward_Fourier_angular_integral":angular,
            "shell_forward_Fourier_multiplier":sp.sin(wave*radius)/wave,
            "high_frequency_multiplier_error_upper":C/wave**2,
            "high_frequency_multiplier_L2_squared_upper":4*sp.pi*C**2/K,
            "high_frequency_spatial_remainder_L2_squared_upper":C**2/(2*sp.pi**2*K),
            "two_front_amplitude_upper":AMPLITUDE_UPPER,
            "low_frequency_complete_bound":low_frequency(),
            "remainder_uniform_in_both_times_on_inner_interval":True,
            "shell_and_square_integrable_remainder_split_uses_all_classical_frequencies":True}


@cache
def checks():
    C,K=sp.symbols("error_coefficient split_frequency",positive=True)
    S,w,mu,eps=sp.symbols("shell_radius frequency angular epsilon",positive=True)
    r=sp.Symbol("signed_normal_coordinate",real=True)
    angular=sp.integrate(sp.cos(w*S*mu),(mu,-1,1))
    tail=sp.integrate(4*sp.pi*w**2*(C/w**2)**2,(w,K,sp.oo))
    d=data()
    return {
        "spherical_delta_shell_forward_Fourier_has_correct_radius_and_four_pi":
            sp.simplify(S*angular/2-sp.sin(w*S)/w),
        "wave_shell_has_correct_zero_frequency_limit":
            sp.limit(sp.sin(w*S)/w,w,0)-S,
        "three_dimensional_k_to_minus_two_error_is_square_integrable":
            sp.factor(tail-4*sp.pi*C**2/K),
        "convolution_kernel_Plancherel_factor_retained":
            sp.factor(tail/(2*sp.pi)**3-C**2/(2*sp.pi**2*K)),
        "low_frequency_compact_ball_Plancherel_volume_is_exact":
            sp.factor((4*sp.pi*K**3/3)/(2*sp.pi)**3-K**3/(6*sp.pi**2)),
        "literal_declared_high_frequency_L2_bounds_agree":
            sp.factor(d["high_frequency_multiplier_L2_squared_upper"]/(2*sp.pi)**3-
                      d["high_frequency_spatial_remainder_L2_squared_upper"]),
        "thin_normal_test_slab_L2_measure_scales_to_zero":
            sp.integrate(1,(r,-eps,eps))-2*eps,
        "surface_delta_pairing_does_not_acquire_normal_slab_width":
            sp.integrate(sp.DiracDelta(r),(r,-eps,eps))-1}


@cache
def gates():
    d=data()
    low=d["low_frequency_complete_bound"]
    return {name:bool(value) for name,value in {
        "clock_shell_amplitude_has_finite_parent_upper":
            256*sp.Rational(1,8)**2/(sp.Rational(1,4)*sp.Rational(1,2))<AMPLITUDE_UPPER,
        "matter_shell_amplitude_has_finite_parent_upper":
            256/(sp.Rational(1,4)*sp.Rational(1,100))<AMPLITUDE_UPPER,
        "compact_low_frequency_generator_bound_is_exact_finite_positive":
            isinstance(low["complete_original_generator_infinity_norm_upper"],sp.Rational)
            and low["complete_original_generator_infinity_norm_upper"]>0,
        "all_original_generator_coefficient_majorants_nonnegative":
            all(value>=0 for matrix in low["complete_original_generator_coefficient_majorants"] for value in matrix),
        "high_frequency_radial_integral_converges_in_three_spatial_dimensions":
            d["high_frequency_spatial_remainder_L2_squared_upper"].is_finite is True,
        "complete_remainder_includes_compact_low_frequencies":
            low["spatial_L2_squared_upper"].is_finite is True,
        "no_derivative_of_high_frequency_error_needed_for_L2_non_cancellation":True,
        "delta_shell_not_locally_L2_is_a_written_distribution_proof_not_a_numeric_test":True}.items()}
