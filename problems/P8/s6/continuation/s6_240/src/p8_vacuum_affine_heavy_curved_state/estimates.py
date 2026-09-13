"""Explicit all-momentum and five-derivative heavy-state stress bounds."""

import math
from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_parent_one_loop import germs

from . import quantum


@cache
def data():
    R = s.Rational
    gates = {}
    gates["complex_Re_inverse_scale_fourth"] = (
        (1 - 6 * R(1, 30) ** 2) / (1 + R(1, 30) ** 2) ** 2
    ) * R(16, 33) ** 4 > R(1, 32)
    gates["complex_v_lower"] = R(4095, 4096) > R(99, 100)
    gates["complex_v_upper"] = R(8321, 4096) < R(33, 16)
    gates["complex_v_phase_ratio"] = R(130, 4095) < R(1, 30)
    gates["complex_inverse_scale_power"] = R(100, 99) ** 6 < 2
    gates["complex_H_upper"] = 4 * R(65, 64) * R(100, 99) < 5
    gates["complex_H_prime_upper"] = 4 * R(8321, 4096) * R(100, 99) ** 2 < 12
    gates["complete_U_upper"] = R(3, 2) * 12 + R(9, 4) * 25 < 100
    gates["complete_lambda_prime_upper"] = 12 * 2 + 2 * 25 * 2 * 3 < 400
    gates["complete_base_Riccati_upper"] = 100 + R(400, 2) + R(100, 4) < 400
    loss = R(1, 1024)
    gates["log_map_difference_constant"] = 3 / loss**2 + 10 / loss < 10**7
    gates["normalized_Riccati_map_constant"] = 4 * 10**7 < 10**8
    gates["auxiliary_epsilon_disk_contraction"] = R(10**8, 10**16) == R(1, 10**8)
    gates["seven_iterations_leave_complex_room"] = R(1, 64) - 7 * loss - 10 * R(
        1, 4096
    ) > R(1, 256)
    gates["full_W6_defect_constant"] = 192 * (10**8) ** 7 < 10**60
    gates["actual_mass_floor"] = germs.MASS2 > 10**196
    gates["actual_mass_ceiling"] = germs.MASS2 < 10**198
    gates["small_full_Volterra_exponent"] = R(4 * 10**60, 10 ** (98 * 13)) < R(1, 2)
    gates["exact_mode_coefficient_error"] = 8 * 10**60 < 10**62

    weight = [216 * 86**j * math.factorial(j) ** 2 for j in range(11)]
    gates["sampling_weight_derivatives_through10"] = max(weight) < 10**36
    gates["all_ten_IBP_coefficients"] = 100 * 8194**10 * max(weight) < 10**79
    gates["full_exact_mode_c2_error"] = 10**70 + 10**79 < 10**80
    gates["complete_c2_over_c1"] = 128 * 10**80 < 10**83
    gates["actual_squeezing_ratio_small"] = R(10**83, 10 ** (98 * 11)) < R(1, 2)

    C = 10**30
    Y = [20 * math.factorial(j) * (2 * C) ** j for j in range(7)]
    F = [
        10**85 * sum(math.comb(j, r) * Y[j - r] for r in range(j + 1)) for j in range(6)
    ]
    Z = [10**65]
    for j in range(5):
        Z.append(C * sum(math.comb(j, r) * Z[j - r] for r in range(j + 1)) + F[j])

    def quadratic(left, right, j):
        return C * sum(
            math.factorial(j)
            // (math.factorial(r) * math.factorial(a) * math.factorial(b))
            * left[a]
            * right[b]
            for r in range(j + 1)
            for a in range(j - r + 1)
            for b in [j - r - a]
        )

    exact_difference = [2 * quadratic(Z, Y, j) for j in range(6)]
    sle_difference = [4 * 10**83 * quadratic(Y, Y, j) for j in range(6)]
    cauchy_tail = 2 * 10**4 * (10**8) ** 6 * 10**20
    gates["all_complete_exact_mode_stress_jets"] = max(exact_difference) < 10**300
    gates["all_complete_SLE_stress_jets"] = max(sle_difference) < 10**300
    gates["complete_W6_adiabatic_tail_all_five_jets"] = cauchy_tail < 10**74
    gates["summed_full_nonlocal_subtraction_mode_bound"] = (
        max(exact_difference) + max(sle_difference) + cauchy_tail < 10**301
    )
    gates["full_radial_integral_constant"] = R(4**3, 3 * 2 * 9) < 2
    gates["complete_renormalized_state_integral_constant"] = 2 * 10**301 < 10**310

    for j, cap in ((0, 5), (1, 6), (2, 12), (3, 30)):
        gates[f"full_H_jet_complex_bound_{j}"] = (
            4 * math.factorial(j) * R(64, 63) ** (j + 1) < cap
        )
    gates["finite_curvature_mass_term"] = R(2, 3) * (3 * 25 + 2 * 6) * 463 < 40000
    gates["finite_curvature_square_energy"] = (
        462 * (6 * 25 * 6 + 2 * 5 * 12 + 6**2) < 10**6
    )
    gates["finite_curvature_square_pressure"] = (
        R(462, 3) * (18 * 25 * 6 + 12 * 5 * 12 + 9 * 6**2 + 2 * 30) < 10**6
    )
    gates["local_nonconstant_five_jet_factor"] = math.factorial(5) * 64**5 < 10**12
    n = germs.MASS2
    local = (462 * n * n + 10**12 * (40000 * n + 10**6)) / 576
    state = R(10**310, 1) / n
    gates["full_local_heat_stress_all_five_jets"] = local < 10**399
    gates["full_state_and_subtraction_below_local_budget"] = state < 10**399
    gates["complete_heavy_stress_all_five_jets"] = local + state < 10**400
    gates["complete_normalized_heavy_stress"] = (local + state) / germs.KAPPA < R(
        1, 10**400
    )
    gates["complete_quantum_profile_nonaffine_four_jet"] = 6 * math.factorial(
        4
    ) * 64**4 * R(9, 55) ** 1024 < R(1, 8)
    gates["combined_finite_classical_and_heavy_tube_budget"] = R(3, 10**400) + R(
        1, 10**2700
    ) < R(1, 10**399)

    gates["closed_full_Riccati_ball_bracket"] = (
        400
        + s.Rational(1, 50) / loss**2
        + s.Rational(1, 10) / loss
        + s.Rational(1, 10000) / loss**2
    ) < 10**5
    gates["closed_auxiliary_Riccati_ball"] = s.Rational(4 * 10**5, 10**16) < s.Rational(
        1, 100
    )
    gates["full_log_rate_on_auxiliary_disk"] = 4 * s.Rational(1, 10**8) / loss < 1
    gates["full_physical_self_pair_amplitude"] = (
        2 * (s.Rational(3, 4 * 10**98) * (384 * 10**8 + 14**2) + s.Rational(14, 2))
        < 100
    )
    gates["complete_c2_exact_comparison_bound"] = 2 * 1000 * 20 * 10**65 < 10**70
    gates["all_complete_matrix_jet_bounds"] = (
        max(10**4 * math.factorial(j) * 4096**j for j in range(6)) < 10**30
    )
    gates["all_complete_defect_jet_bounds"] = (
        max(10**60 * math.factorial(j) * 4096**j for j in range(6)) < 10**85
    )
    gates["full_auxiliary_readout_bound"] = (
        s.Rational(3, 2) * (81 + 192 + s.Rational(196, 10**16)) < 10**4
    )
    gates["joint_time_Cauchy_tail_margin"] = (
        max(math.factorial(j) * 2048**j for j in range(6)) < 10**20
    )
    gates["complete_light_vacuum_profile_extra_bound"] = (
        s.Rational(3, 128 * 9 * germs.KAPPA)
        * 2
        * math.factorial(4)
        * 64**4
        * s.Rational(9, 55) ** 1024
    ) < s.Rational(1, 8 * 10**400)
    gates["full_heavy_profile_and_light_extra_fit_two_B"] = (
        s.Rational(9, 8) + s.Rational(1, 8) + s.Rational(1, 8) < 2
    )
    y, ss, rr = s.symbols("radial_y bump_s squared_squeezing_ratio", real=True)
    checks = {
        "full_all_momentum_radial_integral": s.integrate(
            y * y / (1 + y * y) ** s.Rational(5, 2), (y, 0, s.oo)
        )
        - s.Rational(1, 3),
        "sampling_complex_bump_denominator_margin": s.expand(
            17 * ss / 64 - ss / 4 - ss * ss / 64 - ss * (1 - ss) / 64
        ),
        "sampling_reciprocal_real_part_bound": s.Rational(2, 3) / s.Rational(4, 3) ** 2
        - s.Rational(3, 8),
        "sampling_time_radius_scaling": s.Rational(1, 8) / 8 - s.Rational(1, 64),
        "sampling_derivative_growth": 64 * s.Rational(4, 3) - s.Rational(256, 3),
        "full_uniform_smeared_energy_lower": s.Rational(1, 2 * 4**3)
        - s.Rational(1, 128),
        "complete_SLE_squeezing_magnitude_inequality": s.expand(
            (1 + 2 * rr) ** 2 * (1 - rr) - 1 - rr * (3 - 4 * rr * rr)
        ),
        "complete_actual_quantitative_stress_bound": (local + state) / germs.KAPPA
        - (local / germs.KAPPA + state / germs.KAPPA),
    }
    ad = quantum.adiabatic()
    P2, P4 = ad["P2"].subs(quantum.D, 3), ad["P4"].subs(quantum.D, 3)
    lam = -quantum.H * quantum.z
    U = 3 * quantum.Hd / 2 + 9 * quantum.H**2 / 4
    B2 = quantum.dt(P2) - 2 * lam * P2
    checks["full_Riccati_second_adiabatic_coefficient"] = s.expand(
        2 * P2 + U + quantum.dt(lam) / 2 - lam**2 / 4
    )
    checks["full_Riccati_fourth_adiabatic_coefficient"] = s.expand(
        P2**2 + 2 * P4 + quantum.dt(B2) / 2 - s.Rational(3, 2) * lam * B2
    )
    return {
        "actual_scope": "Every comoving p>=0 and real time[-1,1] at the fixed actual n,kappa0. The global state is defined separately. No numerical momentum cutoff or late-time global bound is introduced.",
        "complex_reference_domain": "Use the union of complex radius1/64 disks about[-1,1]. For v=1+t², Re v>99/100, |v|<33/16, |Im v|/Re v<1/30. The complete inverse-scale expression has Re(v^-4)>1/32, hence Re omega²>=Omega_star²/2 for Omega_star²=n+p²/16. The branch is holomorphic and |omega| lies betweenOmega_star/2 and8Omega_star. Full analytic H,lambda,U bounds are included; the smooth whole QG1 scalar is not assumed analytic.",
        "full_Riccati_iteration": "Normalize Wj=omega*wj and use the full square-root/logarithm map in state. Time domains lose1/1024 at each iteration. On |wj-1|<=1/100, Cauchy proves a complete map Lipschitz boundK|epsilon|²/Omega_star² with K=10^8. On the auxiliary epsilon disk |epsilon|<=Omega_star/10^8 this is at most10^-8. The entire iterates remain in the branch ball. At epsilon1 the jth increment is at most(K/Omega_star²)^j. Seven iterates give the complete W6 oscillator defect below10^60/Omega_star^12.",
        "full_exact_mode_comparison": "The exact normalized mode initialized to W6 at t=-1 differs by full variation-of-constants coefficients below10^62/Omega_star^13 over the entire interval. The exact exponent, not its first term, is bounded. Both comparison phases are retained and the evolution preserves CCR. The actual state remains the independently selected SLE.",
        "complete_sampling_derivative_bounds": dict(enumerate(weight)),
        "full_SLE_state_error": "Ten integrations by parts of the full physical self-pair amplitude, with all compact-sampling derivatives and exact-mode errors included, give |c2|<10^80/Omega_star^10. The Wronskian lower boundc1>=Omega_star/128 gives |beta|<=|c2|/c1<10^83/Omega_star^11. The complex coefficient recursion has sum bound100*8194^j/Omega_star^j; no phase is deleted before these integrations.",
        "full_matrix_and_readout_jet_bound": C,
        "exact_comparison_mode_derivative_coefficients": dict(enumerate(Y)),
        "complete_forced_error_derivative_coefficients": dict(enumerate(Z)),
        "complete_source_error_coefficients": dict(enumerate(F)),
        "complete_exact_mode_stress_difference_coefficients": dict(
            enumerate(exact_difference)
        ),
        "complete_SLE_stress_difference_coefficients": dict(enumerate(sle_difference)),
        "full_adiabatic_subtraction_tail_constant": cauchy_tail,
        "full_tail_proof": "The complete W6 stress readouts, with their phase already cancelled, are holomorphic in the auxiliary epsilon disk and bounded by10^4*Omega_star. Their first three even coefficients are the complete fourth-order adiabatic subtraction. The full tail at epsilon1, with a separate time Cauchy margin1/2048 through five derivatives, is below10^74/Omega_star^5. No complex Cauchy estimate is applied to the oscillatory mode phase.",
        "complete_mode_difference_constant": 10**301,
        "complete_state_dependent_and_subtraction_stress_upper": state,
        "complete_finite_local_heat_stress_upper": local,
        "complete_actual_heavy_stress_all_five_jets_upper": local + state,
        "complete_normalized_heavy_stress_all_five_jets_upper": (local + state)
        / germs.KAPPA,
        "full_radial_integral": "The complete sum is bounded by10^301/Omega_star^5 for every rho/P derivative through5. The full integral with p²dp/(2pi²) is64/(6pi² n)<2/n. A deliberately looser bound10^310/n is used. This proves absolute convergence and differentiation of the full subtracted integral in the stated interval.",
        "finite_local_proof": "Derive all finite terms from the whole dimensionally continued scalar counteraction, including its evanescent metric variation. The rational pole form of H bounds its full complex jets. Separate the constant vacuum term before differentiating: the local bound is[462n²+10^12(40000n+10^6)]/576<10^399. Adding the full state/subtraction term gives stress below10^400, or10^-400 against kappa0, for both components and all derivatives0..5.",
        "new_full_scalar_profile_four_jet_budget": "Use the fixed source-defined profile in quantum, including the extra mass1 light-vacuum constant. Split it into its affine clock part and the full1-T remainder. Only X is complexified for the smooth reference stress. On |u|<=1,7/8<=X<=9/8 the heavy affine four-jets are at most9B/8, B=10^-400. A complex X radius1/64 gives remainder belowB/8 and the light-vacuum term belowB/8; the whole new profile is below2B. Combining it with the S239 finite extensions and S238 classical changes stays below a NEW10^-399 normalized four-jet budget.",
        "not_inferred": "A small coefficient or reference-stress bound is not a bound for the full quantum response or its inverse. No old linear-stability rows are automatically transferred, no omitted interacting loop is controlled and no nonlinear same-state bounce, physical UV, finite-gravity Regge or original P8 closure is inferred.",
        "checks": {key: s.cancel(value) for key, value in checks.items()},
        "gates": {key: bool(value) for key, value in gates.items()},
    }
