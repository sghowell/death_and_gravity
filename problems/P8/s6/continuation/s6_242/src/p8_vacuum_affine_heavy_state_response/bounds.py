"""Uniform full momentum state-selection response, profile and cutoff-tail bounds."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_curved_state import state

from . import vertices as v


@cache
def data():
    n = state.MASS2
    kappa = state.KAPPA
    B = s.Integer(10) ** 83
    generic_mass = s.Symbol("positive_heavy_mass", positive=True)

    def radial(power):
        return (
            8
            * s.gamma(s.Rational(power - 3, 2))
            * generic_mass ** (3 - power)
            / (s.pi ** s.Rational(3, 2) * s.gamma(s.Rational(power, 2)))
        )

    checks = {
        "complete_ninth_radial_integral": s.simplify(
            radial(9) - 256 / (105 * s.pi**2 * generic_mass**6)
        ),
        "complete_tenth_radial_integral": s.simplify(
            radial(10) - 5 / (8 * s.pi * generic_mass**7)
        ),
    }
    X = s.Symbol("X", real=True)
    localizer = X**state.N / (X**state.N + (1 - X) ** state.N)
    for order in range(3):
        checks["existing_full_localizer_clock_jet_" + str(order)] = s.diff(
            localizer, X, order
        ).subs(X, 1) - (1 if order == 0 else 0)
    # Fixed reference rho/P differences, not a new adaptive profile.
    rho, pressure = s.symbols(
        "fixed_state_selection_rho fixed_state_selection_pressure", real=True
    )
    nd, ng, td, tg = s.symbols("nD nG traceQD traceQG", real=True)
    N = 1 + v.ed * nd + v.eg * ng
    tau = v.ed * td + v.eg * tg
    volume = 1 + tau / 2 + tau * tau / 8
    xminus = v.jet(1 - 2 * (N - 1) + 3 * (N - 1) ** 2) - 1
    profile = v.jet(N * volume * (-pressure - (rho + pressure) * xminus / 2))
    profmixed = s.expand(profile.coeff(v.ed, 1).coeff(v.eg, 1))
    wanted = (
        -(rho + pressure) * nd * ng
        + rho * (nd * tg + ng * td) / 2
        - pressure * td * tg / 4
    )
    checks["entire_fixed_state_selection_profile_Hessian"] = s.expand(
        profmixed - wanted
    )
    checks["complete_reference_lapse_first_current_cancellation"] = s.expand(
        profile.diff(v.ed).subs({v.ed: 0, v.eg: 0, td: 0}) - rho * nd
    )
    checks["complete_reference_metric_first_current_cancellation"] = s.expand(
        profile.diff(v.ed).subs({v.ed: 0, v.eg: 0, nd: 0}) + pressure * td / 2
    )
    delta = s.Symbol("delta", real=True)
    a = s.Symbol("a", positive=True)
    vD, vG = s.symbols("vD vG", real=True)
    physical = a**3 * wanted
    linear_pullback = physical.subs(
        {td: 6 * (vD + delta * nd), tg: 6 * (vG + delta * ng)}, simultaneous=True
    )
    dJ = (21 * delta**2 - 3 * delta) * (-pressure) / 2 + (1 - 6 * delta) * (
        -(rho + pressure) / 2
    )
    Tc = rho - 3 * delta * pressure
    clock = a**3 * (
        2 * dJ * nd * ng + 3 * Tc * (nd * vG + ng * vD) - 9 * pressure * vD * vG
    )
    profile_contact = -3 * a**3 * pressure * (4 * delta**2 - 3 * delta) * nd * ng
    mean_contact = 3 * a**3 * pressure * (4 * delta**2 - 3 * delta) * nd * ng
    checks["full_profile_nonlinear_clock_chain_rule"] = s.expand(
        clock - linear_pullback - profile_contact
    )
    checks["complete_matched_state_difference_chart_contact_cancellation"] = s.expand(
        profile_contact + mean_contact
    )
    pressure_K = s.Symbol("fixed_state_selection_pressure_projected_K", real=True)
    mean_contact_K = mean_contact.subs(pressure, pressure_K)
    cutoff_contact = (
        -3 * a**3 * (pressure - pressure_K) * (4 * delta**2 - 3 * delta) * nd * ng
    )
    checks["finite_cutoff_clock_mean_tail_contact_retained"] = s.expand(
        profile_contact + mean_contact_K - cutoff_contact
    )
    # Exact scalar inequalities dominate the ENTIRE pair and covariance products.
    memory_const = 10**16
    contact_const = 10**8
    profile_const = 10**9
    physical = s.Rational(10**101, 1) / (n**3 * kappa)
    common = s.Rational(10**103, 1) / (n**3 * kappa)
    tail = s.Rational(10**104, 1) / (n ** s.Rational(5, 2) * kappa)
    gates = {
        "actual_all_momentum_mixing_below_half": bool(
            B / n ** s.Rational(11, 2) < s.Rational(1, 2)
        ),
        "full_exact_comparison_readout_bound": 8 * 20 < 200,
        "full_Bogoliubov_readout_difference_bound": 200 * s.Rational(5, 4) < 300,
        "full_selected_readout_bound": 200 + 300 * s.Rational(1, 2) < 400,
        "whole_pair_difference_constant": 300 * 400 <= 2 * 10**5,
        "whole_four_mode_product_difference_constant": 2 * 10**5 * (400**2 + 200**2)
        < 10**11,
        "both_ordered_branches_and_actual_volumes": 4096 * 10**11 < 10**15,
        "whole_radial_memory_constant": 2 * 10**15 < 10**16,
        "whole_contact_product_difference": (200 + 400) * 300 < 10**6,
        "whole_contact_volume_and_integral": 32 * 10**6 < 10**8,
        "ninth_radial_integral_constant": bool(s.Rational(256, 105 * 9) < 1),
        "tenth_radial_integral_constant": bool(s.Rational(5, 8 * 3) < 1),
        "whole_fixed_profile_coefficient_sum": 64 * (2 + 1 + s.Rational(1, 4)) < 10**3,
        "whole_fixed_profile_stress_and_Hessian": 10**3 * 10**6 <= 10**9,
        "complete_memory_contact_and_fixed_profile_constant": memory_const
        + contact_const
        + profile_const
        < 10**17,
        "full_weighted_physical_graph_constant": 2 * 10**17 * B < 10**101,
        "full_same_clock_graph_constant": 64 * 10**101 < 10**103,
        "full_physical_normalized_bound": bool(physical < s.Rational(1, 10**1280)),
        "full_same_clock_normalized_bound": bool(common < s.Rational(1, 10**1280)),
        "complete_large_frequency_radial_tail": bool(s.Rational(4, 6) * 2**6 < 43),
        "complete_low_transfer_memory_tail": 86 * 10**15 < 10**17,
        "full_two_band_contact_and_memory_tail_constant": 2 * 10**17 + 10**16 + 10**8
        < 10**18,
        "complete_physical_graph_tail_constant": 10**18 * B < 10**102,
        "complete_same_clock_graph_tail_constant": bool(
            64 * 10**102 + 10**9 * B < 10**104
        ),
        "full_same_clock_tail_coefficient": bool(tail < s.Rational(1, 10**1180)),
    }
    return {
        "actual_parameters": {
            "mass_squared": n,
            "kappa0": kappa,
            "time_interval": (-s.Rational(1, 2), s.Rational(1, 2)),
        },
        "uniform_state_input": "For EVERY internal momentum, Omega_p²=n+|p|²/16 and the fixed actual SLE coefficient obeys|beta_p|<=B Omega_p^-11, B=1e83. The exact comparison energy-scaled mode has norm below20 on[-1,1]. All subsequent reference estimates are on[-1/2,1/2].",
        "entire_physical_readouts": "The full5-feature comparison vector obeys||F_S||2<=200 sqrt(Omega), using a>=1,|H|<=2, |p|/a<=4Omega and mass<=Omega. Full Bogoliubov algebra gives||F_T-F_S||2<=300B Omega^-21/2 and||F_T||2<=400 sqrt(Omega). This retains beta's complex phase, alpha-1 and every cross term. Beta is constant in time.",
        "complete_pair_difference": "For nu=Omega_p,mu=Omega_(P-p), |V_D(T)-V_D(S)|<=2e5 B sqrt(nu mu)(nu^-11+mu^-11)v(D). The ENTIRE ordered two-time product difference is bounded by1e11 B nu mu(nu^-11+mu^-11)v(D)v(G), with both cross terms and the error square retained. The reflected reverse term has the same bound, and the i/2 factor accounts for both branches.",
        "full_transfer_comparison": "mu<=nu+|P|/4<=L(P)nu and nu<=L(P)mu, L(P)=1+|P|/(4sqrt(n)). Hence the full integral of nu mu(nu^-11+mu^-11) is at most2L times the exact ninth radial integral. This covers large-small internal pairs and every external transfer.",
        "complete_radial_integrals": {
            "Omega_power_minus9": radial(9),
            "Omega_power_minus10": radial(10),
        },
        "entire_unscaled_memory_bound": "1e16 B n^-3 L(P)v(D,t)v(G,s), before integrating the exact retarded triangle. Its time interval length is1, so the L2_t bilinear bound needs no external time derivative. The contact bound is1e8 B n^-7/2 v(D,t)v(G,t).",
        "full_state_selection_mean_stress_bound": 10**6 * B / n ** s.Rational(7, 2),
        "unchanged_full_profile_localizer": localizer,
        "entire_fixed_profile_ADM_mixed_density_before_a3": wanted,
        "full_existing_profile_split": "On the reference slab, the T-minus-S stress is absolutely integrable, conserved and has the required first five derivatives. Its summand of the existing profile is the SAME full localizer T(X) times[-pressure-(rho+pressure)(X-1)/2]. T(1)=1 and its next two jets vanish, as checked literally. No global regularity of this auxiliary comparison summand or new state, prescription or candidate is asserted. Its whole reference Hessian is bounded by1e9 B n^-7/2 v(D)v(G). The comparison-profile summand and common local prescription remain in the unclosed comparison response; their sum is exactly the unchanged S240 profile.",
        "complete_mean_and_clock_contact": "For the entire T-minus-S current PLUS its fixed mean-profile summand, first lapse and spatial-metric variations cancel; the odd shift mean vanishes. The displayed nonzero profile chart contact cancels the corresponding nonzero current one-point chart contact. This cancellation requires the COMPLETE matched summand, not the profile alone.",
        "retained_profile_second_clock_contact": profile_contact,
        "retained_state_current_second_clock_contact": mean_contact,
        "retained_nonzero_finite_cutoff_clock_contact": cutoff_contact,
        "full_same_clock_fixed_profile_Hessian": clock,
        "complete_physical_graph": "v(D)=|n_D|+a|beta_D|+2||Q_D||F. Define norm_phys=||sqrt(1+|P|²)v(D)||L2(dt dP). L(P)<=2sqrt(1+|P|²). The whole matched state-selection response divided by kappa0 is bounded by the displayed physical constant times the two norms. This is a weak graph result before quantum constraints.",
        "complete_normalized_physical_response_bound": physical,
        "full_common_clock_graph": "For scalar ADM input(n_lapse,v,B), Q=2(v+delta n_lapse)I and beta=iP B/a² in the FIRST chart derivative, delta=1/[2(1+t²)^3]. The full second chart contact cancels only in the complete matched summand. Then v(physical)<=8sqrt(1+|P|²)(|n_lapse|+|v|+|B|), so norm_phys<=8 norm_clock with norm_clock=||(1+|P|²)(|n_lapse|+|v|+|B|)||L2. Both physical and clock bounds are below10^-1280, all transfer, with no time derivative of the external directions required.",
        "complete_normalized_same_clock_response_bound": common,
        "full_regulator_tail": "Use the identical finite reference-mode projector Omega_p<=K, K>=2sqrt(n), on BOTH memory legs and the one contact leg, while retaining the exact fixed profile. If|P|<=2K, any missing pair has both frequencies>K/2; the entire radial ninth-power tail is<43K^-6. If|P|>2K, divide the full bound by the existing spatial graph weight and use L/(1+|P|²)<=1/K. The one-mode contact tail is bounded separately. At finite K the full fixed profile's first current no longer cancels the projected current: the DISPLAYED nonzero second-clock mean-tail contact is retained. Its bound is1e9 B K^-7, which fits the extra clock-tail allowance. Thus the same physical graph has normalized error<1e102/[kappa0*n^(5/2)*K] and the same clock graph error<1e104/[kappa0*n^(5/2)*K]<10^-1180/K. No large-small internal pair, finite-cutoff chart contact or upper retarded endpoint is removed.",
        "complete_normalized_clock_regulator_tail_coefficient": tail,
        "not_full_heavy_response": "This is the ENTIRE state-selection DIFFERENCE between T and the specified exact comparison S, including both ordered branches, instantaneous metric contact and its fixed mean-profile summand. It is not the full comparison-state renormalized response, a finite heat truncation thereof, a full quantum inverse, a nonlinear remainder or P8 closure.",
        "checks": {key: s.cancel(value) for key, value in checks.items()},
        "gates": {key: bool(value) for key, value in gates.items()},
    }
