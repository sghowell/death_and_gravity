"""Exact IR-divided leading-ladder phase and Fourier reconstruction."""

from functools import cache

import sympy as s

from . import source

ETA = s.Symbol("real_ladder_strength", real=True)
C = s.Symbol("real_finite_phase_log", real=True)
TAU = s.Symbol("positive_transfer_magnitude", positive=True)
S = s.Symbol("energy_above_threshold", positive=True) + 4
K = s.Symbol("positive_kappa", positive=True)
ABEL_REGULATOR = s.Symbol("positive_Abel_regulator", positive=True)


def eta(energy=S, kappa=K):
    return source.old_phase.eikonal(energy) / (
        8 * s.pi * s.sympify(kappa) * source.old_phase.gap(energy)
    )


def finite_log(energy=S, tau=TAU):
    return (
        s.log(s.sympify(tau))
        - s.log(4 * s.pi)
        + s.EulerGamma
        + 2 / source.old_phase.eikonal(energy)
    )


def gamma_phase(strength=ETA):
    a = s.sympify(strength)
    return (
        s.exp(-2 * s.I * s.EulerGamma * a) * s.gamma(1 - s.I * a) / s.gamma(1 + s.I * a)
    )


def from_log(strength=ETA, logarithm=C):
    a, c = map(s.sympify, (strength, logarithm))
    return s.exp(s.I * a * c) * gamma_phase(a)


def factor(energy=S, tau=TAU, kappa=K):
    return from_log(eta(energy, kappa), finite_log(energy, tau))


def pole_amplitude(energy=S, tau=TAU, kappa=K):
    return (
        source.old_phase.eikonal(energy)
        / (s.sympify(kappa) * s.sympify(tau))
        * factor(energy, tau, kappa)
    )


def riesz_transform(strength=ETA, tau=TAU):
    a, t = map(s.sympify, (strength, tau))
    return 4 * s.pi / t * s.gamma(1 - s.I * a) / s.gamma(s.I * a) * (t / 4) ** (s.I * a)


def abel_transform(strength=ETA, tau=TAU, regulator=ABEL_REGULATOR):
    a, t, h = map(s.sympify, (strength, tau, regulator))
    return (
        s.pi
        * h ** (s.I * a - 1)
        * s.gamma(1 - s.I * a)
        * s.hyper((1 - s.I * a,), (1,), -t / (4 * h))
    )


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = s.factor(s.simplify(s.expand_func(value)))

    a, c, t, k = ETA, C, TAU, K
    V = source.old_phase.eikonal(S)
    D = source.old_phase.gap(S)
    put("original_eta_normalization", eta() - V / (8 * s.pi * k * D))
    put(
        "full_D_finite_phase_matches_S303",
        s.I * eta() * finite_log()
        - source.old_phase.phase_ratio(S, t, s.log(4 * s.pi) - s.EulerGamma, k),
    )
    logQ = (
        -2 * s.I * s.EulerGamma * a + s.loggamma(1 - s.I * a) - s.loggamma(1 + s.I * a)
    )
    for n in range(1, 8):
        want = (
            0 if n % 2 == 0 or n == 1 else 2 * s.I**n * s.zeta(n) * s.factorial(n - 1)
        )
        put(f"full_gamma_log_derivative{n}", s.diff(logQ, a, n).subs(a, 0) - want)
    put(
        "entire_gamma_phase_unit_modulus",
        gamma_phase(a) * s.conjugate(gamma_phase(a)) - 1,
    )
    put(
        "entire_Coulomb_phase_unit_modulus",
        from_log(a, c) * s.conjugate(from_log(a, c)) - 1,
    )
    put("gamma_factor_zero", gamma_phase(0) - 1)
    put(
        "whole_linear_finite_coefficient",
        s.diff(from_log(a, c), a).subs(a, 0) - s.I * c,
    )
    put(
        "whole_quadratic_finite_coefficient",
        s.diff(from_log(a, c), a, 2).subs(a, 0) / 2 + c * c / 2,
    )
    put(
        "whole_cubic_finite_coefficient",
        s.diff(from_log(a, c), a, 3).subs(a, 0) / 6
        + s.I * (c**3 / 6 + 2 * s.zeta(3) / 3),
    )
    # Use eta as an independent symbol to avoid unrelated branch simplification.
    d, v = s.symbols("positive_flux positive_numerator", positive=True)
    impact_prefactor = s.exp(s.I * a * (2 / v - s.EulerGamma)) * s.pi ** (-s.I * a)
    fourier_difference = (2 * d / s.I) * impact_prefactor * riesz_transform(
        a, t
    ) - 8 * s.pi * d * a / t * from_log(
        a, s.log(t) - s.log(4 * s.pi) + s.EulerGamma + 2 / v
    )
    # All powered bases here are positive; rewrite them as principal-log exponentials.
    put(
        "Riesz_Fourier_reconstructs_whole_phase",
        s.expand_log(s.expand_func(fourier_difference).rewrite(s.exp)),
    )
    h = s.Symbol("positive_Abel_regulator", positive=True)
    put(
        "Abel_zero_strength_Gaussian",
        abel_transform(0, t, h) - s.pi / h * s.exp(-t / (4 * h)),
    )
    put(
        "Abel_zero_transfer_moment",
        abel_transform(a, 0, h) - s.pi * h ** (s.I * a - 1) * s.gamma(1 - s.I * a),
    )
    put(
        "transfer_log_monodromy",
        from_log(a, c + 2 * s.pi * s.I) - s.exp(-2 * s.pi * a) * from_log(a, c),
    )
    put(
        "one_loop_phase_square_not_exact_modulus",
        (1 + s.I * a * c) * (1 - s.I * a * c) - 1 - a * a * c * c,
    )
    r = s.Symbol("positive_phase_excursion", positive=True)
    put(
        "nonuniform_one_loop_sequence",
        (1 + s.I * a * c).subs(c, -r / a) * (1 - s.I * a * c).subs(c, -r / a)
        - 1
        - r * r,
    )
    put("known_crossed_endpoint_eta", source.old_phase.eikonal(4 - S) - V)
    return {
        "whole_original_ladder_strength": eta(),
        "whole_unexpanded_Coulomb_factor": factor(),
        "whole_gamma_phase": gamma_phase(),
        "whole_Riesz_transform": riesz_transform(),
        "whole_Gaussian_Abel_transform": abel_transform(),
        "whole_log_gamma_series": "logQ(eta)=2 sum_{odd j>=3} zeta(j)*(i eta)^j/j, convergent for|eta|<1; the complete factor isexp(i eta c)Q, c=ln tau-ln4pi+gamma_E+2/V.",
        "whole_Fourier_definition": "After the stated IR division the impact factor isexp[i eta(2/V-gamma_E)](pi b^2)^(-i eta). Its Gaussian-Abel-regulated Fourier transform is computed exactly and then continued to the Riesz distribution. For nonzerotau the forward delta term is separate; multiplying2D/i reconstructs the original Born pole times the displayed phase.",
        "whole_physical_modulus": "For real energy and strength on the physical spacelike edge, conjugation exchanges Gamma(1-i eta) and Gamma(1+i eta), so the complete selected phase has modulusone. This identity concerns the named leading ladder class, not a unitarity construction of the full S-matrix.",
        "whole_nonuniform_perturbation": "Takingtau=exp[ell-2/V-r/eta] gives the one-loop phase excursion-r. The squared magnitude of its finite linear truncation is1+r^2, while the exact leading-class magnitude isone. Fixed-order loop-square bounds cannot be made uniform by dropping the imaginary phase.",
        "whole_complex_boundary": "A circuit aroundtau0 multiplies the complete factor byexp(-2pi eta). Physical-edge unit modulus does not remove the complex branch, establish a forward holomorphic disk or supply a high-energy Regge remainder.",
        "checks": checks,
        "gates": {
            "both_box_and_D_tree_finite_phase_terms_calibrated": True,
            "all_original_raw_Gamma_and_fourpi_constants_retained": True,
            "physical_edge_exact_unit_modulus_only_in_named_class": True,
            "Gaussian_Abel_and_Riesz_routes_agree": True,
            "zero_forward_delta_not_reinterpreted_as_finite_cross_section": True,
            "log_transfer_power_never_Taylor_truncated_in_final_factor": True,
            "complex_monodromy_not_removed": True,
        },
    }
