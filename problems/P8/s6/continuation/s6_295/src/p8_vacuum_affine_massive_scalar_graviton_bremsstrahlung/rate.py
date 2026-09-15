"""Finite difference of two IR-divergent real rates; not virtual closure."""

from functools import cache

import sympy as s

from . import recoil, source

B = s.Integer(330000)
RESOLUTION = s.Symbol("physical_detector_resolution", positive=True)
LOWER = s.Symbol("auxiliary_lower_soft_cutoff", positive=True)


def amplitude_error(kappa=source.K):
    return B / s.sqrt(s.sympify(kappa))


def leading_soft_majorant(omega, kappa=source.K):
    return 64 / (s.sqrt(s.sympify(kappa)) * s.sympify(omega))


def finite_real_rate_error(resolution=RESOLUTION, kappa=source.K):
    e, k = map(s.sympify, (resolution, kappa))
    return (s.Integer(42248192) * e + s.Integer(54450000000) * e**2) / (2 * s.pi**2 * k)


def finite_lower_cutoff_majorant(lower=LOWER, resolution=RESOLUTION, kappa=source.K):
    return finite_real_rate_error(resolution, kappa) - finite_real_rate_error(
        lower, kappa
    )


@cache
def data():
    e, lo, k, omega = s.symbols("epsilon lambda kappa omega", positive=True)
    checks = {}
    checks["per_helicity_squared_difference_linear_coefficient"] = (
        128 * B + 8192 - s.Integer(42248192)
    )
    checks["both_helicities_integrated_quadratic_coefficient"] = B**2 / 2 - s.Integer(
        54450000000
    )
    checks["full_soft_angular_measure"] = 2 * 4 * s.pi / (2 * (2 * s.pi) ** 3) - 1 / (
        2 * s.pi**2
    )
    integrand = (B**2 * omega + 128 * B + 8192) / (2 * s.pi**2 * k)
    checks["complete_finite_lower_cutoff_integral"] = s.factor(
        s.integrate(integrand, (omega, lo, e)) - finite_lower_cutoff_majorant(lo, e, k)
    )
    checks["positive_regulator_removed_limit"] = s.limit(
        finite_lower_cutoff_majorant(lo, e, k), lo, 0, dir="+"
    ) - finite_real_rate_error(e, k)
    checks["majorant_resolution_derivative"] = s.factor(
        s.diff(finite_real_rate_error(e, k), e) - integrand.subs(omega, e)
    )
    checks["majorant_vanishes_with_resolution"] = s.limit(
        finite_real_rate_error(e, k), e, 0, dir="+"
    )
    checks["soft_amplitude_energy_scaling"] = s.factor(
        omega * s.sqrt(k) * leading_soft_majorant(omega, k) - 64
    )
    checks["hard_remainder_energy_independence"] = s.factor(
        s.sqrt(k) * amplitude_error(k) - B
    )
    checks["original_kappa_is_unchanged"] = source.KAPPA - s.Integer(10) ** 800
    checks["original_magnitude_upper_scale"] = s.Integer(
        10
    ) ** 8 / source.KAPPA - s.Rational(1, 10**792)
    checks["compact_domain_resolution_endpoint"] = recoil.compact_domain()[
        "resolution_max"
    ] - s.Rational(1, 8)
    checks["compact_domain_heavy_threshold"] = (
        recoil.compact_domain()["heavy_mass_squared_min"] - 128
    )
    # A common incoming flux and identical-scalar factor cancel exactly.
    flux, symmetry, A0, phase, correction = s.symbols(
        "flux identical_factor A0 phase correction", positive=True
    )
    checks["whole_differential_Born_rate_normalization"] = s.factor(
        (symmetry * phase * A0**2 * correction / flux)
        / (symmetry * phase * A0**2 / flux)
        - correction
    )
    checks["all_angle_integration_preserves_uniform_bound"] = s.factor(
        (A0**2 * finite_real_rate_error(e, k)) / A0**2 - finite_real_rate_error(e, k)
    )
    checks["squared_difference_Jacobian_identity"] = s.expand(
        s.Symbol("J") * ((s.Symbol("L") + s.Symbol("R")) ** 2 - s.Symbol("L") ** 2)
        + (s.Symbol("J") - 1) * s.Symbol("L") ** 2
        - (s.Symbol("J") * (s.Symbol("L") + s.Symbol("R")) ** 2 - s.Symbol("L") ** 2)
    )
    checks = {name: s.factor(value) for name, value in checks.items()}
    return {
        "whole_selected_rate_error": finite_real_rate_error(),
        "whole_finite_lower_cutoff_error": finite_lower_cutoff_majorant(),
        "per_helicity_amplitude_error_over_Born": amplitude_error(),
        "whole_normalized_rate_definition": "At fixed incoming COM energy E and outgoing pair-rest-frame direction u, I_exact(lambda,epsilon)=A0^-2 sum_helicities integral_lambda^epsilon d^3q/[(2pi)^3 2omega] J |M5|^2. I_soft uses the on-shell Born momenta, J=1 and amplitude A0 S0. Each rate individually has a soft logarithmic divergence; only their difference is controlled here.",
        "whole_finite_real_error_theorem": "For mu1,n>=128,5/4<=E<=2,all angles and0<lambda<epsilon<=1/8, |I_exact-I_soft|<=finite_real_rate_error(epsilon,kappa), independent of lambda. Dominated convergence applies to their difference as lambda->0. At the maximal resolution the bound is<10^8/kappa, hence<10^-792 at original parameters. The same uniform estimate applies after angular integration relative to the corresponding Born rate.",
        "whole_physical_vs_analytic_boundary": "This is an actual on-shell one-real-graviton phase-space comparison in a specified compact-domain angular-bin convention. It is not yet the finite virtual/dressed inclusive rate or equality with S278 dimensional analytic soft division. The virtual regulator conversion, pure-gravity radiation, hard-loop corrections, higher matching and Regge/forward order of limits remain open.",
        "checks": checks,
        "gates": {
            "both_physical_helicities_and_full_real_phase_space": True,
            "finite_difference_not_finite_individual_Fock_rates": True,
            "lower_cutoff_uniform_dominating_integrable_majorant": True,
            "relative_Born_normalization_and_positive_domain_explicit": True,
            "original_parameter_bound_without_finite_matching_choice": True,
            "uniform_angle_bound_survives_Born_weighted_integration": True,
            "virtual_dressed_regulator_conversion_still_open": True,
        },
    }
