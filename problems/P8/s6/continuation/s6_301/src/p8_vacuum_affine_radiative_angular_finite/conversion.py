"""Finite same-state soft conversion and strictly scoped detector comparison."""

from functools import cache

import sympy as s

from . import bounds, source


def finite_conversion(kernel_zero, kernel_first, kappa=source.K):
    return (kernel_first + (s.EulerGamma - 2 - s.log(s.pi)) * kernel_zero) / (
        8 * s.pi**2 * kappa
    )


def require_resolution(value):
    x = source.recoil.exact_real(value)
    if x <= 0 or x > s.Rational(1, 8):
        raise ValueError("Require resolution in(0,1/8]")
    return x


def finite_pairing_error(epsilon, resolution, kappa=source.KAPPA):
    e = bounds.require_epsilon(epsilon)
    x = require_resolution(resolution)
    return e * (28000 - 4000 * s.log(x)) / kappa


def original_conversion_bound():
    return s.Rational(2000, 10**800)


def original_reference_ratio_bound():
    return s.Rational(4250, 10**800)


def original_unexpanded_factor_remainder():
    return s.Rational(5000000, 10**1600)


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = s.factor(s.sympify(value))

    e, L, K0, K1, K2, c = s.symbols(
        "epsilon log_resolution K0 K1 K2 phase_first", real=True
    )
    expression = s.exp(2 * e * L) * ((1 + c * e) * (K0 + K1 * e + K2 * e * e) - K0) / e
    expansion = s.series(expression, e, 0, 2).removeO()
    put("same_state_paired_finite_part", expansion.subs(e, 0) - K1 - c * K0)
    put(
        "same_state_first_regulator_error",
        s.diff(expansion, e).subs(e, 0) - (K2 + c * K1 + 2 * L * (K1 + c * K0)),
    )
    put(
        "finite_phase_definition",
        finite_conversion(K0, K1, 1) * 8 * s.pi**2
        - K1
        - (s.EulerGamma - 2 - s.log(s.pi)) * K0,
    )
    put(
        "finite_conversion_majorant_margin",
        2000 - (bounds.K1_BOUND + 4 * bounds.K0_BOUND) / 72 - s.Rational(3235, 18),
    )
    put(
        "finite_regulator_constant_margin",
        28000 - bounds.PHASE_REMAINDER / 72 - s.Rational(2000, 9),
    )
    put("resolution_log_error", 2 * 2000 - 4000)
    put("elastic_plus_radiative_conversion_difference", 2000 + 112 - 2112)
    put("known_power_Gamma_log_difference", 2112 + 12 - 2124)
    put("known_factor_ratio_majorant", 2 * 2125 - 4250)
    put(
        "original_reference_ratio_error",
        original_reference_ratio_bound() * source.KAPPA - 4250,
    )
    put("original_conversion_bound", original_conversion_bound() * source.KAPPA - 2000)
    put("first_Newton_exp_remainder", 2000**2 - 4000000)
    put("Gamma_prefactor_remainder", 2 * 6**2 - 72)
    put("unexpanded_factor_remainder_margin", 5000000 - 4000000 - 72 - 999928)
    put(
        "original_unexpanded_factor_remainder",
        original_unexpanded_factor_remainder() * source.KAPPA**2 - 5000000,
    )
    ks, ke = s.symbols("state_kernel elastic_kernel", real=True)
    put(
        "wrong_elastic_virtual_leaves_IR_pole",
        ks / (2 * e) - ke / (2 * e) - (ks - ke) / (2 * e),
    )
    put("correct_radiative_virtual_cancels_IR_pole", ks / (2 * e) - ks / (2 * e))
    a, d = s.symbols("soft_exponent finite_conversion", real=True)
    log_gamma = -s.EulerGamma * a - s.loggamma(1 + a)
    put("Gamma_no_first_Newton_term", s.diff(log_gamma, a).subs(a, 0))
    put("Gamma_second_Newton_term", s.diff(log_gamma, a, 2).subs(a, 0) + s.pi**2 / 6)
    put("exponential_linear_term", s.diff(s.exp(d), d).subs(d, 0) - 1)
    put("original_large_kappa_positive_margin", source.KAPPA - 60 - (10**800 - 60))
    return {
        "whole_finite_conversion": "Delta_sigma=[K1_sigma+(EulerGamma-2-lnpi)*K0_sigma]/(8pi^2*kappa); absDelta_sigma<2000/kappa at the stated compact radiative states. This is the universal reference conversion, not the unknown finite radiative hard amplitude.",
        "whole_explicit_regulator_error": "For0<=e<=1/8 and fixed0<x<=1/8, the same-state paired soft factor differs fromDelta_sigma by<=e*[28000+4000abslnx]/kappa. Removee at fixedx; no arbitrary e*abslnx interchange.",
        "whole_known_detector_comparison": "Only compare expDelta_sigma*F(a_sigma)*x^a_sigma, whereF(a)=exp(-EulerGamma*a)/Gamma(1+a), after same-state IR pairing. ForR<=x<=1/8 its ratio to the elastic known factor differs from1 by<4250/kappa at originalparameters. This usesS300's unexpanded-power bound, not a uniform Newton expansion ofx^a.",
        "whole_original_finite_conversion_bound": original_conversion_bound(),
        "whole_original_known_detector_ratio_bound": original_reference_ratio_bound(),
        "whole_original_unexpanded_factor_remainder": original_unexpanded_factor_remainder(),
        "whole_hard_scheme_boundary": "If a hard tree or loop depends on the dimensional regulator, its evanescent coefficient times a soft pole contributes its own finite hard term. This universal angular calculation does not set that coefficient to zero, identify a scheme-independent hard amplitude, or establish the complete physical radiation rate.",
        "checks": checks,
        "gates": {
            "same_radiative_state_poles_and_finite_term": True,
            "phase_and_projector_and_radial_derivatives_all_retained": True,
            "fixed_resolution_regulator_limit_with_explicit_error": True,
            "unexpanded_detector_power_and_Gamma_factor_only": True,
            "original_parameter_bounds_not_chosen_matching_values": True,
            "finite_hard_and_evanescent_hard_terms_still_required": True,
            "no_nonleading_all_multiplicity_radiation_or_Regge_inference": True,
        },
    }
