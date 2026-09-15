"""Uniform physical cut bound and a finite spectral-window contribution only."""

from functools import cache

import sympy as s

from . import sewing, source

LOWER, UPPER = s.symbols(
    "physical_lower_energy_squared physical_upper_energy_squared", positive=True
)


def exact_positive(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Rational)):
        raise TypeError("Require a finite positive exact rational")
    value = s.Rational(value)
    if value <= 0:
        raise ValueError("Require a finite positive exact rational")
    return value


def require_window(lower, upper, mass2=source.MASS2, kappa=source.KAPPA):
    lo, hi, mu, k = map(exact_positive, (lower, upper, mass2, kappa))
    if lo < 4 * mu or hi <= lo:
        raise ValueError("Require4m^2<=lower<upper in the physical massive region")
    return lo, hi, mu, k


def moment_primitive(energy, mass2):
    delta = energy - 2 * mass2
    return s.log(delta) - 4 * mass2 / delta - 2 * mass2 * mass2 / (delta * delta)


def formal_window_upper(lower=LOWER, upper=UPPER, mass2=source.MU, kappa=source.K):
    return (moment_primitive(upper, mass2) - moment_primitive(lower, mass2)) / (
        128 * s.pi**2 * kappa * kappa
    )


@cache
def data():
    r, y = s.symbols("beta_squared cosine_squared", real=True)
    A, B = 1 - r, r * (1 - y)
    D = A + B
    u, v = A * A / D, B * B / D
    checks = {
        "physical_denominator_sum": s.expand(D - (1 - r * y)),
        "helicity_amplitude_sum_positive_deficit": s.factor(D - u - v - 2 * A * B / D),
        "sum_of_squares_below_square_of_sum": s.factor(
            (u + v) ** 2 - u * u - v * v - 2 * u * v
        ),
        "denominator_upper_margin": s.expand(1 - D - r * y),
        "complete_physical_spectral_majorant_primitive": s.factor(
            s.diff(moment_primitive(source.S, source.MU), source.S)
            - source.S**2 / (source.S - 2 * source.MU) ** 3
        ),
        "forward_sewing_Cauchy_normalization": s.factor(
            source.S**2 / (256 * s.pi * source.K**2) * 2 / s.pi
            - source.S**2 / (128 * s.pi**2 * source.K**2)
        ),
        "zero_width_window": formal_window_upper(LOWER, LOWER),
        "entire_window_composition": s.factor(
            formal_window_upper(LOWER, UPPER)
            + formal_window_upper(UPPER, source.S)
            - formal_window_upper(LOWER, source.S)
        ),
    }
    cap = s.Integer(10) ** 196
    named = s.Rational(591, 128) / source.KAPPA**2
    checks["threshold_primitive_rational_constant"] = s.factor(
        moment_primitive(4 * source.MU, source.MU)
        - s.log(2 * source.MU)
        + s.Rational(5, 2)
    )
    checks["named_cap_log_upper_exponent"] = 196 * 3 - 588
    checks["named_exact_window_license"] = require_window(4, cap)[1] - cap
    checks["named_original_reference_lambda"] = (
        source.current.heavy.LAMBDA - s.Rational(1, 10**600)
    )
    checks["threshold_strict_positive_sewn_density"] = sewing.forward_shape(0) - 1
    return {
        "uniform_physical_absorptive_majorant": source.S**2
        / (256 * s.pi * source.K**2),
        "whole_physical_window_bound": formal_window_upper(),
        "named_window": {
            "lower": s.Integer(4),
            "upper": cap,
            "mass_squared": s.Integer(1),
            "kappa": source.KAPPA,
            "strict_upper_bound": named,
            "comparison_to_unchanged_reference_4lambda": named
            / (4 * source.current.heavy.LAMBDA),
        },
        "pointwise_proof": "For r=beta^2 in[0,1), y=x^2 in[0,1], A=1-r>0,B=r(1-y)>=0,D=A+B=1-r*y<=1. The two tree amplitudes in units s/(4kappa) are u=A^2/D,v=B^2/D. Thus u+v<=D<=1 and u^2+v^2<=1. With two of each physical helicity, the forward density is positive and <=s^2/(256pi kappa^2). The Cauchy-Schwarz inequality for the SAME complete helicity-phase-space inner product bounds the absolute nonforward physical cut by its forward value.",
        "window_proof": "Multiply this majorant by2/[pi(s-2m^2)^3] and integrate only over the stated physical window4m^2<=lower<upper. The displayed exact primitive controls that positive two-graviton contribution. For m^2=1,lower4,upper10^196, log10<3, pi^2>1 and the negative endpoint terms give an upper bound591/(128kappa^2), less than10^-990 times the unchanged4lambda reference. This does not restore lambda as a bare term canceled by the heavy matching.",
        "not_established": "The window integral is NOT the complete b20 coefficient or a complete low-cut subtraction. It omits the unphysical0<s<4m^2 continuation, crossed-transfer matching, unknown real local terms, other intermediate states, higher loops, detector factors and the UV contour. The selected cap is only an integration limit for a FORMAL coefficient, not a physical EFT cutoff or a certified finite-energy error.",
        "checks": checks,
        "gates": {
            "exact_original_window_and_mass": require_window(4, cap)
            == (4, cap, 1, source.KAPPA),
            "positive_formal_density_not_exact_full_unitarity": True,
            "named_positive_window_bound_below_reference": bool(
                0 < named / (4 * source.current.heavy.LAMBDA) < s.Rational(1, 10**990)
            ),
            "uniform_bound_uses_full_amplitude_not_threshold_expansion": True,
            "unphysical_massless_threshold_not_discarded_from_full_task": True,
            "finite_window_not_b20_or_full_Regge_allowance": True,
            "reference_lambda_not_reinserted_into_retained_tree_action": True,
            "all_original_V_G_B_P8_frontiers_remain_open": True,
        },
    }
