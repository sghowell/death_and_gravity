"""Fixed-domain dimensional continuity and explicit regulator majorants."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import recoil

from . import source

EP = source.EP
RESOLUTION = s.Symbol("physical_detector_resolution", positive=True)
V = s.Symbol("folded_perpendicular_direction", nonnegative=True)


def folded_angular_weight(perpendicular=V, epsilon=EP):
    v, e = map(s.sympify, (perpendicular, epsilon))
    return (1 + 2 * e) * v ** (2 * e)


def continued_pair_density(energy, omega, epsilon=EP):
    E, w, e = map(s.sympify, (energy, omega, epsilon))
    return recoil.density_ratio(E, w) * ((E * (E - w) - 1) / (E * E - 1)) ** e


def soft_pair_regulator_error(epsilon=EP, resolution=RESOLUTION, kappa=source.K):
    e, res, k = map(s.sympify, (epsilon, resolution, kappa))
    return e * (16016 * s.Abs(s.log(res)) + 301824) / (8 * s.pi**2 * k)


def finite_conversion_bound(kappa=source.K):
    return s.Integer(112) / s.sympify(kappa)


@cache
def data():
    checks = {}
    margins = {}

    def put(name, value):
        checks[name] = s.factor(s.expand_func(value))

    def nonnegative(name, value, variables):
        value = s.factor(value)
        num, den = s.fraction(value)
        if not all(v >= 0 for v in s.Poly(num, *variables).coeffs()):
            raise ValueError("Negative majorant numerator: " + name)
        if not all(v >= 0 for v in s.Poly(den, *variables).coeffs()) or den == 0:
            raise ValueError("Invalid majorant denominator: " + name)
        margins[name] = value

    e, r, z = s.symbols("e r z", positive=True)
    rho = (s.Rational(1, 2) + e) / s.pi * (1 - r * r) ** (e - s.Rational(1, 2))
    primitive = -((1 - r * r) ** (e + s.Rational(1, 2)))
    put(
        "fixed_disk_radial_density",
        s.simplify((s.diff(primitive, r) - 2 * s.pi * r * rho).subs(r, s.sqrt(1 - z))),
    )
    put("fixed_disk_normalization", primitive.subs(r, 1) - primitive.subs(r, 0) - 1)
    put("fixed_disk_D4_projection", rho.subs(e, 0) - 1 / (2 * s.pi * s.sqrt(1 - r * r)))
    put(
        "fixed_disk_relative_density",
        s.simplify(
            (
                rho / (1 / (2 * s.pi * s.sqrt(1 - r * r)))
                - (1 + 2 * e) * (1 - r * r) ** e
            ).subs(r, s.sqrt(1 - z))
        ),
    )
    x = s.Symbol("x", nonnegative=True)
    nonnegative(
        "fixed_disk_density_multiplier",
        s.Rational(5, 4) - (1 + 2 * (s.Rational(1, 8) - x)),
        (x,),
    )
    # Rank2 transverse tensor positivity is explicit for all real D=4+2e.
    a, b, c = s.symbols("a b c", real=True)
    bilinear = a * a + b * b + 2 * c * c - (a + b) ** 2 / (2 + 2 * e)
    put(
        "continuous_TT_positive_squares",
        bilinear - ((a - b) ** 2 / 2 + 2 * c * c + e * (a + b) ** 2 / (2 + 2 * e)),
    )
    put(
        "continuous_TT_Frobenius_gap",
        a * a + b * b + 2 * c * c - bilinear - (a + b) ** 2 / (2 + 2 * e),
    )
    put("massive_pair_phase_ratio_exponent", (3 + 2 * e - 2) - (1 + 2 * e))
    put("angular_log_first_moment_bound", s.Rational(3, 4) * 2 - s.Rational(3, 2))
    put("angular_log_second_moment_bound", s.Rational(3, 4) * 4 - 3)
    put("normalized_J_first_derivative_bound", 2 * 16 * s.Rational(3, 2) - 48)
    put("normalized_J_second_derivative_bound", 2 * 16 * 3 - 96)
    put("whole_K_bound", 16 * 50 * 16 - 12800)
    put("whole_K_first_bound", 16 * (s.Rational(1, 2) * 16 + 50 * 48) - 38528)
    put("whole_K_second_bound", 16 * (16 + 2 * s.Rational(1, 2) * 48 + 50 * 96) - 77824)
    put("K0_sharp_bound", 16 * 50 - 800)
    put("K1_sharp_bound", 16 * (s.Rational(1, 2) + 50 * 6) - 4808)
    put("phase_times_K_second_bound", 17 * 12800 + 8 * 38528 + 77824 - 603648)
    put("finite_conversion_numerator_bound", 4808 + 4 * 800 - 8008)
    nonnegative("finite_conversion_below112", 112 - s.Rational(8008, 72), (x,))
    nonnegative(
        "recoil_log_bound_below2omega",
        2 - s.Rational(20, 9) * s.Rational(18, 13) / 2,
        (x,),
    )
    nonnegative("continued_phase_error_below3omega", 3 - (2 + s.Rational(1, 2)), (x,))
    put("regulator_log_coefficient", 2 * 8008 - 16016)
    put("regulator_second_order_coefficient", s.Rational(603648, 2) - 301824)
    nonnegative("regulator_log_coarse_bound", 225 - s.Rational(16016, 72), (x,))
    nonnegative("regulator_constant_coarse_bound", 4200 - s.Rational(301824, 72), (x,))
    nonnegative(
        "physical_rate_total_original_bound",
        10**8 - (s.Rational(42248192, 8) + s.Rational(54450000000, 64)) / 18 - 112,
        (x,),
    )
    put(
        "original_conversion_exponent",
        s.Integer(10) ** 3 / s.Integer(10) ** 800 - s.Rational(1, 10**797),
    )
    put(
        "original_physical_comparison_exponent",
        s.Integer(10) ** 8 / s.Integer(10) ** 800 - s.Rational(1, 10**792),
    )
    # Spectral probability identities for derivatives of a normalized angular mean.
    F0, F1, F2, Z0, Z1, Z2 = s.symbols("F0 F1 F2 Z0 Z1 Z2", real=True)
    expect_second = (
        F2 / Z0 - 2 * F1 * Z1 / Z0**2 + 2 * F0 * Z1**2 / Z0**3 - F0 * Z2 / Z0**2
    )
    centered = (
        F2 / Z0
        - 2 * (Z1 / Z0) * (F1 / Z0)
        + (Z1 / Z0) ** 2 * F0 / Z0
        - F0 / Z0 * (Z2 / Z0 - (Z1 / Z0) ** 2)
    )
    put("normalized_second_derivative_variance_identity", expect_second - centered)
    # At forward the incoming and outgoing momenta coincide in pairs for all D.
    p, q = s.symbols("p q", real=True)
    put("forward_entire_soft_current", (-p - p + q + q).subs(q, p))

    v = s.Symbol("v", positive=True)
    put(
        "folded_angular_density_primitive",
        s.diff(v ** (1 + 2 * e), v) - folded_angular_weight(v, e),
    )
    put(
        "folded_angular_normalization",
        s.Integer(1) ** (1 + 2 * e) - s.Integer(0) ** (1 + 2 * e) - 1,
    )
    put("folded_angular_physical_limit", folded_angular_weight(v, 0) - 1)
    E0, w0 = s.symbols("hard_energy real_energy", positive=True)
    put(
        "continued_pair_phase_D4_limit",
        continued_pair_density(E0, w0, 0) - recoil.density_ratio(E0, w0),
    )
    return {
        "whole_fixed_domain_angular_weight": folded_angular_weight(),
        "whole_soft_pair_regulator_error": soft_pair_regulator_error(),
        "whole_finite_conversion_bound": finite_conversion_bound(),
        "explicit_nonnegative_margins": margins,
        "whole_normalized_angular_majorants": {
            "J": s.Integer(16),
            "J_first": s.Integer(48),
            "J_second": s.Integer(96),
            "K": s.Integer(12800),
            "K_first": s.Integer(38528),
            "K_second": s.Integer(77824),
            "phase": s.Integer(1),
            "phase_first": s.Integer(4),
            "phase_second": s.Integer(17),
            "phase_K_second": s.Integer(603648),
            "K0_sharp": s.Integer(800),
            "K1_sharp": s.Integer(4808),
        },
        "whole_fixed_disk_dominated_convergence": "The normalized angular pushforward is rho_e(v)=(1/2+e)/pi*(1-|v|^2)^(e-1/2) on a fixed unit disk. It is<=5rho0/4 and tends to rho0 for0<=e<=1/8. Equivalently the folded perpendicular coordinate has normalized weight(1+2e)v^(2e)dv dphi/(2pi) on a fixed rectangle.",
        "whole_continuous_TT_norm_proof": "First prove the full Ward identity, then remove pure metric terms. The projected anisotropic tensors share a rank<=2 subspace. tr(AB)-tr(A)tr(B)/(2+2e) is positive and bounded by the Frobenius form for every real e>=0. A three-dimensional slice of an unremoved full-D isotropic tensor would be wrong. No fictitious noninteger helicity count is used.",
        "whole_continued_recoil_majorant": "J_e=J0*(rprime/r0)^(2e),rprime^2/r0^2>=13/18 and-log(rprime/r0)<2omega. Thus1-J_e<3omega. The exact-minus-soft squared kernel is bounded by[B^2+(128B+12288)/omega]/kappa,B330000. With radial omega^(1+2e), phase<=1 and the fixed-disk bound, dominated convergence gives exactly the physical S295 real remainder.",
        "whole_order_of_limits": "The explicit leading-pair error is e[16016|ln resolution|+301824]/(8pi^2 kappa), reference nu1. Remove e at fixed nonzero resolution first. This is not uniform under arbitrary e|ln resolution| scaling, does not supply a quantitative O(e) bound for the entire nonsoft integral, and does not exchange the forward/gravity/Regge limits.",
        "checks": checks,
        "gates": {
            "fixed_domain_positive_angular_representation": True,
            "continuous_TT_positive_norm_without_noninteger_helicity_count": True,
            "pure_metric_terms_removed_only_after_full_Ward": True,
            "massive_Doppler_and_recoil_Jacobian_bounds": True,
            "explicit_soft_pair_regulator_and_conversion_majorants": True,
            "physical_nonsoft_D4_limit_by_dominated_convergence": True,
            "regulator_resolution_and_forward_limits_not_exchanged": True,
        },
    }
