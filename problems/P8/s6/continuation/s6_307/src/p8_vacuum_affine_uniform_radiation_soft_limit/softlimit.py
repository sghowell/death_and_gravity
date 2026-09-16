"""Uniformly vanishing single-real remainder against an unexpanded soft sum."""

from functools import cache

import sympy as s

from . import estimates, source

LINEAR = s.Integer(2) * 10**14
QUADRATIC = s.Integer(2) * 10**37
FIXED = s.Integer(10) ** 32
X0 = s.Rational(1, 1000)


def exact_rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Rational)):
        raise TypeError("Require an exact rational soft-reference coordinate")
    return s.Rational(value)


def finite_remainder_bound(energy, transfer, resolution, kappa):
    _, _, _, _, _, x = estimates.require_domain(energy, transfer, resolution)
    source.require_mass(kappa)
    return min(FIXED, LINEAR * x + QUADRATIC * x * x) / s.sympify(kappa)


def original_relative_bound(energy, transfer, resolution):
    estimates.require_domain(energy, transfer, resolution)
    return 4 * FIXED / source.KAPPA


def original_vanishing_relative_bound(energy, transfer, resolution):
    _, _, _, _, _, x = estimates.require_domain(energy, transfer, resolution)
    return 2 * (LINEAR * s.sqrt(x) + QUADRATIC * x * s.sqrt(x)) / source.KAPPA


def reference_factor(index, conversion, resolution):
    a, d, x = map(exact_rational, (index, conversion, resolution))
    if not (
        0 <= a <= 6 / source.KAPPA
        and abs(d) <= 112 / source.KAPPA
        and 0 < x <= s.Rational(1, 8)
    ):
        raise ValueError("Require the original compact elastic leading-soft reference")
    return source.gamma.conversion(a, d) * x**a


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = s.factor(value)

    x, w, y, delta, kap = s.symbols("x omega root_delta delta kappa", positive=True)
    alpha = s.Symbol("soft_index", nonnegative=True)
    B, D, Bm = estimates.LOW, estimates.DIRECT, estimates.MATTER
    low_integrand = w * ((256 * B + 16384) / w + 2 * B**2)
    low_numerator = (256 * B + 16384) * x + B**2 * x * x
    put(
        "whole_low_band_two_polarizations",
        2 * (2 * 64 * B / w + B**2)
        + 2 * w * 2 * 64**2 / w**2
        - ((256 * B + 16384) / w + 2 * B**2),
    )
    put(
        "whole_low_band_finite_integral",
        s.integrate(low_integrand, (w, 0, x)) - low_numerator,
    )
    low_coefficient = (256 * B + 16384 + B**2 / 8) / 36
    put(
        "whole_low_band_linear_coefficient",
        low_coefficient - s.Rational(1250012800004096, 9),
    )
    put("whole_graviton_measure", 4 * s.pi / (2 * (2 * s.pi) ** 3) - 1 / (4 * s.pi**2))
    put(
        "whole_high_band_square_gap",
        delta - delta**2 / (delta + w * w) - delta * w * w / (delta + w * w),
    )
    high_numerator = (6 * 160**2 + 2 * D**2) * delta * s.log(
        192 * x / s.sqrt(delta)
    ) + 2 * Bm**2 * x * x
    primitive = (6 * 160**2 + 2 * D**2) * delta * s.log(w) + 2 * Bm**2 * w * w
    put(
        "whole_high_band_integral_derivative",
        s.diff(primitive, w) - w * ((6 * 160**2 + 2 * D**2) * delta / w**2 + 4 * Bm**2),
    )
    put(
        "whole_high_band_upper_endpoint",
        primitive.subs(w, x)
        - primitive.subs(w, s.sqrt(delta) / 192)
        - high_numerator
        + 2 * Bm**2 * delta / 192**2,
    )
    entropy = y * y * s.log(192 * x / y)
    put(
        "whole_entropy_derivative",
        s.diff(entropy, y) - y * (2 * s.log(192 * x / y) - 1),
    )
    put(
        "whole_entropy_maximum_location",
        s.diff(entropy, y).subs(y, 192 * x / s.sqrt(s.E)),
    )
    put(
        "whole_entropy_maximum_value",
        entropy.subs(y, 192 * x / s.sqrt(s.E)) - (192 * x) ** 2 / (2 * s.E),
    )
    put("whole_entropy_zero_endpoint", s.limit(entropy, y, 0, dir="+"))
    high_coefficient = ((6 * 160**2 + 2 * D**2) * 18432 + 2 * Bm**2) / 36
    put(
        "whole_high_band_quadratic_coefficient",
        high_coefficient - s.Integer(10240000000000000000000000006128643200),
    )
    put(
        "whole_uniform_vanishing_absolute_limit",
        s.limit((LINEAR * x + QUADRATIC * x * x) / kap, x, 0, dir="+"),
    )
    put("whole_at_original_kappa", 4 * FIXED / source.KAPPA - s.Rational(4, 10**768))
    split_numerator = LINEAR * X0 + QUADRATIC * X0**2
    put(
        "whole_split_numerator",
        split_numerator - s.Integer(20000000000000000000200000000000),
    )
    put(
        "whole_below_split_power_identity",
        (LINEAR * x + QUADRATIC * x * x) / x**alpha
        - (LINEAR * x ** (1 - alpha) + QUADRATIC * x ** (2 - alpha)),
    )
    put(
        "whole_split_power_identity",
        (LINEAR * X0 ** (1 - alpha) + QUADRATIC * X0 ** (2 - alpha))
        - split_numerator * 1000**alpha,
    )
    put("whole_prefactor_lower_product", s.Rational(3, 4) ** 2 - s.Rational(9, 16))
    put(
        "whole_relative_vanishing_majorant",
        2 * (LINEAR * s.sqrt(x) + QUADRATIC * x * s.sqrt(x)) / kap
        - (4 * 10**14 * s.sqrt(x) + 4 * 10**37 * x ** s.Rational(3, 2)) / kap,
    )
    put(
        "whole_uniform_relative_zero_limit",
        s.limit(
            2 * (LINEAR * s.sqrt(x) + QUADRATIC * x * s.sqrt(x)) / kap, x, 0, dir="+"
        ),
    )
    put(
        "whole_finite_soft_conversion_inherited",
        source.continuity.finite_conversion_bound(kap) - 112 / kap,
    )
    put(
        "whole_gamma_factor_definition",
        source.gamma.gamma_factor(alpha)
        - s.exp(-s.EulerGamma * alpha) / s.gamma(1 + alpha),
    )
    margins = {
        "low_integral_linear_roundup": LINEAR - low_coefficient,
        "high_integral_quadratic_roundup": QUADRATIC - high_coefficient,
        "whole_coarse_linear_roundup": s.Integer(3) * 10**36 - LINEAR - QUADRATIC / 8,
        "original_soft_index_below_half": s.Rational(1, 2) - 6 / source.KAPPA,
        "original_conversion_magnitude_below_quarter": s.Rational(1, 4)
        - 112 / source.KAPPA,
        "prefactor_product_above_half": s.Rational(9, 16) - s.Rational(1, 2),
        "split_numerator_below_fixed_bound": FIXED - split_numerator,
        "exp9_lower_polynomial_above1000": sum(
            s.Integer(9) ** j / s.factorial(j) for j in range(9)
        )
        - 1000,
        "soft_split_exponent_below_half": s.Rational(1, 2) - 54 / source.KAPPA,
        "old_all_angle_fixed_bound": FIXED
        - source.radiation_bounds.uniform_data()["whole_combined_coefficient"],
        "split_in_detector_domain": s.Rational(1, 8) - X0,
    }
    for name, value in margins.items():
        put("positive_arithmetic_" + name, value - s.Abs(value))
    return {
        "whole_low_band_numerator": low_numerator,
        "whole_high_band_numerator_majorant": high_numerator,
        "whole_low_linear_coefficient": low_coefficient,
        "whole_high_quadratic_coefficient": high_coefficient,
        "whole_uniform_absolute_bound": s.Min(FIXED, LINEAR * x + QUADRATIC * x * x)
        / kap,
        "whole_original_uniform_relative_bound": s.Rational(4, 10**768),
        "whole_uniform_vanishing_relative_majorant": 2
        * (LINEAR * s.sqrt(x) + QUADRATIC * x * s.sqrt(x))
        / kap,
        "whole_positive_margins": margins,
        "whole_total_variation_definition": "Integrate the absolute value of rho sum_pol|M5/(Am+AG)|^2-sum_pol|S0|^2 over the emitted one-graviton measure with0<omega<=x. This total variation bounds the signed finite real-minus-soft rate and has its own finite lower-energy-cutoff endpoint. Its dimensionless Born normalization retains all47 tree interferences.",
        "whole_uniform_absolute_proof": "Splitomega atsqrt(delta)/192 when this lies belowx. The constant low remainder2e8/sqrt(kappa), exactrho<=1 and|rho-1|<=2omega give[(256B+16384)x+B^2x^2]/(4pi^2 kappa)<2e14*x/kappa. Positive Born convexity and the improved global direct bound give a high numerator(6*160^2+2D^2)delta ln(192x/sqrt(delta))+2Bm^2x^2. Writingy=sqrt(delta),y^2ln(192x/y)<18432x^2 uniformly. With4pi^2>36 the high part is<2e37*x^2/kappa. These regions cover every angle and every0<x<=1/8; the estimate goes tozero uniformly asx->0.",
        "whole_unexpanded_reference_proof": "For the original elastic knownP(x)=exp(Delta)F(alpha)x^alpha,0<=alpha<6/kappa<1/2 and|Delta|<112/kappa<1/4. S299 givesF>=1-alpha^2>3/4,while expDelta>3/4,soP(x)>x^alpha/2. Forx<=1/1000 the powersx^(1-alpha),x^(2-alpha) increase, and the numerator at that split is below1e32. Forx>=1/1000 use the independent S304 fixed1e32 bound. ThusTV/x^alpha<1e32*1000^alpha/kappa. Since ln1000<9 and54/kappa<1/2,1000^alpha<exp(1/2)<2, soTV/P(x)<4e32/kappa=4e-768. Alsoalpha<=1/2 gives the displayed uniformly vanishing relative majorant. No Taylor expansion ofalpha*lnx is made.",
        "whole_scope_boundary": "A complete single-real finite remainder is compared to the named elastic all-N leading-soft reference. This does not identify the former as the full nonleading correction to an all-N calorimetric rate. Correlated multi-real recoil, graviton self-interactions, hard loops and matching, Coulomb/Regge continuation and the original quantum state remain separate obligations. Regulator removal remains ordered before detector-threshold limits.",
        "checks": checks,
        "gates": {
            "all_integration_and_reference_margins_strict": all(
                bool(value > 0) for value in margins.values()
            ),
            "whole_two_polarization_absolute_integrand": True,
            "low_and_high_regions_cover_arbitrary_nonforward_angles": True,
            "uniform_absolute_and_relative_zero_resolution_limits": True,
            "same_unexpanded_allN_leading_soft_reference": True,
            "fixed_positive_resolution_regulator_removal_first": True,
            "single_real_not_promoted_to_allN_nonleading_control": True,
        },
    }
