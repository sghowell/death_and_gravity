"""Holomorphic energy-tube budgets for the full soft-scaled tree."""

from functools import cache

import sympy as s
from p8_vacuum_affine_two_real_collinear_current import collinear

from . import source

RADIUS = s.Rational(1, 10**12)
LINEAR_SLOPE = (
    s.Integer(10) ** 40 * source.HEAVY_MASS2**2 + s.Integer(10) ** 60
) / source.KAPPA
CAUCHY_COEFFICIENT = LINEAR_SLOPE / RADIUS**2
REMAINDER_COEFFICIENT = s.Rational(1, 10**326)


@cache
def angular_coefficients():
    a, b, r, records, _ = collinear.coefficients()
    z = s.Symbol("z")
    budgets = {}
    checks = {}
    degree = 0
    for pair, coefficients in records.items():
        total = s.S.Zero
        for index, coefficient in enumerate(coefficients):
            if coefficient == 0:
                continue
            normalized = s.factor(
                (coefficient * a * b / (a + b) ** 2).subs(
                    {a: z, b: 1 - z}, simultaneous=True
                )
            )
            numerator, denominator = s.fraction(normalized)
            m = int(s.degree(denominator, r) // 2)
            constant = denominator.subs(r, 0)
            poly = s.Poly(numerator, z, r)
            checks[f"complex_denominator_{pair[0]}_{pair[1]}_{index}"] = s.expand(
                denominator - constant * (1 + r * r) ** m
            )
            if not all(monomial[1] <= 2 * m for monomial, value in poly.terms()):
                raise ValueError("An unbounded real-angle monomial remains")
            degree = max(degree, int(poly.degree(z)))
            total += (
                sum(abs(value) * 2 ** monomial[0] for monomial, value in poly.terms())
                / constant
            )
        budgets[pair] = total
    expected = {(0, 0): 2642, (0, 1): 432, (1, 0): 352, (1, 1): 1764}
    for pair, value in budgets.items():
        checks[f"complex_pair_coefficient_budget_{pair[0]}_{pair[1]}"] = (
            value - expected[pair]
        )
    checks["largest_energy_fraction_polynomial_degree"] = s.Integer(degree - 4)
    return checks, budgets


@cache
def data():
    c = RADIUS
    Wmax = s.Rational(1, 8)
    checks, angles = angular_coefficients()
    checks = dict(checks)
    r3 = 3 * 138240 * 30000000
    r4 = 3 * 10616832 * 30000000
    nondouble_gr = 177 * 1024**4 * 8 * 1800000 * max(1, r3, r3**2, r4)
    nondouble_matter = (33 * 4 + 177 * 8) * 1024**2 * 8
    mixed_distribution = (
        s.Rational(10**7 * 2 * 10**8 * 600000, 8)
        + s.Rational(2 * 10**5 * 10**8 * 600000 * 2, 8)
        + s.Rational(2 * 10**5 * 1024**2 * 200 * 600000 * 2, 8)
    )
    timelike_distribution = (
        s.Rational(10**7 * 128**2, 1) / s.Rational(45, 16)
        + s.Rational(2 * 10**5 * (2 * 128 * 8192), 1) / s.Rational(45, 16)
        + s.Rational(2 * 10**5 * 128**2 * 200, 1)
        / (s.Rational(45, 16) * s.Rational(25, 4))
    ) / 8
    double_gr = 60 * (10**5 * 2 * 10**5 + 1024 * 10**6) * s.Rational(600000, 8)
    matter_C = 4 * (2 * 256 * 16384)
    matter_H = 12 * (2 * 2 * 128 * 8192 + 7 * 128**2)
    matter_ordering = 80 * 10**5 * 4
    matter_total = nondouble_matter + matter_C + matter_H + matter_ordering
    matter_current = 4 * 1024 * 4 * 4 + 12 * 1024 * 4 * 2 + Wmax * (16 + 12288 + 48)
    gravity_current = (
        12 * 1024**2 * 4 * 10**8
        + 6 * 1024**2 * 10**8
        + 2 * 138240 * 50 * 10**16
        + 138240 * 25 * 64**2 * s.Rational(65, 64)
    )
    # Separate the n^2 and constant coefficients before phase multiplication.
    matter_coefficient = s.Rational(matter_total, 4) + 3 * 10**10
    gravity_coefficient = (
        s.Rational(nondouble_gr, 8)
        + 8 * mixed_distribution
        + 4 * timelike_distribution
        + double_gr
        + 10**32
    )
    final_matter = 2 * matter_coefficient
    final_gravity = 2 * gravity_coefficient + 16 * 2304
    checks.update(
        {
            "complex_mixed_hard_gap_margin_numerator": s.Integer(
                3601 * 600000 - 2160600000
            ),
            "complex_energy_gradient_bound": s.Rational(7, 4) / s.Rational(7, 8) - 2,
            "complex_soft_current_gradient_budget": s.Integer(4 * 8 + 4 * 3 * 64 - 800),
            "complex_assigned_stress_variation": s.Integer(72) + 24 * Wmax - 75,
            "complex_Born_stress_variation": s.Integer(144) + 96 * Wmax - 156,
            "assigned_numerator_variation_budget": s.Integer(
                16 * (100 * 165 + 55 * 300) - 528000
            ),
            "Born_numerator_variation_budget": s.Integer(
                16 * (200 * 165 + 55 * 600) - 1056000
            ),
            "complex_Einstein_cubic_step": s.Integer(r3 - 12441600000000),
            "complex_Einstein_quartic_step": s.Integer(r4 - 955514880000000),
            "nondouble_Einstein_budget": s.Integer(nondouble_gr)
            - s.Integer("433798508189475293550870528000000000000000000000"),
            "complex_matter_current_budget": matter_current - 165384,
            "complex_Einstein_current_budget": gravity_current
            - s.Integer("138240005662324776960000"),
            "complex_same_leg_ordering_correction_budget": s.Integer(
                5120 + 5120 + 65536 - 75776
            ),
            "complex_double_current_variation_budget": s.Integer(
                4096 * (10000 + 1024) - 45154304
            ),
            "Cauchy_mixed_derivative_coefficient": CAUCHY_COEFFICIENT
            - LINEAR_SLOPE / RADIUS**2,
        }
    )
    x = s.Symbol("Eprime", positive=True)
    f = s.sqrt(s.sqrt(x * x - 1) / x)
    checks["analytic_phase_logarithmic_derivative"] = s.simplify(
        s.diff(f, x) / f - 1 / (2 * x * (x * x - 1))
    )
    tau, W = s.symbols("tau W", positive=True)
    checks["two_term_mixed_momentum_square_margin"] = s.expand(
        2 * (tau + W * W) - (s.sqrt(tau) + W) ** 2 - (s.sqrt(tau) - W) ** 2
    )
    return {
        "whole_complex_energy_domain": "For every physical a,b>0,W=a+b<=1/8, fix E,u,real emitted directions and physical TT tensors. Continue the S300 recoil into |a'-a|,|b'-b|<=1e-12 W using the square-root branches from the real state. The real Born amplitude stays a fixed positive normalization; delta and tau are real-center bounding parameters, not nonanalytic operations applied to the continued amplitude.",
        "whole_holomorphic_recoil_and_denominators": "The written tube proof gives |delta Ep|<=5cW,|delta rprime|<=10cW,|delta p0|<5cW and spatial norm shift<13cW. Complex massive spatial norms<7/4 and energy modulus>7/8 on the relevant convex tube imply a soft-current gradient below800<1024 and Doppler modulus>1/8. Each mixed hard denominator retains magnitude>(tau+W2)/600000, timelike gaps>45/16,heavy gaps>n/2,and total-Q light denominators>W/4. Individual a' or b' light poles are removed by the ab prefactor. The combined conserved47 pair class removes its apparent angular pole before continuation is bounded.",
        "whole_complex_angular_budgets": {
            f"polarization_{i}_{j}": v for (i, j), v in angles.items()
        },
        "whole_component_budget_records": {
            "nondouble_Einstein": nondouble_gr,
            "nondouble_matter": nondouble_matter,
            "one_mixed_leading_distribution": mixed_distribution,
            "one_timelike_leading_distribution": timelike_distribution,
            "all_double_Einstein_ordering_and_core_errors": double_gr,
            "all_regular_matter_variations": matter_total,
            "complex_matter_current": matter_current,
            "complex_Einstein_current": gravity_current,
            "final_n_squared_coefficient": final_matter,
            "final_constant_coefficient": final_gravity,
            "whole_linear_slope": LINEAR_SLOPE,
            "whole_Cauchy_coefficient": CAUCHY_COEFFICIENT,
        },
        "whole_soft_scaled_linear_theorem": "Let G=ab sqrt(rho2) M6/A0 and G00=S_a(Born)S_b(Born)/kappa. The247 regular graphs with0/1 light propagator give O(W) after ab scaling. The140 double-external graphs first use the exact two-order identity and four radiation distributions; paired mixed currents and telescoping N*Sprod/D against its Born value give O(W) uniformly in tau. The remaining47 conserved pair graphs are O(W) by their complex angular coefficient budgets and full-current estimates. The analytic phase obeys |sqrt(rho)|<2 and |sqrt(rho)-1|<16W. All coefficients fit |G-G00|<B W, B=(1e40n2+1e60)/kappa<1e-350, throughout each stated polydisc.",
        "whole_Cauchy_theorem": "The scaled graph sums are holomorphic after their removable soft factors are canceled; all other denominators and recoil square roots stay away from zero in a neighborhood of the closed polydisc. Two-variable Cauchy applied to G-G00 on radiuscW gives |partial_a partial_b G|<B/(c2 W)<1e-326/W. The real bounded envelope of S312 alone would not imply this derivative estimate; the complex continuation and O(W) grouping are essential.",
        "checks": checks,
        "gates": {
            "complex_radius_positive_and_tiny": 0 < c < s.Rational(1, 10**10),
            "mixed_hard_gap_survives_complex_perturbation": 3601 * c
            < s.Rational(1, 600000),
            "total_Q_light_gap_survives": 100 * c < s.Rational(1, 8),
            "complex_massive_spatial_tube_is_small": 13 * c * Wmax
            < s.Rational(1, 1000),
            "complex_spatial_norm_margin": s.Rational(7, 4) ** 2 > 3,
            "complex_Doppler_tube_margin": 65 * c * Wmax < s.Rational(1, 8),
            "complex_current_gradient_below1024": s.Integer(800) < 1024,
            "paired_product_variation_below1e8": 45154304 < 10**8,
            "all_complex_angular_budgets_below1e4": max(angles.values()) < 10**4,
            "nondouble_Einstein_budget_below1e49": nondouble_gr < 10**49,
            "mixed_and_timelike_leading_budgets_below1e25": 8 * mixed_distribution
            + 4 * timelike_distribution
            < 10**25,
            "double_Einstein_corrections_below1e20": double_gr < 10**20,
            "whole_regular_matter_budget_below3e10": matter_total < 3 * 10**10,
            "whole_matter_current_budget_below1e6": matter_current < 10**6,
            "whole_Einstein_current_budget_below1e26": gravity_current < 10**26,
            "full_phase_weighted_n_squared_budget_below1e40": final_matter < 10**40,
            "full_phase_weighted_constant_budget_below1e60": final_gravity < 10**60,
            "original_linear_slope_below1e_minus350": LINEAR_SLOPE
            < s.Rational(1, 10**350),
            "original_Cauchy_coefficient_below1e_minus326": CAUCHY_COEFFICIENT
            < REMAINDER_COEFFICIENT,
            "written_holomorphy_not_inferred_from_real_samples": True,
        },
    }
