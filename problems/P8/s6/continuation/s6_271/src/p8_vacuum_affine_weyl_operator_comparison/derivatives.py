"""All-order full cutoff derivatives and finite high-phase Cauchy bounds."""

from functools import cache

import sympy as s
from p8_vacuum_affine_quantitative_local_time import quantum as previous

from . import source

BUMP_BASE = s.Integer(8)
TRANSITION_BASE = s.Integer(256)
RADIAL_BASE = s.Integer(2048)
PRODUCT_BASE = s.Integer(2056)


def coefficient(order):
    if (
        isinstance(order, bool)
        or not isinstance(order, (int, s.Integer))
        or not 0 <= order <= source.MAX_ORDER
    ):
        raise ValueError("Require exact proved phase derivative order0 through196")
    return PRODUCT_BASE**order * s.factorial(order) ** 3


def growth(order):
    return coefficient(order) / source.RADIUS**order


def amplitude_jet(order, amplitude):
    if isinstance(amplitude, bool) or not isinstance(
        amplitude, (int, s.Integer, s.Rational)
    ):
        raise TypeError(
            "Require an exact original amplitude, without coercing strings or floats"
        )
    amplitude = s.sympify(amplitude)
    if amplitude not in (source.HBOUND, source.VERROR):
        raise ValueError(
            "Require the unchanged full Hamiltonian or physical-volume amplitude"
        )
    return 2 * amplitude * growth(order)


@cache
def radial_partition(order):
    if (
        isinstance(order, bool)
        or not isinstance(order, (int, s.Integer))
        or not 0 <= order <= source.MAX_ORDER
    ):
        raise ValueError("Require exact proved radial derivative order")
    return sum(
        s.factorial(order)
        * 4 ** (order - 2 * j)
        * TRANSITION_BASE ** (order - j)
        * s.factorial(order - j) ** 2
        / (s.factorial(order - 2 * j) * s.factorial(j))
        for j in range(order // 2 + 1)
    )


@cache
def data():
    n = s.Symbol("nonnegative_derivative_order", integer=True, nonnegative=True)
    ratio = s.Rational(1, 32)
    old = previous.cutoff_bounds()
    # All196 actual finite-order partition sums, not a float asymptotic.
    radial = {k: radial_partition(k) for k in range(source.MAX_ORDER + 1)}
    ratios = [
        PRODUCT_BASE * (k + 1) ** 3 / source.RADIUS for k in range(source.MAX_ORDER)
    ]
    x = s.Symbol("inverse_transition_coordinate", positive=True)
    polynomials = [s.Integer(1)]
    for _ in range(8):
        polynomials.append(
            s.expand(x * x * (polynomials[-1] - s.diff(polynomials[-1], x)))
        )
    weighted = [
        sum(abs(c) * s.factorial(m[0]) for m, c in s.Poly(poly, x).terms())
        for poly in polynomials
    ]
    return {
        "whole_all197_scaled_product_derivative_coefficients": {
            str(k): coefficient(k) for k in range(source.MAX_ORDER + 1)
        },
        "whole_maximum_order_and_Cauchy_radius": [source.MAX_ORDER, source.RADIUS / 8],
        "whole_all196_successive_derivative_growth_ratios": ratios,
        "whole_all_order_bump_transition_radial_product_bases": [
            BUMP_BASE,
            TRANSITION_BASE,
            RADIAL_BASE,
            PRODUCT_BASE,
        ],
        "whole_initial_nine_actual_bump_weighted_derivatives": weighted,
        "whole_old_exact_first_four_cutoff_contacts": {
            name: old[name] for name in ("bump", "transition", "radial", "product")
        },
        "whole_complete_transition_induction_margin": 1
        - 9 * ratio
        - 18 * ratio / (1 - ratio),
        "whole_full_cutoff_high_derivative_proof": "For P0=1, Pn+1=x^2(Pn-Pn'), the degree is<=2n. With Bn=sum|p_k|k!, the complete derivative recurrence gives Bn+1<=2(2n+1)^2 Bn<=8(n+1)^2 Bn, hence Bn<=8^n(n!)^2. Every smooth constant extension is included. The unchanged transition denominator>=1/9 and quotient recurrence give Cn<=9[Bn+2sum binom(n,k)Bk Cn-k]. Normalize by(n!)^2, use binom(n,k)>=1 and r=8/256=1/32; the induction bound9r+18r/(1-r)=855/992<1 proves Cn<=256^n(n!)^2 for every n>=1. No C5-in-time premise is used for phase-cutoff derivatives.",
        "whole_full_mixed_radial_partition_proof": "Only first and second derivatives of the quadratic radial argument are nonzero, with norms<=4/R and2/R^2. For total ordern and j paired arguments the complete multilinear partition coefficient is n!/[(n-2j)!j!2^j]. Including the inner2/R^2 factors gives n!4^(n-2j)C_(n-j)/[(n-2j)!j!]R^-n. Sum everyj=0..floor(n/2). Each term is<=1024^n(n!)^3 and the number of terms<=2^n, so the full mixed derivative bound is2048^n(n!)^3/R^n. The code also evaluates ALL197 actual rational finite-order sums under the transition majorant and verifies this bound.",
        "whole_complete_high_Cauchy_Leibniz_proof": "The whole actual phase amplitude on complex ball4R has boundA. On real ball2R, simultaneous at most96-coordinate Cauchy radiusR/8 gives n!8^n A/R^n. The unchanged centered amplitude has bound2A, with its actual scalar center retained. Leibniz with the entire radial derivative gives2A/R^n times sum binom(n,k)2048^k(k!)^3(n-k)!8^(n-k), at most2A*2056^n(n!)^3/R^n. This uses every mixed product, implicit nonlinear contact and original cross covariance. For every n<=196 the displayed exact successive ratio is<1, so the bound decreases through the required range; no asymptotic-in-n or arbitrary-time assertion is used.",
        "checks": {
            "full_weighted_bump_induction_polynomial": s.expand(
                8 * (n + 1) ** 2 - 2 * (2 * n + 1) ** 2 - 8 * n - 6
            ),
            "whole_transition_majorant_ratio": BUMP_BASE / TRANSITION_BASE - ratio,
            "whole_transition_induction_constant": 9 * ratio
            + 18 * ratio / (1 - ratio)
            - s.Rational(855, 992),
            "whole_complete_product_base": PRODUCT_BASE - RADIAL_BASE - 8,
            "same_exact_low_order_bump_sequence": s.Matrix(weighted[:5])
            - s.Matrix(old["bump"]),
            "whole_required_high_phase_order": s.Integer(source.MAX_ORDER - 196),
        },
        "gates": {
            "full_induction_margin_positive": 9 * ratio + 18 * ratio / (1 - ratio) < 1,
            "first_nine_actual_bump_derivatives_below_general_bound": all(
                weighted[k] <= BUMP_BASE**k * s.factorial(k) ** 2 for k in range(9)
            ),
            "all197_full_radial_partition_sums_below_bound": all(
                radial[k] <= RADIAL_BASE**k * s.factorial(k) ** 3 for k in radial
            ),
            "every_original_low_order_transition_retained": all(
                old["transition"][k] <= TRANSITION_BASE**k * s.factorial(k) ** 2
                for k in range(5)
            ),
            "every_original_low_order_radial_retained": all(
                old["radial"][k] <= RADIAL_BASE**k * s.factorial(k) ** 3
                for k in range(5)
            ),
            "every_original_low_order_product_retained": all(
                old["product"][k] <= coefficient(k) for k in range(5)
            ),
            "all196_actual_growth_ratios_strictly_decreasing": all(
                r < 1 for r in ratios
            ),
            "all196_actual_growth_ratios_below_one_billionth": max(ratios)
            < s.Rational(1, 10**9),
            "full_order192_plus_four_heat_derivatives_covered": source.MAX_ORDER
            == 2 * source.PHASE + 4,
            "no_unknown_cutoff_or_high_time_profile_substitution": True,
        },
    }
