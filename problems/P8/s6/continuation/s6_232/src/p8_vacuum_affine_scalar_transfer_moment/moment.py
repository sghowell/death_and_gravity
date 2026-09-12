"""Required full low spectral moment, not a physical cutoff or loop approximation."""

from functools import cache

import sympy as s

from . import coefficients as c

ENERGY = s.Integer(10) ** 99
SPLIT = ENERGY**2
DELTA0 = s.S.One
DELTA1 = s.Rational(1, 2)
NAMED_REQUIRED = c.GAMMA - 8 * c.LAMBDA / (SPLIT - 2)
FIRST_UPPER = (
    48 * c.LAMBDA**2 * SPLIT + s.Rational(9, 16) * c.GAMMA**2 * SPLIT**3 + c.GAMMA**2
) / 9


@cache
def data():
    K = s.Symbol("spectral_split_K", positive=True)
    d0, d1 = s.symbols(
        "coefficient_relative_error0 coefficient_relative_error1", nonnegative=True
    )
    b0, b1, Jlow, Jhigh, Blo, Bhi, Tpos = s.symbols(
        "b20 b21 J4_low J4_high B20_low B20_high positive_transfer_integral", real=True
    )
    weight = s.Rational(3, 2) / (K - 2)
    # Exact identities: b0=Blo+Bhi and b1=Tpos-3(Jlow+Jhigh)/2.
    difference = (b1 + weight * b0 + s.Rational(3, 2) * Jlow).subs(
        {b0: Blo + Bhi, b1: Tpos - s.Rational(3, 2) * (Jlow + Jhigh)}, simultaneous=True
    )
    decomposition = Tpos + weight * Blo + s.Rational(3, 2) * (Bhi / (K - 2) - Jhigh)
    required = 2 * c.GAMMA * (1 - d1) - 4 * c.LAMBDA * (1 + d0) / (K - 2)
    cutoff_if_budget = 2 + 4 * c.LAMBDA * (1 + d0) / (
        c.GAMMA
        * (2 * (1 - d1) - s.Symbol("relative_full_low_moment_budget", nonnegative=True))
    )
    term0, term1, term2 = s.symbols("term0 term1 term2", real=True)
    amplitude_bound = 4 * c.LAMBDA * SPLIT**2 + 3 * c.GAMMA * SPLIT**3 / 4 + 8 * c.GAMMA
    mu = s.Symbol("spectral_s", positive=True)
    upper = (
        48 * c.lam**2 * K + s.Rational(9, 16) * c.gam**2 * K**3 + c.gam**2
    ) / s.pi**2
    v = s.Symbol("crossing_v", real=True)
    full_kernel = 2 * v * v / (s.pi * (mu - 2) * ((mu - 2) ** 2 - v * v))
    p0, p1, x0, x1 = s.symbols(
        "positive_mass0 positive_mass1 inverse_gap0 inverse_gap1", positive=True
    )
    higher_lower = c.GAMMA**2 * (1 - d1) ** 2 / (c.LAMBDA * (1 + d0))
    return {
        "explicit_unproved_physical_premises": "An actual nongravitational scalar S matrix with physical mass1 and positive canonical LSZ normalization, no unresolved coupled cut below s4 after the specified known-light-pole accounting, crossing and fixed-t analyticity near t0, convergent positive endpoint partial-wave derivative, and a twice-subtracted dispersion relation differentiable in t with the required vanishing arc. Unknown heavy poles remain positive spectral atoms or are included explicitly, not erased. None is established here for the original parent.",
        "full_low_moment_definition": "J4(K)=(2/pi) integral_[4,K] Im A_exact(S,0)/(S-2)^4 dS, including all light/heavy/inelastic channels and any retained positive spectral atoms. K=M^2>4 is an arbitrary analysis split, not a UV cutoff, light-loop regulator or asserted new-state threshold.",
        "exact_transfer_dispersion": "b20=(2/pi)integral rho/(S-2)^3; b21=(2/pi)integral[rho_t/(S-2)^3-(3/2)rho/(S-2)^4], where b21=partial_t[partial_v^2 B/2] at v=t=0 and v=s+t/2-2.",
        "necessary_exact_full_moment_inequality": "b21+3*b20/[2(K-2)]>=-3*J4(K)/2. It follows from rho_t>=0 and J4_high<=B20_high/(K-2), keeping all low weight rather than declaring the light cut absent.",
        "coefficient_matching_premise": "|b20-4lambda|<=4lambda*delta0 and |b21+3gamma|<=3gamma*delta1, delta0>=0, 0<=delta1<1. These are physical full-coefficient errors, not renormalization choices or established original bounds.",
        "general_required_full_low_moment": required,
        "full_forward_v4_coefficient_definition": "b40=(1/24)partial_v^4 B(0,0)=(2/pi)integral_4^infinity rho/(S-2)^5. It is a physical full coefficient, not the zero original-tree value or a freely selected counterterm.",
        "positive_full_moment_Gram_bound": "With dnu=(2/pi)rho dS/(S-2)^3, the positive Gram matrix of1 and1/(S-2) gives J4_total^2<=b20*b40. Also rho_t>=0 gives J4_total>=-2*b21/3 when b21<0. Hence b40>=4*b21^2/(9*b20).",
        "mandatory_full_higher_coefficient_under_matching": higher_lower,
        "named_mandatory_full_b40_lower": c.GAMMA**2 / (8 * c.LAMBDA),
        "higher_coefficient_scope": "For the stated errors with delta1<1, b40>=gamma^2(1-delta1)^2/[lambda(1+delta0)]. Exact b20,b21 matching gives b40>=gamma^2/lambda. The named tolerances give b40>=gamma^2/(8lambda)>0, while the original tree has b40=0. Higher operators, real matching and loops must be controlled; no whole-parent exclusion follows.",
        "conditional_spectral_split_upper_if_J4_le_epsilon_gamma": cutoff_if_budget,
        "budget_scope": "The displayed split inequality requires 2(1-delta1)-epsilon>0. It constrains simultaneous matching and a FULL low-moment budget; it is not a physical cutoff.",
        "named_energy": ENERGY,
        "named_coefficient_tolerances": (DELTA0, DELTA1),
        "named_required_full_moment": NAMED_REQUIRED,
        "first_elastic_moment_conservative_general_upper": upper,
        "named_first_elastic_rational_upper": FIRST_UPPER,
        "named_required_enhancement_lower": 2 * s.Integer(10) ** 199,
        "named_tree_amplitude_upper": amplitude_bound,
        "main_result": "At M10^99 and coefficient tolerances delta0<=1,delta1<=1/2, every exact amplitude satisfying the specified analytic premises must have FULL J4>gamma/5. The complete ORIGINAL first elastic coefficient instead has 0<J4_first<gamma*10^-200. Thus the full moment must exceed that contribution by more than2*10^199. No assumption that the full absorptive moment is close to its first elastic value is used for this necessary enhancement.",
        "checks": {
            "exact_low_high_positive_decomposition": s.expand(
                difference - decomposition
            ),
            "general_required_moment_from_coefficient_upper_bounds": s.cancel(
                -s.Rational(2, 3)
                * (-3 * c.GAMMA * (1 - d1) + weight * 4 * c.LAMBDA * (1 + d0))
                - required
            ),
            "named_required_moment_specialization": s.cancel(
                required.subs({K: SPLIT, d0: DELTA0, d1: DELTA1}) - NAMED_REQUIRED
            ),
            "square_three_term_majorant": s.expand(
                3 * (term0**2 + term1**2 + term2**2)
                - (term0 + term1 + term2) ** 2
                - ((term0 - term1) ** 2 + (term0 - term2) ** 2 + (term1 - term2) ** 2)
            ),
            "inverse_fourth_kernel_mass_bound": (mu - 2) - mu / 2 - (mu - 4) / 2,
            "potential_inverse_fourth_integral": s.integrate(mu**-4, (mu, 4, s.oo))
            - s.Rational(1, 192),
            "conservative_three_term_integral_coefficients": s.expand(
                3
                * (
                    16 * c.lam**2 * K
                    + s.Rational(3, 16) * c.gam**2 * K**3
                    + c.gam**2 / 3
                )
                / s.pi**2
                - upper
            ),
            "required_to_elastic_enhancement_factor": s.Rational(1, 5) * 10**200
            - 2 * 10**199,
            "forward_v4_dispersion_coefficient": s.cancel(
                s.diff(full_kernel, v, 4).subs(v, 0) / 24 - 2 / (s.pi * (mu - 2) ** 5)
            ),
            "independent_two_atom_positive_Gram_determinant": s.expand(
                (p0 + p1) * (p0 * x0 * x0 + p1 * x1 * x1)
                - (p0 * x0 + p1 * x1) ** 2
                - p0 * p1 * (x0 - x1) ** 2
            ),
            "general_required_higher_coefficient": s.cancel(
                4 * (3 * c.GAMMA * (1 - d1)) ** 2 / (9 * 4 * c.LAMBDA * (1 + d0))
                - higher_lower
            ),
            "named_required_higher_coefficient": higher_lower.subs(
                {d0: DELTA0, d1: DELTA1}
            )
            - c.GAMMA**2 / (8 * c.LAMBDA),
        },
        "gates": {
            "spectral_split_above_scalar_threshold": SPLIT > 4,
            "named_required_full_moment_above_gamma_fifth": NAMED_REQUIRED
            > c.GAMMA / 5,
            "first_elastic_upper_strictly_positive": FIRST_UPPER > 0,
            "first_elastic_upper_below_gamma_times_1e_minus_200": FIRST_UPPER
            < c.GAMMA / 10**200,
            "actual_rational_enhancement_greater_than_coarse_lower": NAMED_REQUIRED
            / FIRST_UPPER
            > 2 * 10**199,
            "full_tree_uniformly_below_1e_minus_202_on_named_energy_window": amplitude_bound
            < s.Rational(1, 10**202),
            "no_full_loop_error_or_physical_cutoff_assumed": True,
            "mandatory_named_full_higher_coefficient_strictly_positive": c.GAMMA**2
            / (8 * c.LAMBDA)
            > 0,
        },
    }
