"""Complete first light two-point coefficient in an explicit separate-model scheme."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_scalar_tree_matching import matching, model

x = s.Symbol("Feynman_x", real=True)
z = s.Symbol("complex_invariant_s")
M = s.Symbol("positive_heavy_mass", positive=True)
m2 = s.Symbol("heavy_mass_squared", positive=True)
g2 = s.Symbol("positive_coupling_squared", positive=True)
F = (1 - x) ** 2 + x * m2
ALPHA = x * (1 - x) / F
U = ALPHA * (z - 1)
REMAINDER_INTEGRAND = -s.log(1 - U) - U
RADIUS = matching.ENERGY**2
R0 = s.Rational(32, 625)
EPSILON_BOUND = model.G2 * R0 / (864 * model.MASS2 * (1 - R0))
SLOPE_BOUND = model.G2 / (288 * model.MASS2)
MASS_BOUND = (-model.CONTACT + 924 * model.G2) / 288


@cache
def data():
    a, u = s.symbols("positive_alpha log_argument_shift", positive=True)
    b0 = -s.log(x * m2 + (1 - x) - x * (1 - x) * z)
    rem = -s.log(1 - a * (z - 1)) - a * (z - 1)
    polynomial = x * M * M + (1 - x) - x * (1 - x) * z
    q, ds = s.symbols("relative_inverse_correction s_minus_one")
    checks = {
        "exact_mass_one_anchor_denominator": s.expand(
            F - (x * m2 + (1 - x) - x * (1 - x))
        ),
        "anchor_denominator_lower_margin": s.expand(F - 1 - x * (m2 - 2 + x)),
        "anchor_denominator_upper_margin": s.expand(m2 - F - (1 - x) * (m2 - 1 + x)),
        "pointwise_alpha_majorant": s.cancel(
            (1 - x) / m2 - ALPHA - (1 - x) ** 3 / (m2 * F)
        ),
        "first_bubble_anchor_derivative": s.cancel(s.diff(b0, z).subs(z, 1) - ALPHA),
        "second_bubble_anchor_derivative": s.cancel(
            s.diff(b0, z, 2).subs(z, 1) - ALPHA**2
        ),
        "subtracted_integrand_anchor_value": rem.subs(z, 1),
        "subtracted_integrand_anchor_slope": s.diff(rem, z).subs(z, 1),
        "subtracted_integrand_second_coefficient": s.diff(rem, z, 2).subs(z, 1) / 2
        - a * a / 2,
        "log_remainder_derivative": s.cancel(
            s.diff(-s.log(1 - u) - u, u) - u / (1 - u)
        ),
        "alpha_first_integrated_majorant": s.integrate(1 - x, (x, 0, 1))
        - s.Rational(1, 2),
        "alpha_second_integrated_majorant": s.integrate((1 - x) ** 2, (x, 0, 1))
        - s.Rational(1, 3),
        "full_loop_remainder_bound_coefficient": g2
        / (16 * s.pi**2)
        * s.Rational(1, 2)
        * s.Rational(1, 3)
        - g2 / (96 * s.pi**2),
        "threshold_positive_square": s.expand(
            polynomial.subs(z, (M + 1) ** 2) - ((M + 1) * x - 1) ** 2
        ),
        "pseudothreshold_positive_square": s.expand(
            polynomial.subs(z, (M - 1) ** 2) - ((M - 1) * x + 1) ** 2
        ),
        "inverse_plus_self_energy_sign": s.cancel(
            1 / (ds * (1 + q)) - 1 / ds + q / (ds * (1 + q))
        ),
        "propagator_ratio_difference": s.cancel(1 / (1 + q) - 1 + q / (1 + q)),
        "named_radius_matches_complete_tree_window": RADIUS - matching.ENERGY**2,
        "rational_pi_lower_bound_majorant": s.Rational(1, 96 * 9) - s.Rational(1, 864),
        "logarithm_mass_majorant_constant": 198 * s.Rational(7, 3) - 462,
    }
    return {
        "scheme": "V2S-T1 only: dimensional regularization in4-2epsilon, MSbar at mu=1 for UV poles, the heavy one-point-zero condition, and finite local light mass/kinetic subtractions at s=1. Renormalized tree parameters retain S233 values. No finite four-point matching condition is imposed here and no original affine finite prescription is changed.",
        "loop_sign": "Define i Pi as the one-particle-irreducible insertion. Dyson gives i/[s-1+Pi(s)]. The mixed bubble has factor g^2/(16pi^2); the quartic tadpole has -C A0(1)/(32pi^2). A0_MSbar(1)=1 at mu1. The one-point-reducible heavy tadpole is cancelled by the separately displayed one-point condition, not silently discarded.",
        "complete_MSbar_coefficient": "Pi_MS(s)=-C/(32pi^2)+g^2 B0_MS(s;1,M_H^2)/(16pi^2), B0_MS=-integral_0^1 Log[x M_H^2+(1-x)-x(1-x)s-i0]dx.",
        "anchor_F": F,
        "anchor_alpha": ALPHA,
        "complete_on_shell_subtracted_integrand": REMAINDER_INTEGRAND,
        "subtracted_coefficient": "Pi_OS(s)=g^2/(16pi^2) integral[-Log(1-alpha(x)(s-1))-alpha(x)(s-1)]dx. Exactly Pi_OS(1)=Pi_OS'(1)=0. These are renormalization conditions at this loop order, not a proof of exact physical mass or LSZ.",
        "first_sheet_domain": "The mixed one-loop bubble has its physical first-sheet cut at s>=(M_H+1)^2. The pseudothreshold(M_H-1)^2 has a strictly positive Feynman denominator and is not an extra first-sheet cut. This does not locate all cuts of the full quantum model; lower multiparticle light cuts can enter at higher loop orders.",
        "uniform_complex_bound": "For abs(s-1)<=R and r=R/M_H^2<1, abs(Pi_OS)<=g^2 abs(s-1)^2/[96pi^2 M_H^4(1-r)]. The quotient Pi_OS/(s-1) is analytic with value0 at the anchor. This uses alpha<=(1-x)/M_H^2 and the complete complex logarithm remainder, not a real-axis sample.",
        "actual_radius": RADIUS,
        "actual_relative_inverse_upper": EPSILON_BOUND,
        "actual_propagator_comparison": "On the full closed complex disk abs(s-1)<=10^196, abs(Pi_OS/(s-1))<=epsilon_bound and abs((s-1)/(s-1+Pi_OS)-1)<10^-209, with the ratio continued at s1. The reciprocal of this one-loop-truncated inverse has exactly the anchor simple pole with unit residue and no other poles in that disk. Algebraically inverting a truncated inverse does not bound omitted loops or establish the exact propagator.",
        "actual_first_coefficient_mass_and_slope": "At mu1, abs(Pi_MS(1))<10^-7 and0<Pi_MS'(1)<10^-207. The first bound uses1<=F<=M_H^2<10^198, ln10<7/3 and pi>3. The formal first pole shift before the finite on-shell subtraction is -Pi_MS(1). These statements do not control higher-order shifts.",
        "mass_coefficient_upper": MASS_BOUND,
        "slope_coefficient_upper": SLOPE_BOUND,
        "checks": checks,
        "gates": {
            "actual_anchor_denominator_monotone": model.MASS2 > 2,
            "actual_full_disk_below_ratio_bound": 0 < RADIUS / model.MASS2 < R0 < 1,
            "actual_inverse_quotient_small": 0 < EPSILON_BOUND < s.Rational(1, 10**209),
            "actual_propagator_ratio_error_small": EPSILON_BOUND / (1 - EPSILON_BOUND)
            < s.Rational(1, 10**209),
            "actual_first_loop_slope_small": 0 < SLOPE_BOUND < s.Rational(1, 10**207),
            "actual_first_loop_mass_shift_small": 0 < MASS_BOUND < s.Rational(1, 10**7),
            "actual_mass_logarithm_upper_domain": 1 < model.MASS2 < 10**198
            and sum(s.Rational(7, 3) ** k / s.factorial(k) for k in range(9)) > 10,
            "one_loop_condition_not_exact_LSZ_or_all_loop_error": True,
        },
    }
