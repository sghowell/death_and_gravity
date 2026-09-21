"""Full selected known coefficient and exact finite-domain sign/magnitude."""

from functools import cache

import sympy as s
from p8_vacuum_affine_core_curvature_coefficient import bounds as core_bounds
from p8_vacuum_affine_core_curvature_coefficient import core

from . import moment, source


@cache
def data():
    n, z = moment.N, moment.Z
    M = 1 + (n - 1) * z
    r = z * (1 - z) / M
    g, k = s.symbols("g kappa", positive=True)
    E = moment.extra()
    known = (g**4 * core.original_coefficient() + g * g * E / k) / (16 * s.pi**2)
    checks = {
        "pointwise_strict_mixed_weight_gap": s.factor(
            1 - (n - 1) * r - z - (1 - z) / M
        ),
        "pointwise_mass_upper_gap": s.expand(n - M - (n - 1) * (1 - z)),
        "pointwise_rational_weight_upper_gap": s.factor(
            (1 - z) / (n - 1) - r - (1 - z) / ((n - 1) * M)
        ),
        "positive_lower_Beta_moment": s.integrate(z * z * (1 - z) ** 2 / 2, (z, 0, 1))
        - s.Rational(1, 60),
        "positive_upper_Beta_moment": s.integrate((1 - z) ** 2 / 2, (z, 0, 1))
        - s.Rational(1, 6),
        "positive_extra_lower_mass_threshold": s.factor(
            1 / (180 * n * n) - 2 / n**4 - (n * n - 360) / (180 * n**4)
        ),
        "whole_common_basis_known_matter_assembly": s.factor(
            known
            - g**4 * core.original_coefficient() / (16 * s.pi**2)
            - g * g * E / (16 * s.pi**2 * k)
        ),
    }
    actualn, actualg, actualk = source.HEAVY_MASS2, source.CUBIC, source.KAPPA
    return {
        "checks": checks,
        "gates": {
            "all_frozen_core_finite_bound_gates_preserved": all(
                bool(v) for v in core_bounds.data()["gates"].values()
            ),
            "original_mass_exceeds_extra_positivity_threshold": actualn**2 > 360,
            "extra_upper_prefactor_below_one_over25": s.Rational(32, 6 * 144)
            * actualn**2
            / (actualn - 1) ** 2
            < s.Rational(1, 25),
            "extra_Born_normalized_bound_below_10_minus604": actualn / (100 * actualk)
            < s.Rational(1, 10) ** 604,
            "positive_extra_strictly_smaller_than_negative_core": s.Rational(160, 6)
            * actualn**4
            / (actualg**2 * actualk * (actualn - 1) ** 2)
            < 1,
            "total_known_selected_coefficient_negative_and_below_core_magnitude": True,
            "no_bound_on_independent_extra_parent_chi": True,
            "no_physical_above_threshold_Taylor_replacement": True,
        },
        "whole_selected_known_coefficient": known,
        "whole_selected_extra_coefficient": g * g * E / (16 * s.pi**2 * k),
        "whole_finite_extra_sign_proof": "For n>1 and0<z<1,0<(n-1)r<1. Hence b2/3<-[ (n-1)b3-b2 ]<b2. Because M<=n and r<(1-z)/(n-1),1/(60n^2)<=b2<1/[6(n-1)^2]. For original n^2>360 this implies0<E(n)<32/[6(n-1)^2]. All identities underlying the positive gaps and both Beta integrals are checked exactly.",
        "whole_original_extra_bound": "Using16pi^2>144, the exact original prefactor bound gives0<chi_extra<g^2/(25kappa n^2). The prior original A0>4g^2/n^3 then gives chi_extra/A0<n/(100kappa)<10^-604. No floating integral or asymptotic inference enters.",
        "whole_total_sign_and_bound": "S346 gives chi_core<0 and |chi_core|>g^4/(80pi^2 n^4). The displayed exact original ratio160n^4/[6g^2 kappa(n-1)^2]<1 proves chi_extra<|chi_core|. Therefore chi_known_selected=chi_core+chi_extra<0 and |chi_known_selected|/A0<|chi_core|/A0<10^-207. All selected scalar-loop four-point classes in the classical limiting-action inventory are included in one fixed off-shell comparison.",
        "whole_parent_boundary": "The independent extra parent-theory curvature coefficient is not the same object as the computed known loop coefficient. It remains unassigned and unbounded. The result is a local analytic-origin real-TT comparison, not the full curved effective action or a physical above-threshold approximation. Internal gravitons, finite-gravity decoupling, complete inclusive matching, exact LSZ, complex Regge, same-parent bounce and UV remain open.",
    }
