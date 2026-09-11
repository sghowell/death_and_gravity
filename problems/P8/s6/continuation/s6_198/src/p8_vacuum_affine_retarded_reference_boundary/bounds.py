"""Full spatial continuum bounds for the sixth bulk and fifth endpoint."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import reference
from p8_vacuum_affine_spatial_current.bounds import require_band
from p8_vector_state.comparison import AMAX

from . import majorants

DISPLAY = s.Integer(10) ** 48
TAIL = s.Integer(10) ** 52


@cache
def constants():
    m = reference.MASS
    A = AMAX
    r = majorants.RADIUS
    C = majorants.PAIR
    c6 = sum(majorants.recurrence()[6])
    c5 = sum(majorants.recurrence()[5])
    bulk = 9 * C * C * 64 * c6 / r**6
    endpoint = 9 * C * C * 64 * c5 / r**5
    J4 = A**3 / (24 * m)
    return {
        "sixth_bulk_pair_numerator": bulk,
        "fifth_endpoint_pair_numerator": endpoint,
        "sixth_bulk_full_integral": bulk * J4 / 4,
        "fifth_endpoint_full_integral": endpoint * J4 / 4,
        "complete_finite_reference_integral": (bulk + endpoint) * J4 / 4,
        "complete_removed_two_leg_tail_numerator": (bulk + endpoint) * A**4 / 36,
    }


def tail_bound(K):
    return TAIL / require_band(K)


@cache
def data():
    nu, mu = s.symbols("nu mu", positive=True)
    x, K = s.symbols("x K", positive=True)
    c = constants()
    checks = {
        "AM_GM_pair_denominator": s.expand(
            (nu + mu) ** 2 - 4 * nu * mu - (nu - mu) ** 2
        ),
        "sixth_denominator_equal_momentum": (nu * mu / (nu + mu) ** 6).subs(mu, nu)
        - 1 / (64 * nu**4),
        "fifth_denominator_equal_momentum": (nu * mu / (nu + mu) ** 5).subs(mu, nu)
        - 1 / (32 * nu**3),
        "complete_reference_sum": c["complete_finite_reference_integral"]
        - c["sixth_bulk_full_integral"]
        - c["fifth_endpoint_full_integral"],
        "infinite_tail_power": s.integrate(x**-2, (x, K, s.oo)) - 1 / K,
        "removed_tail_at_proof_partition": tail_bound(10**16) - 10**36,
    }
    return {
        "norm": "S6[Gamma]^2=sum_(j=0)^6 ||partial_t^j Gamma||L2(dt dx;F)^2. The detector uses ||D||L2(dt dx;F). No spatial Fourier cutoff is imposed.",
        "sixth_bulk": "The complete nine-pair retarded bulk is bounded by its recorded numerator times nu mu/(nu+mu)^6. Cauchy on the unit time triangle uses S6[Gamma] and ||D||L2 without differentiating the retarded step.",
        "fifth_endpoint": "The j5 equal-time endpoint has one extra g outside L5, so it has the SAME sixth inverse-frequency power. Its pointwise time product is bounded by L2 Cauchy with source jets only through5.",
        "all_momentum": "At fixed external P, nu mu/(nu+mu)^6<=1/(4nu^4), hence internal integral<=J4/4 uniformly in P. Plancherel removes the external integral without a spatial derivative weight for this finite reference part.",
        "complete_bound": "The sum F=endpoint_j5+bulk_6 has |F|<1e48 ||D||L2 S6[Gamma]. Its common both-leg regulator error is<1e52 ||D||L2 S6[Gamma]/K for K>=1000. It is zero when either test is zero.",
        "tail": "The removed union max(|k|,|l|)>K gives J4_tail(K)/2<=Amax^4/(4pi^2 K)<Amax^4/(36K), uniformly in all external momentum.",
        "fifth_step_boundary": "The corresponding fifth-bulk absolute weight at P0 has radial behavior1/(32k), so this absolute majorant is not integrable. This does not assert divergence of the exact fifth-step current after possible cancellations.",
        "constants": c,
        "checks": checks,
        "gates": {
            "sixth_bulk_strict_display": c["sixth_bulk_full_integral"] < 10**48,
            "fifth_endpoint_strict_display": c["fifth_endpoint_full_integral"] < 10**43,
            "complete_finite_reference_strict_display": c[
                "complete_finite_reference_integral"
            ]
            < DISPLAY,
            "complete_finite_reference_tail_display": c[
                "complete_removed_two_leg_tail_numerator"
            ]
            < TAIL,
            "sixth_is_integrable_with_three_dimensional_measure": s.integrate(
                x**-2, (x, 1, s.oo)
            )
            == 1,
            "fifth_absolute_majorant_not_integrable": s.integrate(1 / x, (x, 1, s.oo))
            == s.oo,
        },
    }
