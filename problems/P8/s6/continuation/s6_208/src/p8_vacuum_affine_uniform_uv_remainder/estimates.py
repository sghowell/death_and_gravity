"""Complete far/near UV-subtracted row in the original six-spatial-derivative norm."""

from functools import cache

import sympy as s
from p8_vacuum_affine_reference_spatial_remainder import real
from p8_vacuum_affine_retarded_reference_boundary import majorants
from p8_vacuum_affine_spatial_symbol import extraction

from . import domain

FINITE = 2 * s.Integer(10) ** 48
LOW = 2 * s.Integer(10) ** 38


def row_bound(endpoint):
    j = extraction.require_order(endpoint)
    return 9 * domain.PAIR**2 * sum(majorants.recurrence()[j]) * domain.RADIUS ** (-j)


def jet_bound(endpoint, source_jet):
    j = extraction.require_order(endpoint)
    if (
        isinstance(source_jet, bool)
        or not isinstance(source_jet, (int, s.Integer))
        or not 0 <= source_jet <= j
    ):
        raise ValueError("Exact source time jet0,...,endpoint required")
    r = int(source_jet)
    return (
        9 * domain.PAIR**2 * majorants.recurrence()[j][r] * domain.RADIUS ** (-(j - r))
    )


def spatial_weight(power):
    if (
        isinstance(power, bool)
        or not isinstance(power, (int, s.Integer))
        or not 0 <= power <= 6
    ):
        raise ValueError("Exact spatial degree0,...,6 required")
    return 8 * domain.MASS ** int(power)


@cache
def far_terms():
    return tuple((4 - j, row_bound(j) / (18 * domain.EPS ** (4 - j))) for j in range(5))


@cache
def near_terms():
    e, m = domain.EPS, domain.MASS
    rows = [(4, 4 * real.REAL / (9 * e**4), "complete_real_endpoint_row")]
    for j, d in extraction.slots():
        q = j + d
        if q < 4:
            rows.append(
                (
                    4 - j,
                    row_bound(j) * 2 ** (4 - j - d) / (18 * (4 - j - d) * e ** (4 - j)),
                    f"endpoint_{j}_degree_{d}",
                )
            )
        else:
            rows.append(
                (
                    5 - j,
                    2 * row_bound(j) / (18 * m * e ** (5 - j)),
                    f"endpoint_{j}_degree_{d}_log",
                )
            )
    return tuple(rows)


@cache
def constants():
    far = sum(spatial_weight(q) * c for q, c in far_terms())
    near = sum(spatial_weight(q) * c for q, c, _ in near_terms())
    return {
        "all_five_normalized_endpoint_row_bounds": tuple(
            row_bound(j) for j in range(5)
        ),
        "complete_weighted_far_bound": far,
        "complete_weighted_near_bound": near,
        "unexpanded_low_band_bound": LOW,
        "complete_all_momentum_UV_subtracted_endpoint_bound": LOW + far + near,
    }


@cache
def data():
    r, L, m, P = s.symbols("r L m P", positive=True)
    checks = {}
    for q in range(4):
        checks[f"complete_near_radial_power_{q}"] = s.integrate(
            r ** (3 - q), (r, m, L)
        ) - (L ** (4 - q) - m ** (4 - q)) / (4 - q)
    checks.update(
        {
            "complete_near_logarithmic_radial_power": s.integrate(1 / r, (r, m, L))
            - s.log(L / m),
            "complete_far_fourth_power_radial_integral": s.integrate(
                r**-2, (r, L, s.oo)
            )
            - 1 / L,
            "all_fifteen_UV_slots_retained": len(extraction.slots()) - 15,
            "all_thirty_five_source_time_jets_retained": sum(
                j + 1 for j, d in extraction.slots()
            )
            - 35,
            "maximum_finite_spatial_power": max(q for q, c, label in near_terms()) - 5,
            "complete_spatial_weight_Cauchy_identity": s.expand(
                2 * (m * m + P * P) - (m + P) ** 2 - (m - P) ** 2
            ),
            "complete_mass_weight_monotonicity": s.expand(
                m * m * (1 + P * P) - (m * m + P * P) - (m * m - 1) * P * P
            ),
            "sixth_degree_weight_constant": 2**3 - 8,
        }
    )
    return {
        "coefficient_bounds": "Each normalized source time-jet coefficient has bound9 C^2 c_jr R^(-(j-r)), C=1000. Summing all source jets gives B_j=9 C^2 sum(c_jr)R^-j. The inverse-radius degree-d coefficient is bounded by B_j rho^-d. All nine pairs and every source jet remain.",
        "far": "For r>=L(P)=2(m+|P|)/epsilon, Cauchy gives the complete endpoint-j remainder<=2B_j epsilon^(j-5)(m+|P|)^(5-j)r^-4. Integrating its actual far threshold yields B_j(m+|P|)^(4-j)/(18epsilon^(4-j)), using pi^2>9.",
        "near": "For m<=r<L(P), bound the actual complete raw row by4e27nu<8e27r. Retain every UV coefficient separately. For j+d<4 integrate the complete power; for j+d=4 retain log(L/m), bounded byL/m. The complete near powers are at most5. No odd endpoint or finite oversubtraction is deleted.",
        "low": "The unchanged low bandr<m remains unexpanded and is bounded by2e38 from S203. The improved normalized domain constants are not applied there.",
        "weight": "For every integer q0,...,6, (m+|P|)^q<=8m^q(1+|P|^2)^3 because m>=1. Plancherel and the full time-jet row give |Qnew|,|Qnew,K|<2e48||D||L2 X46[Gamma], with X46^2=sum_r0..4||(1-Delta)^3 partial_t^r Gamma||L2^2.",
        "constants": constants(),
        "checks": checks,
        "gates": {
            "complete_uniform_finite_display": bool(
                constants()["complete_all_momentum_UV_subtracted_endpoint_bound"]
                < FINITE
            ),
            "all_five_far_rows": len(far_terms()) == 5,
            "complete_raw_and_fifteen_near_coefficients": len(near_terms()) == 16,
            "all_source_jets_bounded": all(
                sum(jet_bound(j, r) for r in range(j + 1)) <= row_bound(j)
                for j in range(5)
            ),
            "unchanged_low_band_anchor": bool(
                real.constants()["unexpanded_low_band_integral_upper"] < LOW
            ),
            "at_most_five_finite_spatial_powers": max(q for q, c, label in near_terms())
            <= 5,
            "fixed_mass_weight_at_least_one": bool(domain.MASS >= 1),
        },
    }
