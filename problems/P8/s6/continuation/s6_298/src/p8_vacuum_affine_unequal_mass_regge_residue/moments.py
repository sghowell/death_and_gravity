"""Exact full spectral moment, finite tail and original heavy-quarter bound."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_matter_graviton_vertex import vertex

from . import source

MU, N, G = source.MU, source.N, source.G
Z = s.Symbol("massive_parameter", real=True)
DELTA = s.Rational(1, 10**98)


def anchor(z=Z, mass=MU, heavy=N):
    z, mass, heavy = map(s.sympify, (z, mass, heavy))
    return mass * (1 - z) ** 2 + heavy * z


def moment_integrand(species, z=Z, mass=MU, heavy=N):
    from .waves import require_species

    require_species(species)
    z, mass, heavy = map(s.sympify, (z, mass, heavy))
    weight = z * z * (1 - z) ** 3 if species == "light" else z**3 * (1 - z) ** 2
    return weight / anchor(z, mass, heavy) ** 2


def F1_slope(species, mass=MU, heavy=N, cubic=G):
    return (
        s.sympify(cubic) ** 2
        * s.Integral(moment_integrand(species, Z, mass, heavy), (Z, 0, 1))
        / (96 * s.pi**2)
    )


def moment_bounds(heavy=N, cubic=G):
    heavy, cubic = map(s.sympify, (heavy, cubic))
    return {
        "F1prime_lower": cubic**2 / (2880 * s.pi**2 * heavy**2),
        "F1prime_upper": cubic**2 / (288 * s.pi**2 * heavy**2),
        "unsubtracted_measure_upper": cubic**2 / (32 * s.pi**2 * heavy),
    }


@cache
def data():
    x = Z
    F = anchor()
    a, y, delta = s.symbols("mass_ratio y split_delta", positive=True)
    checks = {
        "both_full_F1_parameter_moments": s.factor(
            moment_integrand("light")
            + moment_integrand("heavy")
            - vertex.self_second_integrand(Z, MU, N)
        ),
        "both_unsubtracted_zero_transfer_measures": s.factor(
            x * x * (1 - x) / F + x * (1 - x) ** 2 / F - x * (1 - x) / F
        ),
        "unsubtracted_measure_pointwise_majorant": s.factor(
            (1 - x) / N - x * (1 - x) / F - MU * (1 - x) ** 3 / (N * F)
        ),
        "whole_Feynman_anchor_upper_n": s.factor(N - F - (1 - x) * (N - MU * (1 - x))),
        "full_lower_weight_integral": s.integrate(x * x * (1 - x) ** 2, (x, 0, 1))
        - s.Rational(1, 30),
        "unsubtracted_upper_weight_integral": s.integrate(1 - x, (x, 0, 1))
        - s.Rational(1, 2),
        "light_large_hierarchy_weight": s.integrate((1 - x) ** 3, (x, 0, 1))
        - s.Rational(1, 4),
        "heavy_large_hierarchy_weight": s.integrate(x * (1 - x) ** 2, (x, 0, 1))
        - s.Rational(1, 12),
        "sum_large_hierarchy_weights": s.Rational(1, 4)
        + s.Rational(1, 12)
        - s.Rational(1, 3),
        "positive_denominator_defect_bound": s.factor(
            2 * y - (1 - (1 + y) ** -2) - y * y * (3 + 2 * y) / (1 + y) ** 2
        ),
        "heavy_defect_full_integral": 2 * s.integrate((1 - x) ** 4, (x, 0, 1))
        - s.Rational(2, 5),
        "total_split_majorant": delta + 2 * delta / 5 - 7 * delta / 5,
        "relative_total_slope_defect": 7 * delta / 5 / (s.Rational(1, 3))
        - 21 * delta / 5,
        "heavy_fraction_upper_defect": s.factor(
            s.Rational(1, 12) / (s.Rational(1, 3) - 7 * delta / 5)
            - s.Rational(1, 4)
            - 21 * delta / (4 * (5 - 21 * delta))
        ),
        "heavy_fraction_lower_defect": s.Rational(1, 4)
        - (s.Rational(1, 12) - 2 * a / 5) * 3
        - 6 * a / 5,
    }
    sigma, U, weight = s.symbols(
        "spectral_transfer cut_cap positive_weight", positive=True
    )
    checks["positive_measure_tail_majorant"] = s.factor(
        weight / (U * sigma) - weight / sigma**2 - weight * (sigma - U) / (U * sigma**2)
    )
    q = s.Symbol("light_threshold_offset", positive=True)
    checks["light_low_window_density_upper"] = s.factor(
        G**2
        * q**2
        / (240 * s.pi**2 * N**3 * (q + 4 * MU) ** 2)
        / (G**2 / (1440 * s.pi**2 * N**2))
        - 6 * q**2 / (N * (q + 4 * MU) ** 2)
    )
    checks["light_low_window_integrand_below_one"] = s.factor(
        1 - q * q / (q + 4 * MU) ** 2 - 8 * MU * (q + 2 * MU) / (q + 4 * MU) ** 2
    )
    original_a = s.S.One / source.HEAVY_MASS2
    bounds = moment_bounds(source.HEAVY_MASS2, source.CUBIC)
    return {
        "whole_light_F1_slope": F1_slope("light"),
        "whole_heavy_F1_slope": F1_slope("heavy"),
        "whole_total_F1_slope": F1_slope("light") + F1_slope("heavy"),
        "whole_spectral_slope_bounds": moment_bounds(),
        "whole_original_slope_bounds": bounds,
        "positive_tail_bound": "For each generated massive F1 Stieltjes measure, integral_U^infinity rho/(pi*T^2) <= Pi_prime(mu)/U. Both graphs give Pi_prime=g^2/(16pi^2) integral z(1-z)/F dz <=g^2/(32pi^2*n). This is the full known formal graph tail, not a physical UV spectral completion.",
        "light_window_bound": "For 4mu<U<=4n only the light cut is open. Its twice-F1 slope weight is <=g^2*(U-4mu)/(240pi^2*n^3), hence its fraction of the full twice-F1 slope is <6U/n. The positive Bose P2 series and a^2-b^2=n(n+q) establish the pointwise bound; no threshold-dominance assumption is used.",
        "whole_original_heavy_fraction_bound": {
            "split_delta": DELTA,
            "mass_ratio": original_a,
            "relative_total_slope_error_upper": 21 * DELTA / 5,
            "absolute_heavy_fraction_error_upper": 2 * DELTA,
            "reference_heavy_fraction": s.Rational(1, 4),
            "low_window_U16_fraction_upper": s.Rational(96) / source.HEAVY_MASS2,
        },
        "hierarchy_proof": "Put a=mu/n and split z at delta with a<=delta^2. The scaled total moment J is between1/3-7delta/5 and1/3: on z<delta its defect is <=delta, and elsewhere the defect <=2a(1-z)^4/z integrates to <=2a/(5delta). The scaled heavy moment H is between1/12-2a/5 and1/12 without a split. Their positive ratio lies within2delta of1/4 for delta<1/10. At the original mass-one hierarchy delta=10^-98, so the heavy cut is a quarter within2*10^-98 and the total leading hierarchy formula has relative error below10^-97.",
        "not_IR_enhanced_equal_mass_formula": "At fixed unequal n the whole slope scales as g^2/(288pi^2*n^2), with the actual proved relative error. The light-particle massless enhancement of the equal-mass self-cubic example cannot be transplanted to this source. A cut restricted to4mu<T<=16 captures less than10^-195 of the complete spin2 slope at the original hierarchy.",
        "checks": checks,
        "gates": {
            "original_mass_ratio_below_split_squared": 0 < original_a < DELTA**2,
            "positive_total_lower_bound": s.Rational(1, 3) - 7 * DELTA / 5 > 0,
            "strict_heavy_upper_fraction_error": 21 * DELTA / (4 * (5 - 21 * DELTA))
            < 2 * DELTA,
            "strict_heavy_lower_fraction_error": 6 * original_a / 5 < 2 * DELTA,
            "strict_original_relative_slope_error": 21 * DELTA / 5
            < s.Rational(1, 10**97),
            "strict_original_low_window_fraction": 0
            < s.Rational(96) / source.HEAVY_MASS2
            < s.Rational(1, 10**195),
            "full_positive_stieltjes_measure_not_full_UV_measure": True,
        },
    }
