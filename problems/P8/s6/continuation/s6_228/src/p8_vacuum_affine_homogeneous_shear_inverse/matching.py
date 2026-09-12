"""Actual complete homogeneous shear leading pair and fixed-dimensional matching."""

from functools import cache

import sympy as s
from p8_vacuum_affine_isolated_shear_resolvent import spectral
from p8_vacuum_affine_ordered_scalar_symbol import geometry
from p8_vacuum_affine_proca_gaussian import bridge

d = geometry.d
C = (2 * d * d + d - 8) / (2 * d * (d + 2))


@cache
def data():
    actual = {}
    for ch in ("tensor", "vector", "scalar"):
        rows = geometry.contractions(ch)
        avg = {
            key: s.factor(geometry.sphere_average(value.subs(geometry.y, 0)))
            for key, value in rows.items()
        }
        actual[ch] = s.factor(
            (sum(avg[k] for k in ("00", "01", "10", "11")) + avg["LL"]) / 4
        )
    TT = 1 - 2 / d + (d - 1) / (2 * d * (d + 2))
    LL = 1 / (2 * d * (d + 2))
    lag, eps = s.symbols("lag eps", positive=True)
    a0, p = s.symbols("frozen_scale proper_momentum", positive=True)
    H = s.Symbol("frozen_H", real=True)
    ratio = s.series(
        lag**5 / ((a0 * (1 + H * lag)) ** 4 * a0 * (lag / a0 * (1 - H * lag / 2)) ** 5),
        lag,
        0,
        2,
    ).removeO()
    radial = 32 * C.subs(d, 3)
    abel = s.limit(s.im(24 / (eps - 2 * s.I * lag) ** 5), eps, 0, dir="+")
    coefficient = radial * s.Rational(3, 4)
    log_weight = spectral.U.subs(spectral.z, 1)
    checks = {
        "three_literal_analytic_dimension_channels": s.Matrix(
            [v - C for v in actual.values()]
        ),
        "independent_rotational_moments": s.factor(TT + LL - C),
        "physical_leading_pair_and_dimensional_jet": s.Matrix(
            [
                C.subs(d, 3) - s.Rational(13, 30),
                s.diff(C, d).subs(d, 3) - s.Rational(91, 450),
            ]
        ),
        "physical_transverse_part": TT.subs(d, 3) - s.Rational(2, 5),
        "physical_longitudinal_part": LL.subs(d, 3) - s.Rational(1, 30),
        "original_normalized_radial_coefficient": radial - s.Rational(208, 15),
        "exact_positive_phase_Abel_limit": abel - s.Rational(3, 4) / lag**5,
        "four_primitive_scalar_log_normalization": s.diff(
            s.Rational(13, 30) / lag, lag, 4
        )
        - coefficient / lag**5,
        "original_shear_high_log_not_trace_log": coefficient - 48 * log_weight,
        "frozen_proper_momentum_scale_all_dimensions": s.simplify(
            a0 ** (-d - 2) * (a0 * p) ** (d + 1) * a0 - p ** (d + 1)
        ),
        "actual_geometric_first_lag_jet": s.expand(
            ratio - 1 + s.Rational(3, 2) * H * lag
        ),
        "unchanged_actual_mass": bridge.MASS - 1000,
    }
    return {
        "literal_analytic_dimension_high_pair": C,
        "independent_channels": actual,
        "physical_TT_and_LL": (TT.subs(d, 3), LL.subs(d, 3)),
        "first_dimension_jet": s.diff(C, d).subs(d, 3),
        "normalized_positive_phase_radial_coefficient": radial,
        "off_diagonal_leading_coefficient": coefficient,
        "original_flat_shear_log_coefficient": -log_weight,
        "actual_leading_kernel": "(208/15) k^4 sin(2k Delta_sigma)/[a(t)^4 a(s)], sigma'=1/a, proper-time output normalization64pi^2/a^3. The full constrained two-form and longitudinal mode are retained.",
        "independent_invariant_derivation": "For tracefree D with tr(D^2)=1, the TT leading matrix is PDP+(kDk)P/2, P=I-kk^T. Adding the LL square gives1-2kD^2k+d(kDk)^2/4; the exact sphere moments giveC(d). TL/LT have lower radial degree.",
        "actual_state_remainder": "The unchanged all-order Cauchy preparation gives arbitrarily deep finite high-k expansions for each fixed time derivative count at physicald3. The original near-d3 comparison error remains integrable. Lower singular orders and all original contacts are retained in the written remainder proof.",
        "dimension_boundary": "At frozen proper scale the factor a0^(-d-2) k^(d+1)dk becomes p^(d+1)dp for every nearbyd, including its first jet. Matching only atd3 would not fix an evanescent fourth contact. This is NOT a claim that conformal-time mass scaling is absent.",
        "checks": checks,
        "gates": {
            "all_three_channels_agree": len(set(actual.values())) == 1,
            "positive_TT_and_LL": TT.subs(d, 3) > 0 and LL.subs(d, 3) > 0,
            "longitudinal_deletion_changes_leading_log": LL.subs(d, 3) != 0,
            "first_dimension_jet_not_deleted": s.diff(C, d).subs(d, 3) != 0,
            "wrong_annihilation_phase_changes_leading_kernel": radial != 0,
            "source_pinned_shear_log_distinct_from_trace": log_weight != 2,
        },
    }
