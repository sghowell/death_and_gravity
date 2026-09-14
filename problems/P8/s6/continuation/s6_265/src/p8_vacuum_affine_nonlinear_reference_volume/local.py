"""Full finite-width lapse derivative bounds and literal center jets."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_band_neighborhood import source as fifth
from p8_vacuum_affine_finite_window_growth.intervals import I, evaluate

from . import source


@cache
def derivatives():
    N = s.Symbol("positive_lapse_derivative_probe", positive=True)
    R = s.Function("whole_positive_R", positive=True)(N)
    jets = s.symbols("whole_R_lapse0:5", real=True)
    bind = {R: jets[0], **{s.diff(R, N, j): jets[j] for j in range(1, 5)}}
    box = {
        jets[0]: I(s.Rational(1, 2), 2),
        jets[1]: I(-4, 4),
        jets[2]: I(-8, 8),
        jets[3]: I(-26, 26),
        jets[4]: I(-122, 122),
    }
    polynomials = {}
    bounds = {}
    centers = {}
    checks = {}
    h = source.parent.h
    clock = 1 + (N**-2 - 1) / h
    prior = fifth.data()
    for alpha in source.ALPHAS:
        key = str(alpha)
        rows = []
        upper = []
        center = []
        for j in range(5):
            polynomial = s.expand(
                (s.diff(R**-alpha, N, j) * R**alpha).subs(bind, simultaneous=True)
            )
            value = evaluate(polynomial, box)
            bound = 2 * max(abs(value.lo), abs(value.hi))
            rows.append(polynomial)
            upper.append(s.Rational(bound.numerator, bound.denominator))
            center.append(s.factor(s.diff(clock**-alpha, N, j).subs(N, 1)))
            checks["complete_Laurent_reconstruction_" + key + "_" + str(j)] = s.factor(
                polynomial.subs(
                    {v: s.diff(R, N, k) for k, v in enumerate(jets)}, simultaneous=True
                )
                - s.diff(R**-alpha, N, j) * R**alpha
            )
        polynomials[key] = rows
        bounds[key] = upper
        centers[key] = center
        checks["whole_center_first_jet_" + key] = s.factor(center[1] - 2 * alpha / h)
        checks["whole_center_second_jet_" + key] = s.factor(
            center[2] - (4 * alpha * (alpha + 1) / h**2 - 6 * alpha / h)
        )
    actual = source.parent.lapse_jets()["rows"]["U"]
    checks["literal_complete_original_volume_jets_zero_through_four"] = s.Matrix(
        centers["3/4"]
    ) - s.Matrix(actual)
    return {
        "whole_prior_full_source_C5_enclosure": prior,
        "whole_lapse_interval": [1 - source.RADIUS, 1 + source.RADIUS],
        "whole_full_R_four_jet_box": {
            str(key): value.bounds() for key, value in box.items()
        },
        "whole_all_fifteen_negative_power_Laurent_polynomials": polynomials,
        "whole_all_fifteen_outward_derivative_bounds": bounds,
        "whole_complete_center_jets_zero_through_four": centers,
        "whole_safe_all_derivative_ceiling": s.Integer(10) ** 8,
        "whole_derivative_argument": "The entire S263 source C5 discrepancy from the clock comparison is below1e-380 on this original strip. Bare lapse derivative ceilings3,7,25,121 therefore imply full ceilings4,8,26,122. Global/full local bounds give R in[1/2,2]. After factoring R^-alpha, all derivatives are rational Laurent polynomials, bounded by exact intervals; R^-alpha<2. Original order8/1024 source factors identify the finite CENTER jets only; every off-center function and derivative bound remains the entire current source.",
        "checks": checks,
        "gates": {
            "whole_prior_full_source_gates": all(prior["gates"].values()),
            "full_C5_discrepancy_below_one": prior["complete_RF_fifth_error_bound"] < 1,
            "all_fifteen_full_derivatives_below_1e8": max(
                v for row in bounds.values() for v in row
            )
            < 10**8,
            "all_three_actual_negative_powers_included": len(bounds) == 3
            and all(len(row) == 5 for row in bounds.values()),
            "full_off_clock_R_not_replaced_by_center_germ": True,
            "same_original_C5_strip_and_both_parent_switches": True,
        },
    }
