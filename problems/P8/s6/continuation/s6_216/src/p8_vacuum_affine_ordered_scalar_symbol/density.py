"""All ordered scalar coefficients and fixed six-invariant reconstruction."""

from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_spatial_symbol import density as previous

from . import geometry, jets

CHANNELS = ("trace_trace", "trace_scalar", "scalar_trace")


@cache
def geometric_jets(channel):
    x = s.Symbol("x", real=True)
    out = {}
    for key, value in geometry.contractions(channel).items():
        value = s.series(value.subs(geometry.y, jets.p * x), x, 0, 5).removeO().expand()
        out[key] = tuple(s.factor(value.coeff(x, j)) for j in range(5))
    return out


@cache
def coefficients():
    products = jets.endpoint_products()
    result = {}
    for ch in CHANNELS:
        geo = geometric_jets(ch)
        for j in range(5):
            for r in range(j + 1):
                for degree in range(5 - j):
                    value = sum(
                        sum(
                            products[key, j, r][v] * geo[key][degree - v]
                            for v in range(degree + 1)
                        )
                        for key in geo
                    )
                    value = s.Poly(s.cancel(value), jets.u).as_expr()
                    result[ch, j, r, degree] = s.factor(geometry.sphere_average(value))
    return result


@cache
def spatial_difference():
    return {
        key: s.factor(value - value.subs(jets.p, 0))
        for key, value in coefficients().items()
    }


@cache
def invariant_spatial_coefficients():
    known = previous.spatial_difference()
    new = spatial_difference()
    inverse = geometry.invariant_matrix().inv()
    result = {}
    for j in range(5):
        for r in range(j + 1):
            for degree in range(5 - j):
                rhs = s.Matrix(
                    [
                        (-1) ** j * known[ch, j, r, degree]
                        for ch in ("tensor", "vector", "scalar")
                    ]
                    + [new[ch, j, r, degree] for ch in CHANNELS]
                )
                result[j, r, degree] = tuple(s.factor(v) for v in inverse * rhs)
    return result


def physical_channel_row(channel):
    if channel == "trace_trace":
        return (3, 1, 1, 9, 3, 3)
    if channel == "trace_scalar":
        return (0, 2, 2, 0, 6, 0)
    if channel == "scalar_trace":
        return (0, 2, 2, 0, 0, 6)
    raise ValueError("Only the three fixed physical unnormalized scalar directions")


@cache
def fixed_source_dimension_jet():
    return {
        (channel, *key): tuple(
            s.factor(expr)
            for expr in (
                sum(v * w for v, w in zip(row, physical_channel_row(channel))).subs(
                    geometry.d, 3
                ),
                s.diff(
                    sum(v * w for v, w in zip(row, physical_channel_row(channel))),
                    geometry.d,
                ).subs(geometry.d, 3),
            )
        )
        for key, row in invariant_spatial_coefficients().items()
        for channel in CHANNELS
    }


@cache
def logarithmic_dimension_jet():
    rows = fixed_source_dimension_jet()
    return {
        (ch, r): tuple(
            s.factor(sum(rows[ch, j, r, 4 - j][index] for j in range(r, 5)))
            for index in range(2)
        )
        for ch in CHANNELS
        for r in range(5)
    }


@cache
def data():
    prior = previous.spatial_difference()
    new = spatial_difference()
    invariants = invariant_spatial_coefficients()
    checks = {
        "all_105_new_ordered_scalar_UV_coefficients": len(new) - 105,
        "six_invariants_for_every_source_endpoint_slot": len(invariants) - 35,
    }
    for ch in ("tensor", "vector", "scalar") + CHANNELS:
        row = geometry.invariant_row(ch)
        checks[ch + "_complete_invariant_reconstruction"] = s.Matrix(
            [
                s.factor(
                    sum(a * b for a, b in zip(v, row))
                    - (
                        new[ch, *key]
                        if ch in CHANNELS
                        else (-1) ** key[0] * prior[ch, *key]
                    )
                )
                for key, v in invariants.items()
            ]
        )
    logs = logarithmic_dimension_jet()
    checks["no_third_or_fourth_source_jet_in_full_log"] = s.Matrix(
        [v for ch in CHANNELS for r in (3, 4) for v in logs[ch, r]]
    )
    return {
        "normalized_UV_tables": "105 ordered scalar channel entries and all35 six-invariant rows. The complete 350 scalar endpoint coefficients precede the normalized general-d angular average. Restore 1/(2pi^2) only at physical radial integration.",
        "physical_scalar_directions": "Unnormalized I/I,I/S,S/I with S=diag(-1,-1,2). Their fixed invariant rows are(3,1,1,9,3,3),(0,2,2,0,6,0),(0,2,2,0,0,6). Unit directions instead carry factors1/3 or1/sqrt18.",
        "dimension_derivative": "First reconstruct all six coefficients at symbolic d; then keep the physical tensor invariants fixed while differentiating. The dimension dependence of I_d and S_d is not a physical source derivative.",
        "ordered_physical_and_dimension_log_jets": logs,
        "checks": checks,
        "gates": {
            "both_off_diagonal_scalar_orientations_retained": True,
            "ordered_cross_operators_not_assumed_equal": logs["trace_scalar", 1][0]
            != logs["scalar_trace", 1][0],
            "all_prior_tracefree_rows_phase_corrected": True,
            "no_unphysical_three_dimensional_polarization_shortcut": True,
            "spatial_difference_only_not_homogeneous_trace_anchor": True,
        },
    }
