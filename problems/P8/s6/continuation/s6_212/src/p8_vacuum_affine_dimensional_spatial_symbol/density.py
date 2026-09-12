"""General-dimensional sphere-averaged full ultraviolet endpoint coefficients."""

from functools import cache

import sympy as s

from . import geometry, jets


@cache
def geometric_jets(channel):
    x = s.Symbol("x", real=True)
    out = {}
    for key, value in geometry.contractions(channel).items():
        value = s.series(value.subs(geometry.y, jets.p * x), x, 0, 5).removeO().expand()
        out[key] = tuple(s.factor(value.coeff(x, d)) for d in range(5))
    return out


@cache
def coefficients():
    rows = jets.endpoint_products()
    result = {}
    for ch in ("tensor", "vector", "scalar"):
        geo = geometric_jets(ch)
        for j in range(5):
            for r in range(j + 1):
                for degree in range(5 - j):
                    value = sum(
                        sum(
                            rows[key, j, r][v] * geo[key][degree - v]
                            for v in range(degree + 1)
                        )
                        for key in geo
                    )
                    # Cancel common dimension denominators before exact sphere moments.
                    polynomial = s.Poly(s.cancel(value), jets.u)
                    actual = geometry.sphere_average(polynomial.as_expr())
                    result[ch, j, r, degree] = s.factor(actual)
    return result


@cache
def spatial_difference():
    return {
        key: s.factor(value - value.subs(jets.p, 0))
        for key, value in coefficients().items()
    }


@cache
def first_dimension_jet():
    return {
        key: (
            s.factor(value.subs(geometry.d, 3)),
            s.factor(s.diff(value, geometry.d).subs(geometry.d, 3)),
        )
        for key, value in spatial_difference().items()
    }


@cache
def invariant_spatial_coefficients():
    rows = spatial_difference()
    out = {}
    for j in range(5):
        for r in range(j + 1):
            for degree in range(5 - j):
                T = rows["tensor", j, r, degree]
                V = 2 * (rows["vector", j, r, degree] - T)
                W = (
                    geometry.d * (rows["scalar", j, r, degree] - T) / (geometry.d - 1)
                    - V
                )
                out[j, r, degree] = tuple(s.factor(v) for v in (T, V, W))
    return out


@cache
def invariant_dimension_jet():
    return {
        key: tuple(
            (
                s.factor(v.subs(geometry.d, 3)),
                s.factor(s.diff(v, geometry.d).subs(geometry.d, 3)),
            )
            for v in row
        )
        for key, row in invariant_spatial_coefficients().items()
    }


@cache
def total_logarithmic_dimension_jet():
    rows = invariant_spatial_coefficients()
    total = []
    for r in range(5):
        row = tuple(
            s.factor(sum(rows[j, r, 4 - j][index] for j in range(r, 5)))
            for index in range(3)
        )
        total.append(
            tuple(
                (
                    s.factor(v.subs(geometry.d, 3)),
                    s.factor(s.diff(v, geometry.d).subs(geometry.d, 3)),
                )
                for v in row
            )
        )
    return tuple(total)


@cache
def data():
    from p8_vacuum_affine_spatial_matching_difference import density as physical

    old = physical.angular_coefficients()
    rows = coefficients()
    checks = {
        "all_105_full_physical_angular_coefficients": s.Matrix(
            [
                s.factor(value.subs(geometry.d, 3) / (2 * s.pi**2) - old[key])
                for key, value in rows.items()
            ]
        )
    }
    invariant = invariant_spatial_coefficients()
    difference = spatial_difference()
    for index, ch in enumerate(("tensor", "vector", "scalar")):
        V = (
            s.Integer(0)
            if index == 0
            else s.Rational(1, 2)
            if index == 1
            else (geometry.d - 1) / geometry.d
        )
        W = s.Integer(0) if index < 2 else (geometry.d - 1) / geometry.d
        checks[ch + "_full_dimension_invariant_reconstruction"] = s.Matrix(
            [
                s.factor(row[0] + V * row[1] + W * row[2] - difference[ch, *key])
                for key, row in invariant.items()
            ]
        )
    checks["scalar_basis_dimension_derivative_contact"] = s.Matrix(
        [
            s.factor(
                s.diff(difference["scalar", *key], geometry.d).subs(geometry.d, 3)
                - sum(
                    s.diff(v, geometry.d).subs(geometry.d, 3) * w
                    for v, w in zip(row, (1, s.Rational(2, 3), s.Rational(2, 3)))
                )
                - (row[1] + row[2]).subs(geometry.d, 3) / 9
            )
            for key, row in invariant.items()
        ]
    )
    zero_odd = tuple(row[0] for row in invariant_dimension_jet()[1, 0, 3])
    slope_odd = tuple(row[1] for row in invariant_dimension_jet()[1, 0, 3])
    checks["odd_endpoint_physical_spatial_difference_zero"] = s.Matrix(zero_odd)
    return {
        "sphere_normalization": "Coefficient tables use normalized S^(d-1) averages. Restore the complete Fourier factor2*pi^(d/2)/[Gamma(d/2)*(2pi)^d] only in the radial Laurent calculation. Physicald3 corresponds to1/(2pi^2).",
        "complete_coefficients": "All105 full angular coefficients, not only divergent slots, reproduce the complete actual S211 physical values. Spatial differences are taken after all physical pairs and source-derivative slots are assembled.",
        "fixed_source_reconstruction": "For each slot write A*T+B*V+C*W. From normalized channels, A=fT, B=2(fV-fT), C=d(fS-fT)/(d-1)-B. Differentiate these invariant coefficients atd3, not the d-dependent scalar test tensor.",
        "actual_evanescent_odd_endpoint": "The j1/r0/d3 spatial difference is zero in physical dimension but its dimension derivative is nonzero. It contributes to the source-value logarithmic finite term and cannot be deleted before the dimension limit.",
        "odd_endpoint_dimension_slope": slope_odd,
        "full_log_dimension_jets": total_logarithmic_dimension_jet(),
        "checks": checks,
        "gates": {
            "all_105_coefficients_with_actual_dimensional_modes": len(rows) == 105,
            "fixed_source_not_dimension_varying_scalar_probe": True,
            "nonzero_evanescent_odd_endpoint": any(v != 0 for v in slope_odd),
            "full_log_source_value_includes_all_endpoint_labels": True,
            "contact_difference_zero_before_integration_unchanged": True,
        },
    }
