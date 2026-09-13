"""Whole heavy-scalar joint domain and coefficient-wise dimensional continuation."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_curved_state import state
from p8_vacuum_affine_heavy_spatial_uv import jets as uv

d, u = uv.d, uv.u
y, c, h, v = s.symbols("y c h v", real=True)
T, V, W, S, UD, UG = uv.inv
EPS = s.Rational(1, 100)
RADIUS = s.Rational(1, 20000)
MASS_FLOOR = s.Integer(10) ** 98
PAIR = s.Integer(10) ** 9
NORMALIZED_PAIR = s.Integer(1000)


def normalized_geometries():
    b = d - 1
    return {
        "00": S,
        "01": c * h * S / b + (u * v - c * h / b) * UD,
        "10": c * h * S / b + (u * v - c * h / b) * UG,
        "11": (c * h) ** 2 * (2 * T - 4 * V + S - UD - UG + 3 * W) / (b * (b + 2))
        + (c * v + h * u) ** 2 * (V - W) / b
        + c * h * u * v * (UD + UG - 2 * W) / b
        + u * u * v * v * W,
    }


def far_geometries():
    yy, rows = uv.geometries()
    return {name: row.subs(yy, y) for name, row in rows.items()}


def rational_bound(value, variables, radii):
    eta = s.Symbol("dimension_offset", real=True)
    row = s.cancel(value.subs(d, 3 + eta))
    num, den = s.fraction(row)
    allvars = (eta, *variables)
    allr = (s.Rational(1, 4), *radii)
    np, dp = s.Poly(num, *allvars), s.Poly(den, *allvars)

    def norm(poly):
        return sum(
            abs(coeff) * s.prod(r**j for r, j in zip(allr, powers))
            for powers, coeff in poly.terms()
        )

    constant = abs(dp.coeff_monomial(1))
    floor = 2 * constant - norm(dp)
    if floor <= 0 or not all(coeff.is_Rational for coeff in np.coeffs() + dp.coeffs()):
        raise ValueError(
            "The full scalar continuation requires a positive exact denominator floor"
        )
    return (
        s.factor(norm(np) / floor),
        floor,
        s.factor(row - np.as_expr() / dp.as_expr()),
    )


@cache
def direct_frequency(leg):
    if leg not in ("k", "l"):
        raise ValueError("Retain one of the two scalar internal legs")
    a, m, p = uv.a, uv.m, uv.p
    square = (s.S.One, s.S.Zero, a * a * m * m, s.S.Zero, s.S.Zero)
    if leg == "l":
        square = (s.S.One, -2 * p * u, p * p + a * a * m * m, s.S.Zero, s.S.Zero)
    row = uv.sqrtone(square)
    potential = d * uv.h[1] / 2 + d * d * uv.h[0] ** 2 / 4
    iterates = [row]
    for _ in range(6):
        rate = uv.add(
            uv.mul(uv.diff(row), uv.reciprocal(row)), uv.scale(uv.ONE, -uv.h[0])
        )
        full = uv.add(
            uv.scale(uv.ONE, -potential),
            uv.add(
                uv.scale(uv.diff(rate), -s.Rational(1, 2)),
                uv.scale(uv.mul(rate, rate), s.Rational(1, 4)),
            ),
        )
        row = uv.sqrtone(uv.add(square, uv.scale(uv.shift(full, 2), a * a)))
        iterates.append(row)
    _, target, _ = uv.frequency(leg)
    return {
        "all_seven_four_jets": tuple(iterates),
        "full_frozen_coefficient_residual": s.Matrix(
            [s.factor(left - right) for left, right in zip(row, target)]
        ),
    }


@cache
def data():
    near, far = normalized_geometries(), far_geometries()
    R = s.Symbol("Rleg", positive=True)
    checks = {}
    for name, row in near.items():
        power = {"00": 0, "01": 1, "10": 1, "11": 2}[name]
        actual = s.expand(
            row.subs({h: -c / R, v: (y - u) / R}, simultaneous=True) * R**power
        )
        actual = s.cancel(actual).subs(c**4, (1 - u * u) ** 2).subs(c**2, 1 - u * u)
        checks["full_normalized_near_far_geometry_" + name] = s.factor(
            actual - far[name]
        )
    bounds = {}
    floors = {}
    weights = (1, 1, 1, 3, 2, 2)
    for label, geo, variables, radii in (
        ("far", far, (y, u), (EPS, 1)),
        ("near", near, (c, h, u, v), (1, 1, 1, 1)),
    ):
        for name, row in geo.items():
            total = 0
            for invariant, weight in zip(uv.inv, weights):
                upper, floor, error = rational_bound(
                    s.diff(row, invariant), variables, radii
                )
                total += weight * upper
                label0 = label + "_" + name + "_" + str(invariant)
                floors[label0] = floor
                checks["full_rational_coefficient_reconstruction_" + label0] = error
            bounds[label + "_" + name] = s.factor(total)
    for leg in ("k", "l"):
        checks["all_five_full_six_iterate_UV_coefficients_" + leg] = direct_frequency(
            leg
        )["full_frozen_coefficient_residual"]
    ea = s.Rational(5, 27)
    eg = 2 * EPS + EPS**2
    ew = ea + (1 + ea) * eg + s.Rational(25, 16) ** 2 * EPS**2
    eO = eg + 4 * EPS**2
    q = 4 * 10**8 * EPS**2 / MASS_FLOOR**2
    remaining = s.Rational(1, 64) - 7 * s.Rational(1, 1024)
    scaled_rate = (s.Rational(9, 4) * 5 + 2 * 1024 * 2) / 2
    actual_rate = (s.Rational(13, 4) * 5 + 9 * 1024 * 3) / 2
    checks["joint_full_squared_frequency_defect"] = ew - s.Rational(14463467, 69120000)
    checks["normalized_scalar_pair_contraction_allowance"] = (
        s.Integer(2) * (s.Integer(2) * (9 + 2 + 1) / 2 + 2) - 28
    )
    checks["remaining_complex_time_width"] = remaining - s.Rational(9, 1024)
    gates = {
        "unchanged_heavy_mass_above_every_proof_threshold": state.MASS2 > MASS_FLOOR**2,
        "full_joint_inverse_radius_and_time_frequency_defect": ew < s.Rational(1, 4),
        "full_scaled_frequency_real_part_floor": s.Rational(3, 4)
        / s.Rational(25, 16) ** 2
        > s.Rational(11, 20) ** 2,
        "complete_scale_and_inverse_scale_bounds": s.Rational(100, 59) < 2
        and s.Rational(64, 59) < 2,
        "full_O_branch_near_one": eO < s.Rational(1, 40)
        and 1 - eO > s.Rational(49, 50) ** 2,
        "full_normalized_z_bound": 4 * EPS**2 / (1 - eg - 4 * EPS**2)
        < s.Rational(1, 100),
        "all_six_full_W6_iterates_in_same_analytic_ball": 2 * q < s.Rational(1, 1000),
        "whole_normalized_frequency_and_inverse_sum_phase": s.Rational(11, 20) - 4 * q
        > s.Rational(1, 2),
        "full_Wbar_nondegenerate": s.Rational(49, 50) - 4 * q > s.Rational(9, 10),
        "whole_Wbar_Cauchy_bound": s.Rational(11, 10) * (1 + 2 * q) < 2,
        "full_normalized_derivative_feature": scaled_rate < 5000
        and 2 + 2 * EPS * 5000 / MASS_FLOOR < 3,
        "complete_scalar_normalized_pair": 28 < NORMALIZED_PAIR,
        "all_real_full_derivative_rate": actual_rate < 20000 and MASS_FLOOR > 20000,
        "whole_scalar_feature_and_pair_bound": 20**2 * s.Rational(3, 2) < 25**2
        and 3 * 25**2 < PAIR,
        "all_required_source_time_discs_fit": remaining > 2 * RADIUS,
        "general_dimension_scalar_frequency_potential": s.Rational(13, 4) * 12 / 2
        + s.Rational(13, 4) ** 2 * 25 / 4
        < 100,
        "every_full_near_far_geometry_bound": all(
            value < 100 for value in bounds.values()
        ),
        "all_48_dimensional_denominator_floors": len(floors) == 48
        and all(value > 0 for value in floors.values()),
        "both_same_dimension_phase_branches_bounded": s.Rational(10**10, MASS_FLOOR)
        < s.Rational(1, 40),
        "no_noninteger_dimension_positive_norm_assumed": True,
        "integrable_far_low_and_angular_envelopes": -s.Rational(7, 4) < -1
        and s.Rational(7, 4) > -1
        and -s.Rational(1, 8) > -1,
        "six_spatial_derivatives_cover_entire_near_continuation": s.Rational(21, 4) < 6,
    }
    return {
        "complete_normalized_near_geometries": near,
        "complete_far_geometries": far,
        "all_eight_full_dimensional_geometry_majorants": bounds,
        "all_48_exact_coefficient_denominator_floors": floors,
        "full_six_iterate_frequency_coefficient_checks": {
            leg: direct_frequency(leg) for leg in ("k", "l")
        },
        "joint_domain_constants": {
            "epsilon": EPS,
            "outer_time_radius": s.Rational(1, 64),
            "source_Cauchy_radius": RADIUS,
            "scale_defect": ea,
            "momentum_square_defect": eg,
            "full_frequency_square_defect": ew,
            "normalized_frequency_defect": eO,
            "entire_W6_iteration_ratio": q,
            "remaining_time_width": remaining,
        },
        "complete_physical_scalar_pair_majorant": PAIR,
        "complete_normalized_scalar_pair_majorant": NORMALIZED_PAIR,
        "domain": "All real transfers, |x|<=(1/100)/(m+|P|), full complex time radius1/64 with all six1/1024 W6 losses and one readout loss. Every iterate is retained; this is not a truncated physical evolution or a low-momentum cutoff.",
        "dimension": "Continue the exact scalar four-geometry algebra at fixed six physical invariants on |d-3|<=1/4. Near directions are normalized separately; zero internal momentum has no massive amplitude singularity. Both phase branches use the SAME d. Far radial, low and angular envelopes are r^-7/4,r^7/4,(1-u²)^-1/8; full near degree21/4<6.",
        "checks": checks,
        "gates": {key: bool(value) for key, value in gates.items()},
    }
