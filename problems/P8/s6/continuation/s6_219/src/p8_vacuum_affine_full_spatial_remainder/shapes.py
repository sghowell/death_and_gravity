"""Full ordered scalar cutoff shapes from every retained physical endpoint."""

from functools import cache

import sympy as s
from p8_vacuum_affine_ordered_scalar_symbol import density, geometry, jets
from p8_vacuum_affine_subleading_band_conversion import angular, centered

t, u, p, m, a = jets.t, jets.u, jets.p, jets.m, jets.a
CHANNELS = density.CHANNELS


@cache
def raw_grades():
    """Azimuthal averages before radial integration, with the original leg."""
    products = jets.endpoint_products()
    out = {}
    for channel in CHANNELS:
        geo = density.geometric_jets(channel)
        for q in range(5):
            for r in range(q + 1):
                value = sum(
                    products[key, j, r][v] * geo[key][q - j - v]
                    for j in range(r, q + 1)
                    for key in geo
                    for v in range(q - j + 1)
                )
                out[channel, q, r] = s.factor(value.subs(geometry.d, 3))
    return out


def predicted_centered(channel):
    q = 3 * u**2 - 1
    ap = s.diff(a, t)
    if channel == "trace_trace":
        return (
            1 / (2 * a),
            -(m**2) * a / 4 - p**2 * (1 + u**2) / (16 * a) + 3 * (t**2 - 1) / 2,
            -ap / 8,
            -a / 8,
        )
    if channel == "trace_scalar":
        return (
            -q / (4 * a),
            -(m**2) * a * q / 4
            - p**2 * (u**2 - 2) / (8 * a)
            - 3 * (3 * t**2 + 1) * q / 4,
            -ap * q / 8,
            a * q / 16,
        )
    if channel == "scalar_trace":
        return (
            -q / (4 * a),
            -(m**2) * a * q / 4 - p**2 * (u**2 - 2) / (8 * a),
            ap * q / 4,
            a * q / 16,
        )
    raise ValueError("Require a fixed ordered scalar channel")


@cache
def shape_rows():
    """Coefficients of K^power, including all grazing value/derivative terms."""
    raw = raw_grades()
    out = {}
    for channel in CHANNELS:
        for power in range(4):
            for r in range(4 - power):
                value = sum(
                    -(p**h)
                    * angular.shape_coefficient(raw[channel, q, r], u, q, h)
                    / (2 * s.pi) ** 3
                    for q in range(r, 4)
                    for h in range(1, 5 - q)
                    if 4 - q - h == power
                )
                out[channel, power, r] = s.factor(value)
    return out


def predicted_shapes(channel):
    ap = s.diff(a, t)
    if channel == "trace_trace":
        cubic = -p / (16 * s.pi**2 * a)
        linear = (
            p * (m**2 * a / 32 + p**2 / (64 * a) + 3 * (1 - t**2) / 16) / s.pi**2,
            p * ap / (64 * s.pi**2),
            p * a / (64 * s.pi**2),
        )
    elif channel == "trace_scalar":
        cubic = p / (64 * s.pi**2 * a)
        linear = (
            p
            * (m**2 * a / 64 - 5 * p**2 / (256 * a) + 3 * (3 * t**2 + 1) / 64)
            / s.pi**2,
            p * ap / (128 * s.pi**2),
            -p * a / (256 * s.pi**2),
        )
    elif channel == "scalar_trace":
        cubic = p / (64 * s.pi**2 * a)
        linear = (
            p * (m**2 * a / 64 - 5 * p**2 / (256 * a)) / s.pi**2,
            -p * ap / (64 * s.pi**2),
            -p * a / (256 * s.pi**2),
        )
    else:
        raise ValueError("Require a fixed ordered scalar channel")
    return cubic, linear


def unit_factor(channel):
    if channel not in CHANNELS:
        raise ValueError("Require a fixed ordered scalar channel")
    return s.Rational(1, 3) if channel == "trace_trace" else 1 / s.sqrt(18)


@cache
def data():
    raw, shapes = raw_grades(), shape_rows()
    checks = {}
    for channel in CHANNELS:
        f0, *f2 = predicted_centered(channel)
        actual_f2 = [
            centered.centered_second(raw[channel, 0, 0], raw[channel, 2, 0], u, p),
            raw[channel, 2, 1],
            raw[channel, 2, 2],
        ]
        checks[channel + "_full_centered_leading_and_second"] = s.Matrix(
            [s.factor(raw[channel, 0, 0] - f0)]
            + [s.factor(x - y) for x, y in zip(actual_f2, f2)]
        )
        checks[channel + "_uncentered_first_grade_from_translation"] = s.factor(
            raw[channel, 1, 0] + p * (u * f0 + (1 - u**2) * s.diff(f0, u)) / 2
        )
        cubic, linear = predicted_shapes(channel)
        checks[channel + "_complete_original_mask_shapes"] = s.Matrix(
            [
                s.factor(
                    value
                    - (
                        cubic
                        if power == 3 and r == 0
                        else linear[r]
                        if power == 1 and r <= 2
                        else 0
                    )
                )
                for (ch, power, r), value in shapes.items()
                if ch == channel
            ]
        )
        absolute = lambda value: 2 * s.integrate(u * value, (u, 0, 1))
        checks[channel + "_independent_centered_shape_integrals"] = s.Matrix(
            [s.factor(-p * absolute(f0) / (8 * s.pi**2) - cubic)]
            + [
                s.factor(
                    absolute(
                        (p**3 * (3 - 5 * u**2) * f0 / 16 if r == 0 else 0)
                        - p * f2[r] / 2
                    )
                    / (4 * s.pi**2)
                    - linear[r]
                )
                for r in range(3)
            ]
        )
    tt = predicted_shapes("trace_trace")[1]
    left = predicted_shapes("trace_scalar")[1]
    right = predicted_shapes("scalar_trace")[1]
    checks["complete_ordered_linear_shape_formal_transpose"] = s.Matrix(
        [
            s.factor(tt[1] - s.diff(tt[2], t)),
            s.factor(right[2] - left[2]),
            s.factor(right[1] - 2 * s.diff(left[2], t) + left[1]),
            s.factor(right[0] - left[0] + s.diff(left[1], t) - s.diff(left[2], t, 2)),
        ]
    )
    return {
        "full_scalar_shapes": {
            channel: {
                "cubic": predicted_shapes(channel)[0],
                "linear_source_jets": predicted_shapes(channel)[1],
                "unit_channel_factor": unit_factor(channel),
            }
            for channel in CHANNELS
        },
        "derivation": "Use every full corrected j/degree/source-jet product from S216 before angular integration. Sum j+degree=q0,...,3 and apply the original uncentered shell hemisphere plus exact positive grazing strip. All quadratic and finite terms vanish in the complete sum, not in individual rows. Independently integrate the centered q0/q2 formulas against abs(u).",
        "ordered_scalar_terms": "The LC,CL,CC odd nonlog products are nonzero. In particular the two cross source-jet1 terms differ and their source-jet0 difference is the required formal-transpose derivative. S217's tracefree odd-nonlog zero is not used for these channels.",
        "scope": "These are the unnormalized I/I,I/S,S/I channels, S=diag(-1,-1,2), with unit factors1/3,1/sqrt18,1/sqrt18. The q4 shell stays in the decaying error and its bulk logarithm stays in S216. No contact value or finite prescription is removed.",
        "checks": checks,
        "gates": {
            "all_45_scalar_radial_source_jet_grades": len(raw) == 45,
            "all_thirty_scalar_cutoff_power_source_jet_shapes": len(shapes) == 30,
            "ordered_cross_first_jets_differ": left[1] != right[1],
            "trace_trace_odd_nonlog_endpoint_not_zero": jets.endpoint_products()[
                "CC", 1, 0
            ][1]
            != 0,
            "full_quadratic_and_finite_cancellation_checked": all(
                value == 0
                for (ch, power, r), value in shapes.items()
                if power in (0, 2)
            ),
            "all_shape_angular_degrees_at_most_eight": all(
                s.Poly(value, u).degree() <= q + 4 for (ch, q, r), value in raw.items()
            ),
        },
    }
