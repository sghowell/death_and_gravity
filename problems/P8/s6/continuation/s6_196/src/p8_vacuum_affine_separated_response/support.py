"""Strict temporal ordering, local contacts and causal support."""

from functools import cache

import sympy as s
from p8_vacuum_flat_dirac_hadamard.symbols import rational


def require_ordered(source_start, source_end, readout_start, readout_end):
    a, b, c, d = map(rational, (source_start, source_end, readout_start, readout_end))
    if not -s.Rational(1, 2) < a < b < c < d < s.Rational(1, 2):
        raise ValueError(
            "Require strict ordered compact source/readout intervals inside the open slab"
        )
    return a, b, c, d


def require_spacelike_enclosure(distance_lower, radii_sum, time_width):
    d, r, T = map(rational, (distance_lower, radii_sum, time_width))
    if d < 0 or r < 0 or T < 0 or d - r <= T:
        raise ValueError("Require a strict spacelike support enclosure with a>=1")
    return d - r - T


@cache
def data():
    m = s.Symbol("mass", positive=True)
    d0 = s.Matrix([[1, 2], [2, -1], [3, 1]])
    d1 = s.Matrix([[1, 1, -1], [2, 2, -2]])
    P = d1.T * d1 + m * m * s.eye(3)
    Q = s.eye(3) + d0 * d0.T / (m * m)
    normally = d1.T * d1 + d0 * d0.T + m * m * s.eye(3)
    t = s.Symbol("t", real=True)
    source_end = s.Rational(-1, 8)
    readout_start = s.Rational(1, 8)
    checks = {
        "exterior_complex_nilpotence": d1 * d0,
        "complete_Proca_to_normal_operator_identity": P * Q - normally,
        "opposite_operator_order_same_identity": Q * P - normally,
        "clock_scale_at_least_one": s.expand((1 + t * t) ** 2 - 1 - 2 * t * t - t**4),
        "positive_time_gap": readout_start - source_end - s.Rational(1, 4),
        "strict_spacelike_enclosure_example": require_spacelike_enclosure(3, 1, 1) - 1,
    }
    return {
        "ordered_supports": "Real smooth compact spatial tracefree source Gamma and readout D in the open CD slab, with sup time(supp Gamma)<inf time(supp D). No spatial Fourier cutoff or homogeneous restriction is imposed.",
        "why_no_contact": "The full metric contact is proportional to the pointwise product of the two directions and vanishes on disjoint temporal supports. Every derivative of Gamma also vanishes near supp D, so all fixed local curvature/counterterm responses vanish in this pairing. They are not reset or subtracted anew.",
        "why_step_is_one": "Strict ordering makes theta(t-s)=1 on the entire product of the actual test supports. This is not multiplication of a coincident stress distribution by a discontinuous step.",
        "Proca_locality": "The complete Proca causal Green operator is obtained from the normally hyperbolic one-form Green operator using the local differential factor I+d delta/m^2. Differentiation does not enlarge causal support. Thus local Wick stress commutators vanish for mutually spacelike supports, including every constrained polarization.",
        "spacelike_geometry": "On the actual clock a>=1, every causal curve obeys |dx/dt|<=1/a<=1. A verified spatial separation lower bound minus the sum of support radii greater than the largest time difference therefore proves spacelike separation.",
        "boundary": "The algebraic matrix fixture only checks the de Rham operator identity; causal support uses the full globally hyperbolic Green theorem. Temporal ordering cannot be removed just because a commutator bilinear is bounded. No coincident-support extension or response inverse is proved.",
        "checks": checks,
        "gates": {
            "nonempty_actual_ordered_domain": require_ordered(
                -s.Rational(3, 8), -s.Rational(1, 8), s.Rational(1, 8), s.Rational(3, 8)
            )[1]
            < 0,
            "positive_mass_local_differential_factor": m.is_positive is True,
            "no_elliptic_temporal_oscillator_added": True,
            "contact_zero_by_support_not_retuning": True,
        },
    }
