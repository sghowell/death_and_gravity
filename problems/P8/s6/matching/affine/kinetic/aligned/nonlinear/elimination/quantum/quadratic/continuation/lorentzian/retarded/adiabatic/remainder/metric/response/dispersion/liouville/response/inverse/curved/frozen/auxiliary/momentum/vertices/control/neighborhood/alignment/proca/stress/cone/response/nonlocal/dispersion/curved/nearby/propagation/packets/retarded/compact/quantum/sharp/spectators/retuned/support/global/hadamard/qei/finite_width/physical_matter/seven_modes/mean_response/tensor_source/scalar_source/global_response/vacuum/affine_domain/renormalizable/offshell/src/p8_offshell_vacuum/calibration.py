"""Exact spectral-window and canonical-jet-cube calibration."""

from fractions import Fraction

import sympy as sp
from p8_exceptional_vacuum import heavy
from p8_polynomial_vacuum import calibration as parent

from . import jets, norms


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational")
    return sp.Rational(value)


def point(radius=heavy.CHANNEL_RADIUS, jet_bound=1):
    r, B = rational(radius), rational(jet_bound)
    d = parent.point()
    D = d["D"]
    if not 0 < r < D:
        raise ValueError(
            "Require a nonempty spectral window strictly within the heavy gap"
        )
    if B < 0:
        raise ValueError("Require nonnegative canonical jet bound")
    C = norms.data()["actual_C_R"]
    CL = norms.data()["actual_C_truncated_quartic"]
    change = 3**6 * C * B**3
    field_error = sp.Rational(37, 2) * C * C * B**6 + CL * B**4 * (
        (1 + 3**6 * C * B * B) ** 4 - 1
    )
    return {
        "radius": r,
        "jet_bound": B,
        "actual_D": D,
        "spectral_quartic_action_error_coefficient": d["cubic_squared"]
        * r**4
        / (8 * D**5 * (1 - r / D)),
        "six_derivative_jet_change_bound": change,
        "pointwise_field_redefinition_error": field_error,
        "norms_are_separate": "Do not add the action L2-squared coefficient and pointwise jet-cube error without a common function-space/volume bound",
    }


def bad_cases():
    invalid = (True, False, 1.0, sp.Float(1), "1", sp.I, sp.oo, -sp.oo, sp.zoo, sp.nan)
    rows = []
    for i, v in enumerate(invalid):
        rows += [
            ("radius_" + str(i), point, (v,)),
            ("jet_bound_" + str(i), point, (1, v)),
        ]
    D = parent.point()["D"]
    rows += [
        ("spectral_domain_" + str(i), point, (r,))
        for i, r in enumerate((0, -1, D, D + 1))
    ]
    rows += [("negative_jet_bound", point, (1, -1))]
    rows += [
        ("axis_" + str(i), jets.derivative, (jets.PHI, axis))
        for i, axis in enumerate(
            (True, False, 0.0, sp.Integer(0), "0", -1, 4, None, sp.I)
        )
    ]
    rows += [("jet_order_exceeded", jets.derivative, (jets.JETS[(6, 0, 0, 0)], 0))]
    return rows
