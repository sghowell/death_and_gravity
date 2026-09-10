"""Pointwise real-momentum routing; no complex contour translation."""

from functools import cache

import sympy as sp

from . import halfplane


def channel_margin(center, radius, sigma=1):
    c, R, sig = map(halfplane.exact, (center, radius, sigma))
    if c < 0 or R < 0 or sig < 1:
        raise ValueError("Require nonnegative real center/radius and sigma>=1")
    line_real_upper = (c + R) / 4
    return {
        "external_real_center": c,
        "external_disc_radius": R,
        "line_invariant_real_part_upper": line_real_upper,
        "halfplane_sigma": sig,
        "halfplane_margin": sig - line_real_upper,
        "selected_symmetric_bubble_routing_covered": bool(line_real_upper <= sig),
        "scope": "Only q plus or minus half of this channel energy, q real; not arbitrary graph momentum routings.",
    }


@cache
def data():
    a, b, q0, q1, q2, q3 = sp.symbols("a b q0 q1 q2 q3", real=True)
    E = a + sp.I * b
    external = 4 * E**2
    spatial = q1 * q1 + q2 * q2 + q3 * q3
    lines = {}
    checks = {}
    for sign, label in ((1, "plus"), (-1, "minus")):
        line = -((q0 + sign * sp.I * E) ** 2 + spatial)
        target = -((q0 - sign * b) ** 2) - spatial + a * a
        lines[label] = line
        checks[label + "_real_part_complete_square"] = sp.expand(sp.re(line) - target)
        checks[label + "_nonnegative_momentum_margin"] = sp.expand(
            a * a - sp.re(line) - (q0 - sign * b) ** 2 - spatial
        )
    modulus = 4 * (a * a + b * b)
    checks.update(
        {
            "external_modulus_squared": sp.expand(
                sp.re(external) ** 2 + sp.im(external) ** 2 - modulus**2
            ),
            "energy_real_part_from_invariant": sp.expand(
                a * a - (modulus + sp.re(external)) / 8
            ),
            "unit_forward_disc_line_upper": sp.Rational(3, 4) - (3 + 3) / sp.Integer(8),
            "unit_forward_disc_halfplane_margin": 1
            - sp.Rational(3, 4)
            - sp.Rational(1, 4),
        }
    )
    v = sp.symbols("forward_displacement")
    checks["crossed_forward_disc_center"] = sp.expand((4 - (2 + v)) - 2 + v)
    return {
        "energy": E,
        "external_channel_invariant": external,
        "internal_line_invariants": lines,
        "real_momentum_variables": [q0, q1, q2, q3],
        "selected_external_disc": "|s-2|<=1, u=4-s in the same disc; E=sqrt(s)/2. The sign of the square root only exchanges the two routed lines.",
        "pointwise_line_halfplane": "Re(s_line)=-(q0 minus/plus Im E)^2-|q_spatial|^2+(Re E)^2 <= (|s|+Re s)/8 <=3/4<1.",
        "scope": "The symmetric one-loop bubble backbone is covered at all real internal q, including |q| arbitrarily larger than mF. This is not a claim about every multi-loop routing or a complex shift of q.",
        "checks": checks,
    }
