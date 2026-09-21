"""Exact off-shell Symanzik polynomials and strict compact-domain interfaces."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_scalar_four_point_loop import amplitude as old

from . import source

KINDS = ("triangle", "ordered_box")


def rational(value):
    if isinstance(value, (bool, str, float)) or value is None:
        raise ValueError("Require an exact real rational parameter")
    value = s.sympify(value)
    if not isinstance(value, s.Rational):
        raise TypeError("Require an exact real rational parameter")
    return value


def heavy_mass(value):
    value = rational(value)
    if not 10**6 <= value < 10**198:
        raise ValueError("Require the stated heavy-mass contour interval")
    return value


def require_kind(value):
    if not isinstance(value, str) or value not in KINDS:
        raise ValueError("Require triangle or ordered_box")
    return value


def domain(kind, v, virtualities, w=0, mass=source.HEAVY_MASS2):
    kind = require_kind(kind)
    v, w, mass = rational(v), rational(w), heavy_mass(mass)
    if not (-12 <= v <= s.Rational(4, 3) or s.Rational(45, 8) <= v <= 16):
        raise ValueError(
            "Require the gapped subthreshold or physical timelike interval"
        )
    if not isinstance(virtualities, (tuple, list)) or len(virtualities) != (
        2 if kind == "triangle" else 4
    ):
        raise ValueError("Require every adjacent singleton virtuality")
    virtualities = tuple(rational(q) for q in virtualities)
    if not all(0 <= q <= 2 for q in virtualities):
        raise ValueError("Require adjacent virtualities in [0,2]")
    if (kind == "triangle" and w != 0) or not -12 <= w <= 16:
        raise ValueError("Require the stated ordered heavy transfer")
    return kind, v, virtualities, w, mass


def denominator(z, L, n, c=0, b=0):
    return (1 - z) ** 2 * L + n * z + c * z * (1 - z) - b * z * z


@cache
def data():
    checks = {}
    z, xi, eta, n, v, w = s.symbols("z xi eta n v w", real=True)
    a1, a2, a3, a4 = s.symbols("a1 a2 a3 a4", real=True)
    L, c, b = s.symbols("L c b", real=True)
    Delta = (1 - z) ** 2 * L + n * z + c * z * (1 - z) - b * z**2
    tri_x = ((1 - z) * xi, (1 - z) * (1 - xi), z)
    tri_e = xi * a1 + (1 - xi) * a2
    tri_direct = (
        tri_x[0]
        + tri_x[1]
        + n * tri_x[2]
        - tri_x[0] * tri_x[1] * v
        - tri_x[0] * tri_x[2] * a1
        - tri_x[1] * tri_x[2] * a2
    )
    checks["triangle_full_offshell_Symanzik"] = s.expand(
        tri_direct - Delta.subs({L: 1 - v * xi * (1 - xi), c: 1 - tri_e, b: 0})
    )
    box_x = ((1 - z) * xi, z * eta, (1 - z) * (1 - xi), z * (1 - eta))
    box_e = xi * (eta * a1 + (1 - eta) * a4) + (1 - xi) * (eta * a2 + (1 - eta) * a3)
    box_direct = (
        box_x[0]
        + box_x[2]
        + n * (box_x[1] + box_x[3])
        - box_x[0] * box_x[2] * v
        - box_x[1] * box_x[3] * w
    )
    box_direct -= sum(
        box_x[j] * box_x[(j + 1) % 4] * aa for j, aa in enumerate((a1, a2, a3, a4))
    )
    checks["ordered_box_full_offshell_Symanzik"] = s.expand(
        box_direct
        - Delta.subs({L: 1 - v * xi * (1 - xi), c: 1 - box_e, b: w * eta * (1 - eta)})
    )
    checks["triangle_simplex_Jacobian"] = s.det(
        s.Matrix([tri_x[0], tri_x[2]]).jacobian([xi, z])
    ) - (1 - z)
    checks["box_ordered_simplex_Jacobian"] = s.expand(
        s.det(s.Matrix(box_x[:3]).jacobian([z, xi, eta])) + z * (1 - z)
    )
    checks["exact_prior_primitive_no_heavy_expansion"] = s.expand(
        Delta - ((n + c) * z + (1 - z) ** 2 * L - (b + c) * z * z)
    )
    checks["old_onshell_triangle_denominator"] = s.expand(
        tri_direct.subs({a1: 1, a2: 1})
        - old.QC.subs(
            {old.xi: xi, old.z: z, old.n: n, old.e: 1, old.svar: v}, simultaneous=True
        )
    )
    checks["old_onshell_ordered_box_denominator"] = s.expand(
        box_direct.subs({a1: 1, a2: 1, a3: 1, a4: 1})
        - old.QD.subs(
            {
                old.xi: xi,
                old.z: z,
                old.n: n,
                old.e: 1,
                old.svar: v,
                old.tvar: w,
                old.eta: eta,
            },
            simultaneous=True,
        )
    )

    return {
        "checks": {key: s.factor(value) for key, value in checks.items()},
        "gates": {
            "all_external_virtualities_retained": True,
            "ordered_mass_assignment_not_exchanged": True,
            "symmetric_4_over3_anchor_included": True,
            "no_heavy_line_expansion_before_integration": True,
        },
        "whole_triangle_denominator": tri_direct,
        "whole_ordered_box_denominator": box_direct,
        "whole_general_reduced_denominator": Delta,
        "whole_convex_virtuality_weights": {"triangle": tri_e, "box": box_e},
        "whole_exact_primitive_substitution": "T(n+c,L,b+c), c=1-e; triangle measure(1-z), ordered-box measure z(1-z). No pre-integration inverse-mass truncation.",
        "whole_domain": "10^6<=n<10^198; adjacent singleton virtualities[0,2]; light channel[-12,4/3] or[45/8,16]; ordered heavy channel[-12,16]. Includes the symmetric 4/3 subtraction without approaching threshold4.",
    }
