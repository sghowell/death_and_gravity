"""Exact compact forward-disc points and actual heavy-mass applicability."""

from fractions import Fraction
from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model

from . import forests, graphs


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational parameter")
    return sp.Rational(value)


def point(kind, choices, weights, real=0, imaginary=0, mass=None):
    graph = graphs.refinement(kind, choices)
    if type(weights) is not tuple or len(weights) != len(graph["edges"]):
        raise TypeError("Require a native parameter tuple of the exact edge count")
    alpha = tuple(rational(value) for value in weights)
    if any(value <= 0 for value in alpha) or sum(alpha) != 1:
        raise ValueError("Require the strictly positive normalized parameter simplex")
    re, im = rational(real), rational(imaginary)
    if re * re + im * im > 1:
        raise ValueError("Require the closed unit complex forward disc")
    M = data()["actual_heavy_mass_squared"] if mass is None else rational(mass)
    if M < 16:
        raise ValueError("Require heavy mass squared at least sixteen")
    d = forests.data(kind, choices)
    substitution = dict(zip(d["parameters"], alpha))
    U = d["U"].subs(substitution)
    H = d["H"] if type(d["H"]) is int else d["H"].subs(substitution)
    F = d["F"].subs(substitution).subs({d["M"]: M, d["v"]: re + sp.I * im})
    real_F, imaginary_F = sp.expand(F).as_real_imag()
    lower = U * (sp.Rational(1, 4) + (M - 16) * H)
    derivative = abs(d["crossing_derivative_polynomial"].subs(substitution))
    return {
        "kind": kind,
        "choices": choices,
        "normalized_parameters": alpha,
        "forward_real_increment": re,
        "forward_imaginary_increment": im,
        "heavy_mass_squared": M,
        "actual_parent_mass": M == data()["actual_heavy_mass_squared"],
        "U": U,
        "real_F": real_F,
        "imaginary_F": imaginary_F,
        "proved_real_F_lower": lower,
        "absolute_crossing_derivative_upper": derivative,
        "relative_crossing_derivative_upper": derivative / lower,
        "scope": "An interior denominator enclosure, not an integrated or renormalized two-loop amplitude bound.",
    }


@cache
def data():
    p = model.data()["actual_parameters"]
    M = p["heavy_mass_squared"]
    U, H, total, v = sp.symbols(
        "positive_U heavy_parameter_sum total_parameter_sum real_disc_increment",
        real=True,
    )
    Ps, Pt, Pu = sp.symbols(
        "positive_s_forest positive_t_forest positive_u_forest", nonnegative=True
    )
    return {
        "actual_heavy_mass_squared": M,
        "minimum_supported_heavy_mass_squared": sp.Integer(16),
        "uniform_interior_simplex_gap": sp.Rational(1, 4),
        "uniform_relative_crossing_derivative_bound": sp.Integer(4),
        "homogeneous_real_F_lower": U * (total / 4 + (M - 16) * H),
        "simplex_real_F_lower": U * (sp.Rational(1, 4) + (M - 16) * H),
        "disc_coarse_kinematic_gap": (1 - v) * Ps + 3 * Pt + (1 + v) * Pu,
        "scope": "Full heavy exchanges in the canonical parent. No ultraviolet boundary estimate, finite counterterm or complete two-loop forward coefficient is supplied by this denominator bound.",
        "checks": {
            "simplex_homogeneous_lower_normalization": sp.expand(
                (U * (total / 4 + (M - 16) * H)).subs(total, 1)
                - U * (sp.Rational(1, 4) + (M - 16) * H)
            ),
            "uniform_relative_derivative_bound_from_gap": sp.Rational(1, 1)
            / sp.Rational(1, 4)
            - 4,
        },
        "bounds": {
            "actual_heavy_mass_in_supported_domain": M > 16,
            "actual_heavy_mass_lift_positive": M - 16 > 0,
            "simplex_gap_strictly_positive": sp.Rational(1, 4) > 0,
        },
    }


def bad_cases():
    base = ("double_bubble", (0, 0, 0), (sp.Rational(1, 4),) * 4)
    rows = []
    for module in (graphs, forests):
        call = module.refinement if module is graphs else module.data
        for i, kind in enumerate((True, None, "sunset", 1)):
            rows.append(
                (
                    module.__name__.split(".")[-1] + "_kind_" + str(i),
                    call,
                    (kind, (0, 0, 0)),
                )
            )
        for i, choices in enumerate(
            (
                [0, 0, 0],
                (0, 0),
                (0, 0, 0, 0),
                (True, 0, 0),
                (sp.Integer(1), 0, 0),
                (4, 0, 0),
                (-1, 0, 0),
                (1.0, 0, 0),
            )
        ):
            rows.append(
                (
                    module.__name__.split(".")[-1] + "_choices_" + str(i),
                    call,
                    ("double_bubble", choices),
                )
            )
    for i, bad in enumerate((True, 1.0, sp.Float(1), "1", None, sp.oo, sp.I, sp.nan)):
        rows.append(
            (
                "inexact_parameter_" + str(i),
                point,
                (base[0], base[1], (bad,) + base[2][1:]),
            )
        )
        rows.append(("inexact_disc_" + str(i), point, base + (bad, 0)))
    rows += [
        ("wrong_simplex_container", point, (base[0], base[1], list(base[2]))),
        ("wrong_simplex_length", point, (base[0], base[1], base[2][:3])),
        ("unnormalized_simplex", point, (base[0], base[1], (1,) * 4)),
        ("boundary_simplex", point, (base[0], base[1], (0, 0, 0, 1))),
        ("negative_simplex", point, (base[0], base[1], (-1, 1, 0, 1))),
        ("outside_real_disc", point, base + (2, 0)),
        ("outside_complex_disc", point, base + (1, 1)),
        ("unsupported_heavy_mass", point, base + (0, 0, 15)),
        ("inexact_heavy_mass", point, base + (0, 0, 16.0)),
        ("boolean_heavy_mass", point, base + (0, 0, True)),
    ]
    return rows
