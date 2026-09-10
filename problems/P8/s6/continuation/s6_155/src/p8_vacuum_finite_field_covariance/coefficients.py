"""Minimal-pole projection after a full regulated finite Phi field change."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_proper_references import conversion as proper
from p8_vacuum_full_one_loop_matching import scheme as first_scheme


def field_series(weight):
    h, e = s.symbols("h epsilon")
    X, p1, p22, p21, k0, k1, k2, t0, t1 = s.symbols("X p1 p22 p21 k0 k1 k2 t0 t1")
    w = s.Rational(weight)
    k = k0 + e * k1 + e**2 * k2
    t = t0 + e * t1
    total = X + h * p1 / e + h**2 * (p22 / e**2 + p21 / e)
    factor = 1 - h * w * k + h**2 * (w * (w + 1) * k**2 / 2 - w * t)
    transformed = s.series(total * factor, h, 0, 3).removeO().expand()
    second = transformed.coeff(h, 2)
    finite_second = s.expand(second).coeff(e, 0)
    expected = (w * (w + 1) * k0**2 / 2 - w * t0) * X - w * k1 * p1
    return {
        "weight": w,
        "transformed_regulated_total_coefficient": transformed,
        "finite_second_coefficient": finite_second,
        "expected_finite_second_coefficient": expected,
        "second_poles": p22 / e**2 + (p21 - w * k0 * p1) / e,
        "checks": {
            "first_finite_map": transformed.coeff(h, 1).coeff(e, 0) + w * k0 * X,
            "second_finite_map": s.expand(finite_second - expected),
            "second_double_pole": second.coeff(e, -2) - p22,
            "second_simple_pole": second.coeff(e, -1) - p21 + w * k0 * p1,
        },
    }


@cache
def data():
    L, G, y, Y, a, Cf, N, Q, M = s.symbols("L G y Y a Cf N Q M")
    k0, k1, t0 = s.symbols("k0 k1 new_second_slope")
    weights = {"L": 2, "G": 1, "y": s.Rational(1, 2), "M": 0}
    parameters = {"L": L, "G": G, "y": y, "M": M}
    residues = {
        "L": (3 * L**2 / 2 - 24 * N * Y**2) / Q,
        "G": L * G / (2 * Q),
        "y": (Y - 4 * a * Cf) * y / Q,
        "M": G**2 / (2 * Q),
    }
    second = {
        key: (s.Rational(w) * (w + 1) * k0**2 / 2 - w * t0) * parameters[key]
        - w * k1 * residues[key]
        for key, w in weights.items()
    }
    first = {key: -w * k0 * parameters[key] for key, w in weights.items()}
    g2 = 2 * G * second["G"] + first["G"] ** 2
    Y2 = 2 * y * second["y"] + first["y"] ** 2
    g = s.Symbol("g")
    tree = 2 * g / (M - 2) ** 3
    comm_g = -k1 * L * g / Q
    checks = {}
    old = first_scheme.data()
    Ibar = old["reference_symbols"]["Ibar"]
    replacements = {
        old["parameters"][key]: value
        for key, value in (("L", L), ("G", G), ("M", M), ("Y", Y))
    }
    for key in ("L", "G", "M"):
        normalized = s.factor(
            old["MS_total_vertex_counterterms"][key] * 16 * s.pi**2 / Ibar
        )
        checks["same_parent_total_pole_" + key] = s.factor(
            normalized.subs(replacements) - Q * residues[key].subs(N, 6)
        )
    old_y = proper.data()["total_Yukawa_UV_counterterm_over_y"]
    symbol_map = {"Y": Y, "a": a, "Cf": Cf, "Q": Q}
    checks["same_parent_total_y_pole"] = s.factor(
        old_y.xreplace({v: symbol_map[str(v)] for v in old_y.free_symbols}) * y
        - residues["y"]
    )
    for key, w in weights.items():
        checks.update(
            {
                key + "_" + name: value
                for name, value in field_series(w)["checks"].items()
            }
        )
    checks.update(
        {
            "fundamental_G_square_cross_term_retained": s.expand(
                g2 - G**2 * (3 * k0**2 - 2 * t0 - k1 * L / Q)
            ),
            "fundamental_y_square_cross_term_retained": s.expand(
                Y2 - y**2 * (k0**2 - t0 - k1 * (Y - 4 * a * Cf) / Q)
            ),
            "forward_pole_commutator_relative": s.factor(
                comm_g * s.diff(tree, g) / tree + k1 * L / Q
            ),
            "quartic_commutator_is_constant_in_forward_invariant": s.diff(
                -2 * k1 * residues["L"], s.Symbol("nu"), 2
            ),
            "heavy_mass_unchanged_by_Phi_field_only": second["M"],
            "inert_yukawa_remains_zero": second["y"].subs({y: 0, Y: 0}),
        }
    )
    return {
        "Phi_field_weights": weights,
        "actual_one_loop_total_vertex_pole_residues": residues,
        "first_finite_coefficient_map": first,
        "second_finite_coefficient_map": second,
        "cubic_squared_second_coefficient": g2,
        "Yukawa_squared_second_coefficient": Y2,
        "local_MS_commutator_tree_b2_relative": -L * k1 / Q,
        "checks": checks,
        "scope": "The map is the finite part of the regulated total coefficients after multiplication by K_D^-w, at fixed MS input coordinates. It contains the undetermined new second slope t0 symbolically. Its finite pole commutator is a coordinate representation term, not a second addition to the fully renormalized field-rescaled amplitude.",
    }
