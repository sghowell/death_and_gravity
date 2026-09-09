"""Actual affine sampler transport, scalar reference and outgoing index form."""

from functools import cache

import sympy as sp
from p8a_double_null import curved


@cache
def data():
    a = sp.Symbol("a", positive=True)
    a1, a2, a3, a4 = sp.symbols("a1 a2 a3 a4", real=True)
    v, v1, v2, j, j1 = sp.symbols("v v1 v2 j j1", real=True)
    H = a1
    Hd = a * a2
    Hdd = a * a1 * a2 + a * a * a3
    Hthird = a * a1 * a1 * a2 + a * a * a2 * a2 + 3 * a * a * a1 * a3 + a**3 * a4
    parent = curved.data()
    symbols = {
        s.name: s
        for value in [
            parent["proper_second_time_operator"],
            parent["proper_mixed_operator"],
            parent["proper_state_weight_operator"],
            parent["physical_reference_null_numerator"],
        ]
        for s in value.free_symbols
    }
    values = {
        "a": a,
        "H": H,
        "Hdot": Hd,
        "Hddot": Hdd,
        "Hthird": Hthird,
        "f": v * j,
        "ft": a * v1 * j,
        "ftt": (a * a * v2 + a * a1 * v1) * j,
        "fz": v * j1,
        "ftz": a * v1 * j1,
    }
    substitution = {symbols[k]: val for k, val in values.items() if k in symbols}
    L = sp.factor(
        parent["proper_second_time_operator"].subs(substitution, simultaneous=True)
    )
    M = sp.factor(parent["proper_mixed_operator"].subs(substitution, simultaneous=True))
    G = sp.factor(
        parent["proper_state_weight_operator"].subs(substitution, simultaneous=True)
    )
    P = a * v2 - 2 * a1 * v1 + 2 * (a1 * a1 / a - a2) * v
    U = a * v1 - 2 * a1 * v
    V = v1 - 2 * a1 * v / a
    beta = sp.Symbol("beta_S", real=True)
    reference = sp.factor(
        parent["physical_reference_null_numerator"].subs(
            substitution, simultaneous=True
        )
    )
    reference_wanted = -4 * a1 * a1 * a2 / a + beta * (
        48 * a1 * a1 * a2 / a + 84 * a2 * a2 + 72 * a1 * a3 + 12 * a * a4
    )
    radius = sp.Symbol("coordinate_outgoing_radius", positive=True)
    b = a * radius
    b1 = a1 * radius + 1 / a
    b2 = a2 * radius
    h, h1 = sp.symbols("h h1", real=True)
    k, k1, k2 = sp.symbols("b b1 b2", real=True)
    boundary_derivative = (k2 / k - k1 * k1 / k**2) * h * h + 2 * k1 / k * h * h1
    return {
        "P": P,
        "Mraw": U,
        "V": V,
        "reference_null_numerator": reference,
        "outgoing_screen_scale": b,
        "outgoing_screen_first_derivative": b1,
        "outgoing_screen_second_derivative": b2,
        "checks": {
            "actual_affine_second_sampler_operator": sp.factor(L / a - P * j),
            "actual_affine_mixed_sampler_operator": sp.factor(M - U * j1),
            "actual_affine_state_sampler_operator": sp.factor(
                G / a - V * j - v * j1 / a**2
            ),
            "actual_affine_scalar_reference_all_four_jets": sp.factor(
                reference - reference_wanted
            ),
            "proper_to_affine_second_Hubble_derivative": sp.expand(
                Hdd - a * (a1 * a2 + a * a3)
            ),
            "proper_to_affine_third_Hubble_derivative": sp.expand(
                Hthird
                - a
                * (
                    sp.diff(Hdd, a) * a1
                    + sp.diff(Hdd, a1) * a2
                    + sp.diff(Hdd, a2) * a3
                    + sp.diff(Hdd, a3) * a4
                )
            ),
            "outgoing_screen_first_derivative": sp.factor(
                b1 - sp.diff(b, a) * a1 - sp.diff(b, radius) / a**2
            ),
            "outgoing_screen_second_derivative_cancels_coordinate_radius_terms": sp.factor(
                b2
                - sp.diff(b1, a) * a1
                - sp.diff(b1, a1) * a2
                - sp.diff(b1, radius) / a**2
            ),
            "actual_outgoing_Jacobi_equation": sp.factor(b2 / b - a2 / a),
            "exact_null_index_square_identity": sp.expand(
                h1 * h1 + k2 / k * h * h - (h1 - k1 / k * h) ** 2 - boundary_derivative
            ),
            "null_Ricci_affine_change": sp.factor(2 * Hd / a**2 - 2 * a2 / a),
            "physical_plane_and_null_affine_measures_match": sp.factor(a / a - 1),
        },
    }
