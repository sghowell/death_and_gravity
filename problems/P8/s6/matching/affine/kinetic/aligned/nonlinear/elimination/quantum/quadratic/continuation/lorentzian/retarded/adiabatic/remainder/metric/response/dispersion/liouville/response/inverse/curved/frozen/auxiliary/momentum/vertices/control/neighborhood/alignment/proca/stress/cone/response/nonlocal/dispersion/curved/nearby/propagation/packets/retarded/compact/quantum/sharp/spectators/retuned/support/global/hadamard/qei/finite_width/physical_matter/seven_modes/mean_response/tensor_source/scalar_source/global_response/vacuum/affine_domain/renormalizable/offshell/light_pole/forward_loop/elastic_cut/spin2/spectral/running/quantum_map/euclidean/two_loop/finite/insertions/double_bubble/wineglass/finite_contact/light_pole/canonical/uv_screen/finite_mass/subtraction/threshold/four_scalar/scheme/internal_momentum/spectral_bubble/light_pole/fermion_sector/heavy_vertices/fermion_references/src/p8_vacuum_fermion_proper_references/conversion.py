"""Proper UV counterterms and finite zero-momentum canonical dictionary."""

from functools import cache

import sympy as sp
from p8_vacuum_gauge_yukawa_screen import flow


@cache
def data():
    Y, a, Cf, N, Q = sp.symbols("Y a Cf N Q", positive=True)
    J0, J1, R, w = sp.symbols("J0 J1 R scalar_Phi_finite_kinetic", real=True)
    h = sp.Symbol("loop_marker", real=True)
    z = (-Y * J1 + a * Cf / 2) / Q
    eta = (Y * J0 + 2 * a * Cf) / Q
    ups = (Y * (J0 + 2 * R) - 6 * a * Cf) / Q
    psi_uv = -(Y / 2 + a * Cf) / Q
    phi_uv = -2 * N * Y / Q
    mass_total = (Y - 4 * a * Cf) / Q
    vertex_total = (Y - 4 * a * Cf) / Q
    bare_mass = mass_total - psi_uv
    bare_y = vertex_total - psi_uv - phi_uv / 2
    mass_ratio = (1 + h * eta) / (1 + h * z)
    yukawa_ratio = (1 + h * ups) / ((1 + h * z) * sp.sqrt(1 + h * w))
    y = sp.Symbol("y", positive=True)
    old = flow.data()["beta_Yukawa_times_loop_denominator"]
    old = old.xreplace(
        {s: a if str(s) == "positive_gauge_squared" else y for s in old.free_symbols}
    )
    checks = {
        "mass_and_kinetic_UV_poles_cancel": mass_total + (-Y + 4 * a * Cf) / Q,
        "fermion_kinetic_UV_pole_cancels": psi_uv + (Y / 2 + a * Cf) / Q,
        "proper_local_Yukawa_UV_pole_cancels": vertex_total + (-Y + 4 * a * Cf) / Q,
        "bare_mass_counterterm_includes_fermion_field": sp.factor(
            bare_mass - (3 * Y / 2 - 3 * a * Cf) / Q
        ),
        "bare_Yukawa_counterterm_includes_both_fields": sp.factor(
            bare_y - ((N + sp.Rational(3, 2)) * Y - 3 * a * Cf) / Q
        ),
        "actual_Yukawa_beta_matches_independent_frozen_matrix_result": sp.factor(
            (2 * Q * bare_y * y).subs({N: 6, Cf: sp.Rational(4, 3), Y: y * y}) - old
        ),
        "finite_mass_conversion_first_order": sp.diff(mass_ratio, h).subs(h, 0)
        - (eta - z),
        "finite_Yukawa_conversion_first_order": sp.diff(yukawa_ratio, h).subs(h, 0)
        - (ups - z - w / 2),
        "finite_mass_reference_is_canonical": sp.factor(
            mass_ratio * (1 + h * z) - (1 + h * eta)
        ),
        "finite_Yukawa_both_external_fermions_and_one_Phi": sp.simplify(
            yukawa_ratio * (1 + h * z) * sp.sqrt(1 + h * w) - (1 + h * ups)
        ),
        "inert_flavor_has_no_scalar_kinetic_correction": z.subs(Y, 0)
        - a * Cf / (2 * Q),
        "inert_flavor_retains_gauge_mass_correction": eta.subs(Y, 0) - 2 * a * Cf / Q,
        "inert_flavor_scalar_Yukawa_vertex_is_zero": (y * ups)
        .subs(Y, y * y)
        .subs(y, 0),
        "scalar_Yukawa_tadpoles_cancel_opposite_flavors": y + (-y),
    }
    return {
        "finite_fermion_kinetic_z": z,
        "finite_mass_relative_eta": eta,
        "finite_Yukawa_relative_upsilon": ups,
        "fermion_total_kinetic_UV_counterterm": psi_uv,
        "fermion_total_mass_UV_counterterm_over_m": mass_total,
        "total_Yukawa_UV_counterterm_over_y": vertex_total,
        "bare_mass_UV_counterterm_over_m": bare_mass,
        "bare_Yukawa_UV_counterterm_over_y": bare_y,
        "selected_zero_momentum_mass_ratio": mass_ratio,
        "selected_zero_momentum_Yukawa_ratio": yukawa_ratio,
        "first_order_mass_conversion": eta - z,
        "first_order_Yukawa_conversion": ups - z - w / 2,
        "scalar_field_dictionary": "w=r_scalar-fp from the complete one-loop physical Phi pole; not a new zero-momentum scalar field normalization",
        "counterterm_ownership": "These proper fermion kinetic, mass and scalar-Yukawa counterterms enter -Tr[D_F^-1 D_F^(1)] in S6.137. Any occurrence assigned to a renormalized proper subgraph must be removed from a separate insertion list.",
        "scope": "The ratios are the selected one-loop local functional's finite dictionary; their formal expansion is used for matching. MS and zero-momentum canonical couplings are not declared identical. No all-order resummation, charged-fermion pole, or complete two-loop amplitude is inferred.",
        "checks": checks,
    }
