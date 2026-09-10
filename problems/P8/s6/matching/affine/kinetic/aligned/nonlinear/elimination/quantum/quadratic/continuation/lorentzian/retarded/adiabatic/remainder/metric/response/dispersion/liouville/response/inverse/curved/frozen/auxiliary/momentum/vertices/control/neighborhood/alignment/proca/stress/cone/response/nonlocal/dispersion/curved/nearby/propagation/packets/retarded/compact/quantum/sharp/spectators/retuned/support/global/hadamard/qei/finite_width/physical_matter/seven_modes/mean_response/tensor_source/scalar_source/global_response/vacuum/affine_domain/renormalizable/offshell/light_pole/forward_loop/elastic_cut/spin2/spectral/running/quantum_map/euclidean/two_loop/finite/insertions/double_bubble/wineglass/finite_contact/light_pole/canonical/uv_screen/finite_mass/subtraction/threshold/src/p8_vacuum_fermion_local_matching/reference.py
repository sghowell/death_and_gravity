"""Explicit vacuum, mass and canonical-field reference ledger."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model as scalar_parent


@cache
def data():
    s, h = sp.symbols("channel_s loop_order")
    f0, f1, fp = sp.symbols("fermion_f_zero fermion_f_one fermion_slope_one", real=True)
    L, G, M, y, v4, V0 = sp.symbols(
        "quartic_reference cubic_reference heavy_mass_squared Yukawa_reference fermion_quartic_threshold fermion_vacuum_constant",
        real=True,
    )
    phi, H = sp.symbols("canonical_Phi heavy_H", real=True)
    f = sp.Function("fermion_kernel")
    kappa = 1 - fp
    fR = f(s) - f1 - (s - 1) * fp
    kernel = (1 - f1 - s + f(s)) / kappa
    canonical_mass = (1 - f1 + f0) / kappa
    Gcan = G / kappa
    Lcan = (L + v4) / (kappa * kappa)
    d = L - 3 * G * G / M
    dcan = (d + v4) / (kappa * kappa)
    potential = (
        canonical_mass * phi * phi / 2
        + Lcan * phi**4 / 24
        + Gcan * H * phi * phi / 2
        + M * H * H / 2
    )
    square = M * (H + Gcan * phi * phi / (2 * M)) ** 2 / 2
    old_potential = scalar_parent.data()["canonical_potential"]
    old_symbols = {str(symbol): symbol for symbol in old_potential.free_symbols}
    old_cubic = sp.diff(
        old_potential,
        old_symbols["canonical_heavy_field"],
        old_symbols["canonical_light_field"],
        2,
    )
    checks = {
        "unchanged_parent_heavy_cubic_sign": sp.factor(
            sp.diff(potential, H, phi, 2) / Gcan
            - old_cubic / old_symbols["positive_cubic_coupling"]
        ),
        "canonical_kernel_exact_reference_identity": sp.factor(
            kernel - (1 - s) - fR / kappa
        ),
        "canonical_mass_remainder_identity": sp.factor(
            canonical_mass - 1 - (f0 - f1 + fp) / kappa
        ),
        "canonical_quartic_after_heavy_elimination": sp.factor(
            Lcan - 3 * Gcan * Gcan / M - dcan
        ),
        "full_local_potential_completed_square": sp.factor(
            potential - square - canonical_mass * phi * phi / 2 - dcan * phi**4 / 24
        ),
        "all_flavor_vacuum_constant_counterterm_cancellation": -V0 + V0,
        "first_order_cubic_field_dictionary": sp.expand(
            sp.series(G / (1 - h * fp), h, 0, 2).removeO() - G - h * fp * G
        ),
        "first_order_quartic_field_and_threshold_dictionary": sp.expand(
            sp.series((L + h * v4) / (1 - h * fp) ** 2, h, 0, 2).removeO()
            - L
            - h * (v4 + 2 * fp * L)
        ),
        "first_order_Yukawa_field_dictionary": sp.expand(
            sp.series(y / sp.sqrt(1 - h * fp), h, 0, 2).removeO() - y - h * fp * y / 2
        ),
    }
    return {
        "vacuum_energy_counterterm_increment": -V0,
        "mass_reference_parameter": 1 - f1,
        "kinetic_normalization": kappa,
        "canonical_field_map": "Phi_can=sqrt(kappa) Phi_reference; kappa=1-f'(1)>0",
        "canonical_two_point_kernel": kernel,
        "canonical_potential_mass_squared": canonical_mass,
        "canonical_heavy_cubic": Gcan,
        "canonical_Yukawa": y / sp.sqrt(kappa),
        "canonical_zero_momentum_potential_quartic": Lcan,
        "canonical_heavy_eliminated_quartic": dcan,
        "completed_heavy_square": square,
        "scope": "These are exact constant-field rescalings of the selected tree-scalar plus one-loop-fermion functional. Their first-order expansions are separately checked. They are not an all-loop resummation theorem or a claim that canonical couplings remain equal to the prospective MSbar fixed-flow inputs. The old scalar-loop reference and bounce state are not silently changed.",
        "checks": checks,
    }
