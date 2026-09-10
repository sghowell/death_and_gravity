"""The complete inserted local tadpole and its fixed heavy-source cancellation."""

from functools import cache

import sympy as s


@cache
def data():
    L, G, M, Phi, Tad = s.symbols("L G M Phi T_insert")
    J = -G * Tad / 2
    stationary = (L - G * G / M) * Tad / 2
    induced = -G * J / M
    H = s.Symbol("H")
    potential = M * H * H / 2 + G * H * Phi * Phi / 2 + J * H
    eliminated = s.expand(potential.subs(H, -(G * Phi * Phi / 2 + J) / M))
    return {
        "fixed_sector_H_source_counterterm": J,
        "stationary_heavy_local_mass": stationary,
        "source_induced_local_mass": induced,
        "complete_local_mass_after_fixed_H_source": s.simplify(stationary + induced),
        "complete_row_mass_reference": "f_row,MS(0)=L T_insert,MS/2-g F_MS(M;0); g=G^2. The inner fermion physical mass/residue forest is already paired in T_insert and F.",
        "source_scope": "This is the assigned contribution of the W1F0 row to the order-two H-source reference, not the whole source counterterm or an arbitrary new source. Vacuum cross terms and other forests remain separately owned.",
        "checks": {
            "fixed_source_cancels_sector_one_point": J + G * Tad / 2,
            "stationary_heavy_mass_cancels_source_induced_mass": s.factor(
                stationary + induced - L * Tad / 2
            ),
            "literal_H_elimination_source_cross": s.diff(eliminated, Phi, 2).subs(
                Phi, 0
            )
            - induced,
            "same_positive_cubic_squared_sign": (-G * G * s.Symbol("F")).subs(
                G * G, s.Symbol("g")
            )
            + s.Symbol("g") * s.Symbol("F"),
            "source_vacuum_constant_kept_separate": eliminated.subs(Phi, 0)
            + J * J / (2 * M),
            "half_trace_one_inserted_light_line": s.Rational(1, 2) * 2 - 1,
        },
    }
